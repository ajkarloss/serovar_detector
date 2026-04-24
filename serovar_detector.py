#!/usr/bin/env python3

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def validate_snakemake(tool_dir, debug=False):
    """
    Validate that workflow/Snakefile exists relative to this script.
    This makes the script work from any current working directory.
    """
    workflow_dir = tool_dir / "workflow"
    snakefile = workflow_dir / "Snakefile"

    if not workflow_dir.is_dir():
        eprint(f"ERROR: No workflow directory detected at: {workflow_dir}")
        eprint("Software is probably corrupt. Consider redownloading.")
        sys.exit(1)

    if not snakefile.is_file():
        eprint(f"ERROR: No Snakefile detected at: {snakefile}")
        eprint("Software is probably corrupt. Consider redownloading.")
        sys.exit(1)

    if debug:
        print(f"Snakefile detected: {snakefile}")

    return snakefile


def ensure_writable_dir(path):
    """
    Create a directory if needed and verify it is writable.
    """
    path.mkdir(parents=True, exist_ok=True)
    test_file = path / ".write_test"
    try:
        with open(test_file, "w") as handle:
            handle.write("ok\n")
        test_file.unlink()
    except OSError as exc:
        raise OSError(f"Directory is not writable: {path}") from exc


def resolve_conda_prefix(args_conda_prefix, tool_dir, outdir, debug=False):
    """
    Resolve conda cache location in a portable way.

    Priority:
      1. --conda-prefix
      2. SEROVAR_DETECTOR_CONDA_PREFIX environment variable
      3. <tool_dir>/.snakemake/conda   (if writable)
      4. ~/.cache/serovar_detector/snakemake/conda   (if writable)
      5. <outdir>/.snakemake/conda

    This gives a reusable cache where possible, while still working in
    read-only installations and shared environments.
    """
    env_conda_prefix = os.environ.get("SEROVAR_DETECTOR_CONDA_PREFIX")

    candidates = []

    if args_conda_prefix:
        candidates.append(Path(args_conda_prefix).expanduser().resolve())
    elif env_conda_prefix:
        candidates.append(Path(env_conda_prefix).expanduser().resolve())
    else:
        candidates.append((tool_dir / ".snakemake" / "conda").resolve())

        home_dir = Path.home()
        if str(home_dir) not in ("", "."):
            candidates.append((home_dir / ".cache" / "serovar_detector" / "snakemake" / "conda").resolve())

        candidates.append((Path(outdir).resolve() / ".snakemake" / "conda").resolve())

    last_error = None
    for candidate in candidates:
        try:
            ensure_writable_dir(candidate)
            if debug:
                print(f"Using conda prefix: {candidate}")
            return str(candidate)
        except OSError as exc:
            last_error = exc
            if debug:
                print(f"Conda prefix not usable: {candidate} ({exc})")

    eprint("ERROR: Could not create a usable conda prefix directory.")
    if last_error:
        eprint(str(last_error))
    sys.exit(1)


def generate_configfile(database, outdir, threshold, append_results, threads, debug, tmpdir):
    """
    Write config into the job outdir so Snakemake (running with --directory outdir)
    reads config/config.yaml from that job-specific working directory.
    """
    outdir_path = Path(outdir).resolve()
    config_file = outdir_path / "config" / "config.yaml"
    config_file.parent.mkdir(parents=True, exist_ok=True)

    config = {
        "database": str(Path(database).resolve()),
        "outdir": str(outdir_path).rstrip("/"),
        "threshold": threshold,
        "append_results": append_results,
        "threads": threads,
        "debug": debug,
        "tmpdir": str(Path(tmpdir).resolve()),
    }

    with open(config_file, "w") as config_yaml:
        yaml.safe_dump(config, config_yaml)

    return str(config_file)


def screen_files(directory, file_type):
    """
    Scan input directories for reads or assemblies and return a metadata dataframe.
    """
    if not directory:
        if file_type == "Reads":
            return pd.DataFrame(columns=["sample_name", "mate", "file", "type"])
        if file_type == "Assembly":
            return pd.DataFrame(columns=["sample_name", "file", "type"])
        return pd.DataFrame(columns=["sample_name", "mate", "file", "type"])

    directory = str(Path(directory).resolve())

    if file_type == "Reads":
        fastqs = glob.glob(f"{directory}/*.fastq*")
        fqs = glob.glob(f"{directory}/*.fq*")
        read_files = sorted(set(fastqs + fqs))

        if len(read_files) == 0:
            return pd.DataFrame(columns=["sample_name", "mate", "file", "type"])

        pattern = (
            r"^(?P<file>\S+\/(?P<sample_name>\S+?)"
            r"((_S\d+)(_L\d+))?_(?P<mate>[Rr]?[12])(_\d{3})?"
            r"\.(?P<ext>((fastq)?(fq)?)?(\.gz)?))$"
        )
        search = [re.search(pattern, read_file) for read_file in read_files]
        search = [s for s in search if s is not None]

        if len(search) == 0:
            return pd.DataFrame(columns=["sample_name", "mate", "file", "type"])

        metadata_raw = pd.DataFrame(
            {
                "sample_name": [s.group("sample_name") for s in search],
                "mate": [s.group("mate") for s in search],
                "file": [s.group("file") for s in search],
                "type": file_type,
            }
        )

        metadata = metadata_raw.sort_values(by=["sample_name", "mate"]).reset_index(drop=True)

    elif file_type == "Assembly":
        fastas = glob.glob(f"{directory}/*.fasta")
        fas = glob.glob(f"{directory}/*.fa")
        assembly_files = sorted(set(fastas + fas))

        if len(assembly_files) == 0:
            return pd.DataFrame(columns=["sample_name", "file", "type"])

        pattern = r"^(?P<file>\S+\/(?P<sample_name>\S+)\.(fasta|fa))$"
        search = [re.search(pattern, assembly_file) for assembly_file in assembly_files]
        search = [s for s in search if s is not None]

        if len(search) == 0:
            return pd.DataFrame(columns=["sample_name", "file", "type"])

        metadata_raw = pd.DataFrame(
            {
                "sample_name": [s.group("sample_name") for s in search],
                "file": [s.group("file") for s in search],
                "type": file_type,
            }
        )

        metadata = metadata_raw.sort_values(by="sample_name").reset_index(drop=True)

    else:
        return pd.DataFrame(columns=["sample_name", "mate", "file", "type"])

    return metadata


def create_symlinks(metadata, outdir):
    """
    Create normalized symlinks into the temporary working directory.
    """
    outdir_path = Path(outdir).resolve()

    for _, sample_metadata in metadata.iterrows():
        sample_name = sample_metadata["sample_name"]

        if sample_metadata["type"] == "Reads":
            reads_dir = outdir_path / "reads"
            reads_dir.mkdir(parents=True, exist_ok=True)

            mate = sample_metadata["mate"]
            sample_file = Path(sample_metadata["file"]).resolve()
            sample_link = reads_dir / f"{sample_name}_{mate}.fastq.gz"
        else:
            assemblies_dir = outdir_path / "assemblies"
            assemblies_dir.mkdir(parents=True, exist_ok=True)

            sample_file = Path(sample_metadata["file"]).resolve()
            sample_link = assemblies_dir / f"{sample_name}.fasta"

        if sample_link.exists() or sample_link.is_symlink():
            sample_link.unlink()

        os.symlink(str(sample_file), str(sample_link))

    print("Files successfully linked!")
    return True


def make_sample_sheet(table):
    print("Generating sample sheet", end="... ")

    if len(table.index) > 0:
        sample_subset = table[["sample_name", "type"]]
        sample_sheet = sample_subset.drop_duplicates().reset_index(drop=True)
    else:
        eprint("ERROR: Sample sheet has no rows.")
        sys.exit(1)

    print(f"Success: A total of {len(sample_sheet.index)} samples have been annotated!")
    return sample_sheet


def write_sample_sheet(sample_sheet, pepdir):
    sample_file = Path(pepdir) / "sample_sheet.csv"
    sample_exists = sample_file.exists()

    print("Writing sample sheet", end="... ")
    sample_sheet.to_csv(sample_file, index=False)

    if sample_exists:
        print(f"Success: Overwriting {sample_file}")
    else:
        print(f"Success: Written to {sample_file}")


def write_subsample_sheet(subsample_sheet, pepdir):
    subsample_file = Path(pepdir) / "subsample_sheet.csv"
    subsample_exists = subsample_file.exists()

    print("Writing subsample sheet", end="... ")
    subsample_sheet.to_csv(subsample_file, index=False)

    if subsample_exists:
        print(f"Success: Overwriting {subsample_file}")
    else:
        print(f"Success: Written to {subsample_file}")


def write_pep(pepdir):
    pep_file = Path(pepdir) / "project_config.yaml"
    pep_exists = pep_file.exists()

    pep_text = (
        "pep_version: 2.1.0\n"
        "sample_table: 'sample_sheet.csv'\n"
        "subsample_table: 'subsample_sheet.csv'\n"
    )

    print("Writing project configuration file", end="... ")
    with open(pep_file, "w") as config_file:
        config_file.write(pep_text)

    if not pep_exists:
        print(f"Success: Written to {pep_file}")
    else:
        print(f"Success: Overwriting {pep_file}")


def generate_sheets(reads_dir, assembly_dir, enable_blacklist, outdir, blacklist_file, tmpdir):
    outdir_path = Path(outdir).resolve()
    tmpdir_path = Path(tmpdir).resolve()
    blacklist_path = Path(blacklist_file).resolve()

    outdir_path.mkdir(parents=True, exist_ok=True)
    tmpdir_path.mkdir(parents=True, exist_ok=True)

    reads_metadata = screen_files(directory=reads_dir, file_type="Reads")
    assembly_metadata = screen_files(directory=assembly_dir, file_type="Assembly")

    metadata = pd.concat([reads_metadata, assembly_metadata], ignore_index=True, sort=True)

    if blacklist_path.exists():
        print("Blacklist file detected, reading!")
        blacklist = pd.read_csv(blacklist_path, sep="\t")

        if "file" not in blacklist.columns:
            eprint(f"ERROR: Blacklist file exists but lacks required 'file' column: {blacklist_path}")
            sys.exit(1)

        exclude_samples = blacklist["file"].values
        metadata = metadata[~metadata["file"].isin(exclude_samples)].reset_index(drop=True)

        if not enable_blacklist:
            print("Warning: Blacklist was included, yet current samples will NOT be added to the blacklist file!")

    if len(metadata.index) > 0:
        create_symlinks(metadata, tmpdir_path)

        pepdir = outdir_path / "schemas"
        pepdir.mkdir(parents=True, exist_ok=True)

        sample_sheet = make_sample_sheet(metadata)
        write_sample_sheet(sample_sheet, pepdir)

        subsample_cols = [c for c in ["sample_name", "mate", "file"] if c in metadata.columns]
        subsample_sheet = metadata[subsample_cols]
        write_subsample_sheet(subsample_sheet, pepdir)

        sample_files = metadata["file"]
        write_pep(pepdir)
    else:
        print("No new samples detected after blacklist filtering.")
        sample_files = pd.Series(dtype=str)

    return sample_files


def update_blacklist(enable_blacklist, blacklist_file, blacklist_clean, sample_files):
    blacklist_path = Path(blacklist_file).resolve()
    blacklist_exists = blacklist_path.is_file()
    include_header = (not blacklist_exists) or blacklist_clean

    if enable_blacklist:
        mode = "w" if blacklist_clean else "a"
        if blacklist_exists and not blacklist_clean:
            print("Appending to existing blacklist file")
        if blacklist_exists and blacklist_clean:
            print("Overwriting preexisting blacklist file")
    elif blacklist_clean:
        mode = "w"
        if blacklist_exists:
            print("Overwriting preexisting blacklist file")
    else:
        return False

    sample_files.to_csv(blacklist_path, sep="\t", index=False, mode=mode, header=include_header)
    return True


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Screen read files and assemblies for serovar biomarker genes to provide "
            "suggestions for isolate serovar. Currently supports Actinobacillus "
            "pleuropneumoniae."
        )
    )
    parser.add_argument("-r", dest="reads_dir", help="Input path to reads directory", required=False)
    parser.add_argument("-a", dest="assembly_dir", help="Input path to assembly directory", required=False)
    parser.add_argument(
        "-D",
        dest="database",
        help="Path and prefix to kmer-aligner database",
        default="./db",
    )
    parser.add_argument(
        "-o",
        dest="outdir",
        help="Output path to results and temporary files directory",
        required=True,
    )
    parser.add_argument(
        "-T",
        dest="threshold",
        help="Cutoff threshold of match coverage and identity. Ignore threshold by setting to 0. Default: 98",
        default=98,
    )
    parser.add_argument(
        "-R",
        dest="append_results",
        help="Append to existing results file",
        action="store_true",
    )
    parser.add_argument(
        "-b",
        dest="enable_blacklist",
        help="Update existing blacklist file with new samples. Creates a blacklist file if none exists",
        action="store_true",
    )
    parser.add_argument(
        "-B",
        dest="blacklist_clean",
        help="Ignore and overwrite existing blacklist file. Creates a blacklist if none exists",
        action="store_true",
    )
    parser.add_argument(
        "-k",
        dest="keep_tmp",
        help="Preserve temporary files such as KMA result files",
        action="store_true",
    )
    parser.add_argument(
        "-t",
        dest="threads",
        help="Number of threads to allocate for the pipeline. Default: 3",
        default=3,
    )
    parser.add_argument(
        "-F",
        dest="force",
        help="Force rerun of all tasks in pipeline",
        action="store_true",
    )
    parser.add_argument(
        "-n",
        dest="dry_run",
        help="Perform a dry run with Snakemake to see jobs without executing them",
        action="store_true",
    )
    parser.add_argument(
        "-d",
        dest="debug",
        help="Enable debug mode",
        action="store_true",
    )
    parser.add_argument(
        "--conda-prefix",
        dest="conda_prefix",
        help=(
            "Path to Snakemake conda environment cache. "
            "Can also be set via SEROVAR_DETECTOR_CONDA_PREFIX. "
            "Default: auto-resolved portable cache location."
        ),
        required=False,
    )
    parser.add_argument(
        "--conda-frontend",
        dest="conda_frontend",
        choices=["conda", "mamba"],
        default=None,
        help="Optional Snakemake conda frontend to use if supported by installed Snakemake",
    )
    parser.add_argument(
        "--latency-wait",
        dest="latency_wait",
        type=int,
        default=240,
        help=(
            "Snakemake latency wait in seconds for shared/network filesystems. "
            "Default: 240 (tuned for HPC shared FS). Lower on local disks if desired."
        ),
    )

    return parser.parse_args()


def main():
    args = parse_arguments()

    tool_dir = Path(__file__).resolve().parent
    outdir = str(Path(args.outdir).resolve())
    tmpdir = f"{outdir}/tmp"
    blacklist_file = f"{outdir}/blacklist.tsv"

    conda_prefix = resolve_conda_prefix(args.conda_prefix, tool_dir, outdir, debug=args.debug)

    db_arg = args.database
    if os.path.isabs(db_arg):
        database = str(Path(db_arg).resolve())
    else:
        database = str((tool_dir / db_arg).resolve())

    try:
        threshold = int(args.threshold)
    except Exception:
        eprint(f"ERROR: Invalid threshold value: {args.threshold}")
        sys.exit(1)

    try:
        threads = int(args.threads)
        if threads < 1:
            raise ValueError
    except Exception:
        eprint(f"ERROR: Invalid thread count: {args.threads}")
        sys.exit(1)

    append_results = args.append_results
    enable_blacklist = args.enable_blacklist
    blacklist_clean = args.blacklist_clean
    keep_tmp = args.keep_tmp
    force = args.force
    dry_run = args.dry_run
    debug = args.debug

    reads_dir = str(Path(args.reads_dir).resolve()) if args.reads_dir else None
    assembly_dir = str(Path(args.assembly_dir).resolve()) if args.assembly_dir else None

    if not reads_dir and not assembly_dir:
        eprint("ERROR: At least one of -r/--reads_dir or -a/--assembly_dir must be provided.")
        sys.exit(1)

    snakefile = validate_snakemake(tool_dir, debug)

    if enable_blacklist and blacklist_clean and os.path.isfile(blacklist_file):
        eprint("Blacklist file detected, and BOTH blacklist update (-b) and blacklist clean (-B) were selected.")
        eprint("Please choose either update existing blacklist (-b) OR overwrite (-B), not both.")
        sys.exit(1)

    generate_configfile(
        database=database,
        outdir=outdir,
        threshold=threshold,
        append_results=append_results,
        threads=threads,
        debug=debug,
        tmpdir=tmpdir,
    )

    outdir_path = Path(outdir).resolve()
    config_dir = outdir_path / "config"
    config_dir.mkdir(parents=True, exist_ok=True)

    profiles_src = tool_dir / "config" / "serovar_profiles.yaml"
    profiles_dst = config_dir / "serovar_profiles.yaml"

    if not profiles_src.is_file():
        eprint(f"ERROR: Missing serovar profiles file: {profiles_src}")
        sys.exit(1)

    shutil.copy2(profiles_src, profiles_dst)

    if debug:
        print(f"Staged serovar profiles: {profiles_dst}")

    sample_files = generate_sheets(
        reads_dir=reads_dir,
        assembly_dir=assembly_dir,
        enable_blacklist=enable_blacklist,
        outdir=outdir,
        blacklist_file=blacklist_file,
        tmpdir=tmpdir,
    )

    if len(sample_files) == 0:
        print("Nothing to do, exiting!")
        sys.exit(0)

    snakemake_cmd = [
        "snakemake",
        "--snakefile",
        str(snakefile),
        "--directory",
        str(outdir_path),
        "--use-conda",
        "--conda-prefix",
        conda_prefix,
        "--cores",
        str(threads),
    ]

    if args.latency_wait < 0:
        eprint("ERROR: --latency-wait must be >= 0")
        sys.exit(1)
    snakemake_cmd.extend(["--latency-wait", str(args.latency_wait)])

    if args.conda_frontend:
        snakemake_cmd.extend(["--conda-frontend", args.conda_frontend])

    if force or blacklist_clean:
        snakemake_cmd.append("-F")
    if dry_run:
        snakemake_cmd.append("-n")

    if debug:
        print("Running command:")
        print(" ".join(snakemake_cmd))

    results_file = outdir_path / "serovar.tsv"
    do_append = append_results and results_file.is_file()

    if do_append:
        results_tmp = results_file.with_suffix(".tmp")
        print(f"Copying {results_file} to {results_tmp}")
        shutil.copy(results_file, results_tmp)

    if shutil.which("snakemake") is None:
        eprint("ERROR: 'snakemake' command not found in PATH.")
        sys.exit(1)

    result = subprocess.run(snakemake_cmd, check=False)
    if result.returncode != 0:
        eprint("Something went wrong while executing Snakemake.")
        sys.exit(result.returncode)

    update_blacklist(
        enable_blacklist=enable_blacklist,
        blacklist_file=blacklist_file,
        blacklist_clean=blacklist_clean,
        sample_files=sample_files,
    )

    if do_append:
        print("Appending new results to existing results")
        results_tmp = results_file.with_suffix(".tmp")
        serovar_new = pd.read_csv(results_file, sep="\t")
        shutil.move(results_tmp, results_file)
        serovar_new.to_csv(results_file, sep="\t", index=False, mode="a", header=False)

    if not keep_tmp:
        print("Cleaning up temporary files.")
        shutil.rmtree(tmpdir, ignore_errors=True)

    print("All Done!")


if __name__ == "__main__":
    main()