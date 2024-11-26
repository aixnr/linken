import shutil
import pwd
import pathlib
import os

def lunar_doctor() -> None:
    """
    Lunar Doctor checks for environment
    """
    _tools = ["curl", "python", "bwa", "samtools", "bwa-mem2", "minimap2", "fastqc",
              "cutadapt", "trimgalore", "lofreq", "bcftools", "gatk", "muscle", "ivar",
              "eevee", "canu"]

    _tools_dict: dict[int, list] = {}
    for _i, _t in enumerate(_tools):
        _location = shutil.which(_t)
        _tools_dict[_i] = [_t, str(_location)]

    # Print current environment
    print(f"  Current user: {pwd.getpwuid(os.getuid()).pw_name}")
    print(f"  Current working directory ($PWD): {os.getcwd()}")
    print("")

    # Check if raw_reads directory is available
    _raw_reads_dir = pathlib.Path("raw_reads")
    if _raw_reads_dir.exists():
        print("  The raw_reads directory exists at this location. Please proceed...")
        print("")
    else:
        print("  Warning! The raw_reads directory is not present at this current location.")
        print("")

    # Check if index directory is available
    _index_dir = pathlib.Path("index")
    if _index_dir.exists():
        print("  The index directory exists at this location.")
        print("  Were the index files for a reference genome already generated?")
        print("")
    else:
        print("  The index directory not detected. Please generate index files before mapping.")
        print("")

    # Print header
    print("  {:<15} {:<15}".format("tool", "path"))
    print("  {:<15} {:<15}".format("========", "========"))

    for _, _v in _tools_dict.items():
        print("  {:<15} {:<15}".format(*_v))
