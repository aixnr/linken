from urllib.request import urlopen
from pathlib import Path
import json
import sys
import subprocess
import shutil
import pathlib
import os


def read_linken_contrib(timestamp: dict, source_data:dict) -> None:
    """
    Fetch index from central repo at 'github.com/aixnr/linken-contrib' and print the list
    """
    location_source = "https://raw.githubusercontent.com/aixnr/linken-contrib/main/source.json"
    _response = urlopen(location_source)
    _response = json.loads(_response.read())
    timestamp["data"] = _response["updated_on"]
    source_data["data"] = _response["refs"]


def pretty_print_list(source_data: dict) -> None:
    """
    Pretty print list of reference genomes available for download.
    """
    _data_ref: list[dict] = source_data["data"]

    _data_print: dict[int, list] = {}
    for _index, _item in enumerate(_data_ref):
        _data_print[_index] = [_item["id"], _item["description"], _item["source"]]

    # Print data header
    print("  {:<15} {:<42} {:<10}".format("id", "description", "source"))
    print("  {:<15} {:<42} {:<10}".format("==========", "==========", "=========="))

    # Print data
    for _, _v in _data_print.items():
        print("  {:<15} {:<42} {:<10}".format(*_v))


def pull_create_index(source_data: dict, ref_id: str = None, ref_path: str = None) -> None:
    _data_ref: list[dict] = source_data["data"]
    _data_id_key: dict[str, str] = {}

    for _item in _data_ref:
        _data_id_key[_item["id"]] = _item["location"]

    if ref_id:
        if ref_id not in _data_id_key.keys():
            print(f"  {ref_id} does not exist in our collection.")
            print("  Program exits...")
            sys.exit()

        _download_url = _data_id_key[ref_id]
        print(f"  {ref_id} selected...")
        print(f"  Downloading ref genome from {_download_url}")
        generate_index(location_index=_download_url, retrieve="net")
    elif ref_path:
        generate_index(location_index=ref_path, retrieve="local")
    else:
        print("  Reference genome not supplied...")
        print("  Exiting...")
        sys.exit()

    print(f"  Index files for the reference genome were successfully generated!")


def generate_index(location_index: str, retrieve: str = "net"):
    _index_dir = pathlib.Path("index")
    if _index_dir.exists():
        print("  Warning! The index directory already exists! I can't overwrite it")
        print("  Exiting...")
        print("")
        sys.exit()
    else:
        os.mkdir("index")

    if retrieve == "net":
        reference_download(link=location_index)
    elif retrieve == "local":
        reference_move(file_path=location_index)

    proc_bwa = subprocess.run(["bwa", "index", "./index/ref.fasta", "-p", "./index/ref"], capture_output=True, text=True)
    if proc_bwa.returncode == 0:
        print("  'bwa index' completed!")
    else:
        print("  'bwa index' failed...")
        sys.exit()

    proc_samtools = subprocess.run(["samtools", "faidx", "./index/ref.fasta"], capture_output=True, text=True)
    if proc_samtools.returncode == 0:
        print("  'samtools faidx' completed!")
    else:
        print("  'samtools faidx' failed...")
        sys.exit()

    proc_picard = subprocess.run(["java", "-jar", "/usr/local/lib/picard.jar", "CreateSequenceDictionary",
                                    "R=./index/ref.fasta", "O=./index/ref.dict"], capture_output=True, text=True)
    if proc_picard.returncode == 0:
        print("  'java -jar picard.jar CreateSequenceDictionary' completed!")
    else:
        print("  'java -jar picard.jar CreateSequenceDictionary' failed...")
        sys.exit()


def reference_download(link: str):
    proc_download = subprocess.run(["curl", "-L", link, "--output", "./index/ref.fasta"], capture_output=True, text=True)
    if proc_download.returncode == 0:
        print("  Download completed!")
    else:
        print("  Download failed...")
        sys.exit()


def reference_move(file_path: str):
    if Path(file_path).exists():
        print(f"  Using {file_path} as the reference genome...")
        print("  Moving the reference genome to ./index/ref.fasta")
        shutil.move(file_path, "./index/ref.fasta")
    else:
        print(f"  File {file_path} does not exist!")
        print("  Exiting....")
        sys.exit()
