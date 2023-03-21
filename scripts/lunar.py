#!/usr/bin/env python
# Lunar module for interfacing with user


## -----------------------------------------------------------------------------
# Import modules
import argparse
import shutil
from urllib.request import urlopen
import json
import os
import pwd
import pathlib
import sys
import subprocess
import pprint
from typing import Optional


## -----------------------------------------------------------------------------
# Fancy banner
banner_lunar = r"""
      ___       ___           ___           ___           ___     
     /\__\     /\__\         /\__\         /\  \         /\  \    
    /:/  /    /:/  /        /::|  |       /::\  \       /::\  \   
   /:/  /    /:/  /        /:|:|  |      /:/\:\  \     /:/\:\  \  
  /:/  /    /:/  /  ___   /:/|:|  |__   /::\~\:\  \   /::\~\:\  \ 
 /:/__/    /:/__/  /\__\ /:/ |:| /\__\ /:/\:\ \:\__\ /:/\:\ \:\__\
 \:\  \    \:\  \ /:/  / \/__|:|/:/  / \/__\:\/:/  / \/_|::\/:/  /
  \:\  \    \:\  /:/  /      |:/:/  /       \::/  /     |:|::/  / 
   \:\  \    \:\/:/  /       |::/  /        /:/  /      |:|\/__/  
    \:\__\    \::/  /        /:/  /        /:/  /       |:|  |    
     \/__/     \/__/         \/__/         \/__/         \|__|    

  Linken's Lunar Utility Script
  https://github.com/aixnr/linken
  https://github.com/aixnr/linken-contrib
  """


## -----------------------------------------------------------------------------
# Global variables
timestamp: dict = {}
source_data: dict = {}


## -----------------------------------------------------------------------------
# Functions
def read_linken_contrib() -> None:
    """
    Fetch records
    """
    location_source = "https://raw.githubusercontent.com/aixnr/linken-contrib/main/source.json"
    _response = urlopen(location_source)
    _response = json.loads(_response.read())
    timestamp["data"] = _response["updated_on"]
    source_data["data"] = _response["refs"]


def pretty_print_list() -> None:
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


def lunar_doctor() -> None:
    """
    Lunar Doctor checks for environment
    """
    _tools = ["curl", "python", "nextflow", "bwa", "samtools", "bwa-mem2", "minimap2", "fastqc",
              "cutadapt", "trimgalore", "lofreq", "bcftools", "gatk", "deeptools"]

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


def pull_create_index(ref_id: str) -> None:
    """
    Pull reference genome fasta and generate index files.
    
    TODO: Handle stdout/stderr if returncode =/= 0 for subprocess calls.
    """
    _data_ref: list[dict] = source_data["data"]
    _data_id_key: dict[str, str] = {}

    for _item in _data_ref:
        _data_id_key[_item["id"]] = _item["location"]

    if ref_id not in _data_id_key.keys():
        print(f"  {ref_id} does not exist in our collection.")
        print("  Program exits...")
        sys.exit()

    _download_url = _data_id_key[ref_id]
    print(f"  {ref_id} selected...")
    print(f"  Downloading ref genome from {_download_url}")

    def generate_index():
        _index_dir = pathlib.Path("index")
        if _index_dir.exists():
            print("  Warning! The index directory already exists! I can't overwrite it")
            print("  Exiting...")
            print("")
            sys.exit()
        else:
            os.mkdir("index")

        proc_download = subprocess.run(["curl", "-L", _download_url, "--output", "./index/ref.fasta"], capture_output=True, text=True)
        if proc_download.returncode == 0:
            print("  Download completed!")
        else:
            print("  Download failed...")
            sys.exit()

        check_bwa_to_exit = shutil.which("bwa")
        if check_bwa_to_exit == None:
            print("  Error! The 'bwa' tool does not exist...")
            print("  Exiting...")
            print("")
            sys.exit()

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

    generate_index()
    print(f"  Index files for ref genome {ref_id} were successfully generated!")


def nextflow_list_scripts(stage: Optional[str] = None) -> None:
    with open("/nextflow/pipeline_list.json") as json_file:
        _data_nextflow: list[dict] = json.load(json_file)["pipelines"]

    _data_full: dict[str, dict] = {}
    for _, _item in enumerate(_data_nextflow):
        _data_full[_item["id"]] = {
            "Path": _item["path"],
            "Date Created": _item["dateCreated"],
            "Date Modified": _item["dateModified"],
            "Description": _item["description"],
            "Tasks": _item["tasks"]
        }
    
    if not stage:
        pprint.pprint(_data_full)

    else:
        if stage not in _data_full.keys():
            print(f"  [ERR] Pipeline {stage} not found...")
            print("  [ERR] Please run lunar flows to see available pipelines...")
            print("  [ERR] Exiting...")
            sys.exit()

        _chosen_pipeline = _data_full[stage]
        _pipeline_path = _chosen_pipeline['Path']

        _cwd = os.getcwd()

        _proc_copy_pipeline = subprocess.run(["cp", _pipeline_path, _cwd], capture_output=True, text=True)
        if _proc_copy_pipeline.returncode == 0:
            print(f"  [INFO] Copied pipeline {stage}!")
        else:
            print("  [INFO] Failed to copy pipeline...")
            sys.exit()


## -----------------------------------------------------------------------------
# Argparse CLI frontend
def cli():
    """
    Lunar has five subcommands:
      lunar list
        - List all available references from aixnr/linken-contrib
      lunar doctor
        - Doctor checks for environment and for tools
      lunar create --ref <ref> [--local-ref <local reference fasta>]
        - TODO: update the implementation
        - Create reference from <ref>
        - Pulls data from aixnr/linken-contrib
      lunar flows
        - List all available nextflow scripts. Nextflow script definitions
          are available at $PROJECTROOT/nextflow/pipeline_list.json
        - The whole $PROJECTROOT/nextflow directory is copied to /nextflow in the final container
      lunar stage --flow <nextflow script>
        - Lunar copies nextflow script from /nextflow/*.nf to $PWD
      lunar clean
        :TODO pending implementation
        - Clear current directory; removes work and .nextflow*

    By default, it creates a folder in $PWD called "index" where it downloads ref.fasta.
    Then, it calls bwa, samtools, and picard.jar to generate the index files in $PWD/index.
    """
    # Create top-level parser
    parser = argparse.ArgumentParser()

    # Initialize sub-parsers under "subcommand"
    subparsers = parser.add_subparsers(dest="subcommand")

    # Add parser for subcommands that do not accept additional arguments
    subparsers.add_parser("list", help="For listing all available reference genomes")
    subparsers.add_parser("doctor", help="For checking if associated tools exist in $PATH")
    subparsers.add_parser("flows", help="List all available nextflow scripts")

    # Generate the "stage" subcommand parser
    subparser_stage = subparsers.add_parser("stage", help="Stage the specified nextflow script in current directory")
    subparser_stage.add_argument("--flow", type=str, required=True,
                                 default="pipeline_001", help="Choose a pipeline from the available list")

    # Generate the "create" subcommand parser
    subparser_create = subparsers.add_parser("create", help="For generating index files of the reference genome")
    subparser_create.add_argument("--ref", type=str, required=True,
                                  default="IAV-PR8-WG", help="Reference genome to use")

    # Complete the activation of ArgumentParser()
    args = parser.parse_args()

    # Conditional switching
    if args.subcommand == "list":
        print("")
        read_linken_contrib()
        pretty_print_list()
        print("")
    elif args.subcommand == "doctor":
        print("")
        lunar_doctor()
        print("")
    elif args.subcommand == "create":
        print("")
        read_linken_contrib()
        pull_create_index(ref_id=args.ref)
        print("")
    elif args.subcommand == "flows":
        nextflow_list_scripts()
    elif args.subcommand == "stage":
        nextflow_list_scripts(stage=args.flow)
    else:
        print(banner_lunar)
        parser.print_help()


## -----------------------------------------------------------------------------
# main
if __name__ == "__main__":
    cli()

