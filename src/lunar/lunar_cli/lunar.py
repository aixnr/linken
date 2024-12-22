## -----------------------------------------------------------------------------
# Import modules
import sys
import argparse
from .doctor import lunar_doctor
from .index_manager import read_linken_contrib, pretty_print_list, pull_create_index


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
# Argparse CLI frontend
def cli():
    """
    Lunar has the following subcommands:
      lunar list
        - List all available references from aixnr/linken-contrib
      lunar doctor
        - Doctor checks for environment and for tools
      lunar create --ref <ref> [--local-ref <local reference fasta>]
        - Create reference from <ref>
        - Pulls data from aixnr/linken-contrib

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

    # Generate the "create" subcommand parser
    subparser_create = subparsers.add_parser("create", help="For generating index files of the reference genome")
    subparser_create.add_argument("--ref", type=str, required=False,
                                  help="Reference genome to use")

    subparser_create.add_argument("--local-ref", type=str, required=False,
                                  help="User-supplied reference genome to use")

    # Complete the activation of ArgumentParser()
    args = parser.parse_args()

    # Conditional switching
    if args.subcommand == "list":
        print("")
        read_linken_contrib(timestamp, source_data)
        pretty_print_list(source_data)
        print("")
    elif args.subcommand == "doctor":
        print("")
        lunar_doctor()
        print("")
    elif args.subcommand == "create":
        print("")
        read_linken_contrib(timestamp, source_data)
        if args.ref:
            pull_create_index(source_data, ref_id=args.ref)
        elif args.local_ref:
            pull_create_index(source_data, ref_path=args.local_ref)
        else:
            print("  No valid reference supplied...")
            print("  Exiting...")
            sys.exit()
        print("")
    else:
        print(banner_lunar)
        parser.print_help()
