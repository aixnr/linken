#!/usr/bin/env python
# Eevee module for working on analyzing data and producing interpretable outputs


## -----------------------------------------------------------------------------
# Import modules
import argparse
from coverage import Coverage
from substitution import Substitution
from .eevee_bulk import config_read, config_doctor, substitution_bulk


## -----------------------------------------------------------------------------
# Fancy banner
banner_eevee = r"""
   ╔N
   █▓▌      Linken's Eevee Utility Script
   ╠▓▓▓     https://github.com/aixnr/linken
    █▓█▌  , ,      .,╓▄▄▄@`
     ▀█▒▒▒▒▒▒▒▒,,▄██▓▓▓▓╜ 
      ║░▒▒▒▒▒▒▒▒█████▀"  .«░░ 
      ▒██▒▒▒▒██▒╫  ╓╗╢╢╢╣@╣░░ 
     »▒▒▒▒░▒▒██╢, ╢╢╢╢╢╢╢╢▌╓░ 
   -  "╚▒▒▒▒▒╢▒▒▒░▓▓╢╢╢╢╢╢▓▓` 
   *░   ░░"░░░░░░░▒╢▓▓▓▓▓▓▓'
     "░░:░  ▒░░░▒▒╣╣▀▀▀▀" 
       ▒▒░▒░▒▒▒▒▓╣╢╢Γ  
       ╘╢▒▒╢╢▒▓▓ ▒▒▒
       ,▒▒▒▒▒╩▀ ó▒▒" 
          ╙╙" 
"""


## -----------------------------------------------------------------------------
# main() function
def main():
    # initialize top-level parser and 'command' subparser
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand")

    # Add parser for subcommand 'coverage'
    subparser_coverage = subparsers.add_parser("coverage",
                                               help="Provides analysis for coverage & variants.")
    subparser_coverage.add_argument("--sample", type=str, required=True,
                                    help="Sample name.")
    subparser_coverage.add_argument("--cov_tsv", type=str, required=True,
                                    help="Path to coverage (tsv from bam).")
    subparser_coverage.add_argument("--vcf", type=str, required=False, default=False,
                                    help="Path to variant table (tsv from vcf).")
    
    # Add parser for subcommand 'substitute'
    subparser_substitute = subparsers.add_parser("substitute",
                                                 help="Convert nucleotide variant to residue change.")
    subparser_substitute.add_argument("--sample", type=str, required=True,
                                      help="Sample name.")
    subparser_substitute.add_argument("--gene", type=str, required=True,
                                      help="Gene ID to map the variant as gene:start-stop; must match the header of the refseq.")
    subparser_substitute.add_argument("--ref", type=str, required=True,
                                      help="Path to the reference fasta.")
    subparser_substitute.add_argument("--vcf", type=str, required=True,
                                      help="Path to variant table (tsv from vcf).")
    subparser_substitute.add_argument("--csv_dir", type=str, required=True,
                                      help="Output directory of the final CSV.")
    
    # Add parser for subcommand 'substitute-bulk'
    subparser_sbulk = subparsers.add_parser("substitute-bulk", help="Running 'eevee substitute' in bulk.")
    subparser_sbulk.add_argument("--check", type=str, required=False,
                                 help="Check if the TOML file either valid or not.")
    subparser_sbulk.add_argument("--run", type=str, required=False,
                                 help="Run eevee substitute in bulk mode.")

    # Complete the activation of ArgumentParser()
    args = parser.parse_args()

    # Conditional switching
    if args.subcommand == "coverage":
        if args.vcf:
            Coverage(sample=args.sample, path_tsv=args.cov_tsv).plot_individual_coverage_variant(variant_table=args.vcf)
        else:
            Coverage(sample=args.sample, path_tsv=args.cov_tsv).plot_individual_coverage()

    elif args.subcommand == "substitute":
        gene_name = str(args.gene).split(":")[0]
        gene_info_str = str(args.gene).split(":")[1].split("-")

        kws_gene = {
            "gene_id": gene_name,
            "gene_info": [int(i) for i in gene_info_str]
        }

        delta = Substitution(reference_fasta_file=args.ref).tabularize(**kws_gene)
        delta.map_variant(variant_table=args.vcf).collect_residue_context()
        delta.write_csv(f"{args.csv_dir}/{args.sample}_residue_table_{gene_name}.csv")
        delta.write_genbank(f"{args.csv_dir}/{args.sample}_variant_{gene_name}.gb")

    elif args.subcommand == "substitute-bulk":
        if args.check:
            config_toml = config_read(args.check)
            config_doctor(config_data=config_toml)
        elif args.run:
            config_toml = config_read(args.run)
            substitution_bulk(config_data=config_toml)

    else:
        print(banner_eevee)
        parser.print_help()


## -----------------------------------------------------------------------------
# main
if __name__ == "__main__":
    main()
