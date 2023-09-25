import sys
import pandas as pd
from pyfaidx import Fasta
from Bio.Seq import Seq
from .reference_table import ReferenceSequenceToPandas
from .residue_variant import ResidueVariant, TranscriptMeta
from .genbank_writer import GBWriter
from typing import Optional, List


class Substitution:
    """Substitution class converts nucleotide variants to amino acid

    Relies on .reference_table.ReferenceSequenceToPandas class and
      .residue_variant.ResidueVariant class
    """
    def __init__(self, reference_fasta_file: str):
        """Open a reference fasta file with pyfaidx

        Args:
          reference_fasta_file (str): path to the reference fasta
        """
        reference_gene = Fasta(reference_fasta_file)
        self.gene = reference_gene

    def tabularize(self, gene_id: Optional[str]=None, gene_info: Optional[List[int]]=None):
        """Uses ReferenceSequenceToPandas class.
        When gene_id is supplied, only valid for 1 gene ID per run.
        Of note, the gene from a Fasta class is not a str type but pyfaidx.FastaRecord.
        Use str() to convert it to python string primitive

        Args:
          gene_id (str): only process for the supplied gene_id
          gene_info (list): a list of 2 ints for the first [0] and last [1] of the 
            transated nucleotide.
        
        Return:
          Self:
            self.nuc_table (pd.DataFrame class)
            self.residue_table (pd.DataFrame class)
            self.gene_data (TranscriptMeta class)
        """
        if gene_id not in self.gene.keys():
            print(f"[ERR] Error, the gene '{gene_id}' is not in the reference genome...")
            print("[ERR] Program exited.")
            sys.exit()

        if not gene_id:
            self.run_mode = "multi"

            print("Currently not supported")
            sys.exit()

        if gene_id:
            self.run_mode = "single"

            gene_string_sequence = str(self.gene[gene_id])
            gene_string_length = len(gene_string_sequence)

            gene_meta = {
                "reference": gene_id,
                "reference_length": gene_string_length,
                "tx_nuc_begin": gene_info[0],
                "tx_nuc_end": gene_info[1],
            }

            tx_meta = TranscriptMeta(**gene_meta)

            gene = {
                "sequence": gene_string_sequence,
                "tx_nuc_first": gene_info[0],
                "tx_nuc_last": gene_info[1],
            }

            nuc_table = (ReferenceSequenceToPandas(**gene)
                         .tabularize()
                         .position_transcript().position_offset().position_residue()
                         .empty_unassigned()
                         .translate_transcript())

            # Save to self: nuc_table and residue_table
            self.nuc_table = nuc_table.df
            self.residue_table = nuc_table.residue
            self.gene_data = tx_meta
            return self

    def map_variant(self, variant_table: str):
        """Map nucleotide variant to amino acid

        Args:
          variant_table (str):
        """
        # match the column names to the self.nuc_table for merging later
        column_mapper = {
            "CHROM": "Chr", "POS": "ref_pos", "REF": "nuc", "ALT": "alt", "AF": "freq"
        }

        _df_variant = pd.read_csv(variant_table, sep="\t").rename(columns=column_mapper)

        if self.run_mode == "single":
            self.residue_variant = collect_variant(self.gene_data.reference, _df_variant, self.nuc_table)

        elif self.run_mode == "multi":
            print("Not yet implemented")
            sys.exit()

        return self
    
    def collect_residue_context(self):
        """Amino acid itself does not make it easier to identify the region.
        A context (string of 3 residues, e.g., CFP where the change is in the middle, 'Y')
          does help to locate where the actual change is.

        Return:
          Updates the self.residue_table (ResidueVariant class) for
            tx_residue_n1 and tx_residue_n3 fields.
        """
        residue_table = self.residue_table

        for variant in self.residue_variant:
            pos = variant.tx_residue_pos
            pos_1 = pos - 1
            pos_3 = pos + 1

            residue_n1 = residue_table.query(f" position == {pos_1} ")["residue"].to_list()[0]
            residue_n3 = residue_table.query(f" position == {pos_3} ")["residue"].to_list()[0]

            variant.residue_context_ref = residue_n1 + variant.tx_residue_ref + residue_n3
            variant.residue_context_alt =  residue_n1 + variant.tx_residue_alt + residue_n3

        return self

    def write_csv(self, path_csv: str):
        """Write the result into CSV

        Args:
          path_csv (str): output location
        
        Return:
          None
        """
        collect_residue_dict = {}
        for variant, n in zip(self.residue_variant, range(len(self.residue_variant))):
            collect_residue_dict[n] = variant.model_dump()

        df = pd.DataFrame.from_dict(collect_residue_dict, orient="index")
        df.to_csv(path_csv, index=None)

    def write_genbank(self, path_gb: str):
        """Write into annotated genbank file.
        Writing for the variant sequence, not for the reference sequence.

        Return:
          None
        """
        (GBWriter(self.nuc_table, self.residue_variant)
         .create_variant_record(self.gene_data)
         .add_annotation()
         .write_to_file(path_gb))


def collect_variant(gene_id: str, df_vcf: pd.DataFrame, df_nuc_table: pd.DataFrame = None) -> List[ResidueVariant]:
    """Get nucleotide codon for a given nucleotide variant

    Args:
      gene_id (str): the reference gene ID
      df_vcf (pd.DataFrame): tabular data from VCF file containing called SNPs
      df_nuc_table (pd.DataFrame): nucleotide reference table from ReferenceSequenceToPandas

    Return:
      ResidueVariant class
    """
    df_vcf_gene = (df_vcf
                   .query(f" Chr == '{gene_id}' ")
                   .reset_index(drop=True)
                   .drop(columns="Chr"))

    # only collect those in the coding sequence
    df_merged = df_nuc_table.merge(df_vcf_gene, how="outer").dropna(subset="Tx_pos")

    substitution_collected: List[ResidueVariant] = []
    for i, row in df_merged.iterrows():
        # pd.isnull() to only iterate over variants, where row["alt"] contains value
        if not pd.isnull(row["alt"]):
            rv = collect_substitution(df=df_merged, i=i, offset=row["Offset"])
            substitution_collected.append(rv)

    return substitution_collected


def collect_substitution(df: pd.DataFrame, i: int, offset: int) -> ResidueVariant:
    """Collect the amino acid substitution and return the data into ResidueVariant class.
    
    This function uses a sub-function extractor() to reduce code duplication.
    """

    def extractor(offset: int, i: int, df: pd.DataFrame, p1: str, p2: str, p3: str) -> ResidueVariant:
        """Extracts information from the called variants/SNPs to dump into ResidueVariant class

        Args:
          offset (int): codon offset (1, 2, or 3)
          i (int): index from the iterator, called by pd.DataFrame.iterrows()
          df (pd.DataFrame): the dataframe which contains merged nucleotide table and variant table
          p1 (str): nucleotide at codon position 1
          p2 (str): nucleotide at codon position 2
          p3 (str): nucleotide at codon position 3

        Return:
          ResidueVariant class
        """
        alt = df.at[i, "alt"]
        codon_ref = p1 + p2 + p3
        snp_ref = ""
        codon_alt = ""

        if offset == 1:
            snp_ref = p1
            codon_alt = alt + p2 + p3
        elif offset == 2:
            snp_ref = p2
            codon_alt = p1 + alt + p3
        elif offset == 3:
            snp_ref = p3
            codon_alt = p1 + p2 + alt

        tx_residue_ref = Seq(codon_ref).translate()
        tx_residue_alt = Seq(codon_alt).translate()

        substitution_data = {
            "ref_pos": df.at[i, "ref_pos"],
            "snp_ref": snp_ref,
            "snp_alt": alt,
            "snp_freq": df.at[i, "freq"],
            "tx_nuc_pos": df.at[i, "Tx_pos"],
            "tx_nuc_offset": df.at[i, "Offset"],
            "tx_codon_ref": codon_ref,
            "tx_codon_alt": codon_alt,
            "tx_residue_ref": str(tx_residue_ref),
            "tx_residue_alt": str(tx_residue_alt),
            "tx_residue_pos": df.at[i, "Residue"]
        }

        return ResidueVariant(**substitution_data)

    if offset == 1:
        codon_pos_1 = df.at[i, "nuc"]
        codon_pos_2 = df.at[i + 1, "nuc"]
        codon_pos_3 = df.at[i + 2, "nuc"]

        return extractor(1, i, df, codon_pos_1, codon_pos_2, codon_pos_3)

    elif offset == 2:
        codon_pos_1 = df.at[i - 1, "nuc"]
        codon_pos_2 = df.at[i, "nuc"]
        codon_pos_3 = df.at[i + 1, "nuc"]

        return extractor(2, i, df, codon_pos_1, codon_pos_2, codon_pos_3)

    elif offset == 3:
        codon_pos_1 = df.at[i - 2, "nuc"]
        codon_pos_2 = df.at[i - 1, "nuc"]
        codon_pos_3 = df.at[i, "nuc"]

        return extractor(3, i, df, codon_pos_1, codon_pos_2, codon_pos_3)
