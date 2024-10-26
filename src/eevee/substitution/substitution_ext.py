from .residue_variant import ResidueVariant
from Bio.Seq import Seq
import pandas as pd
from typing import List


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
