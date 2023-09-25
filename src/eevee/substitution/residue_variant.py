from pydantic import BaseModel
from typing import Optional


class TranscriptMeta(BaseModel):
    reference: str        # reference sequence gene ID
    reference_length: int # length of the reference sequence
    tx_nuc_begin: int     # transcript nucleotide start
    tx_nuc_end: int       # transcript nucleotide end


class ResidueVariant(BaseModel):
    """A validated data model when calling the amino acid variant
    """
    ref_pos: int         # variant position
    snp_ref: str         # reference nucleotide
    snp_alt: str         # alternate nucleotide variant
    snp_freq: float      # variant frequency

    tx_nuc_pos: int      # nucleotide position on the transcript reference
    tx_nuc_offset: int   # nucleotide offset (possible value: 1, 2, or 3)
    tx_codon_ref: str    # reference affected codon
    tx_codon_alt: str    # alternate affected codon

    tx_residue_ref: str  # reference residue
    tx_residue_alt: str  # alternate residue
    tx_residue_pos: int  # residue position

    # addditional residue contect information
    residue_context_ref: Optional[str]=None
    residue_context_alt: Optional[str]=None
