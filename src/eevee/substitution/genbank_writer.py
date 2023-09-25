from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pandas as pd
from .residue_variant import ResidueVariant, TranscriptMeta
from typing import List


class GBWriter:
    def __init__(self, nuc_table: pd.DataFrame, rv: List[ResidueVariant]):
        """
        """
        nuc_list = nuc_table["nuc"].to_list()

        for v in rv:
            nuc_list[v.ref_pos - 1] = v.snp_alt

        self.variant_gene = "".join(nuc_list)
        self.rv = rv

    def create_variant_record(self, txm: TranscriptMeta):
        """
        """
        sequence_object = Seq(self.variant_gene)
        record = SeqRecord(sequence_object, annotations={"molecule_type": "DNA"})

        # be careful with offsetting the index
        transcript = self.variant_gene[txm.tx_nuc_begin - 1: txm.tx_nuc_end]
        translated = Seq(transcript).translate()

        # be careful with offsetting the index
        translate_feature = SeqFeature(
            FeatureLocation(start=txm.tx_nuc_begin - 1, end=txm.tx_nuc_end),
            type="CDS",
            qualifiers={"translation": translated})

        record.features.append(translate_feature)

        self.record = record
        return self

    def add_annotation(self):
        for v in self.rv:
            q = {"note": v.tx_residue_ref + str(v.tx_residue_pos) + v.tx_residue_alt}
            pos_start, pos_end = v.ref_pos - 1, v.ref_pos
            f = SeqFeature(FeatureLocation(start=pos_start, end=pos_end), type="misc_feature", qualifiers=q)
            self.record.features.append(f)

        return self

    def write_to_file(self, path_gb: str):
        SeqIO.write(self.record, path_gb, "genbank")
