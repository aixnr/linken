import pandas as pd
import numpy as np
import math
from Bio.Seq import Seq


class ReferenceSequenceToPandas:
    def __init__(self, sequence: str, tx_nuc_first: int, tx_nuc_last: int):
        """Turn reference sequence string into Pandas DataFrame.

        Requires 3 inputs from caller:
        - reference sequence (string),
        - position of first translated nucleotide
        - position of last translated nucleotide

        Positional numbering is 1-based.

        Args:
            sequence (str): The reference sequence.
            tx_nuc_first (int): First translated nucleotide.
            tx_nuc_last: Last translated nucleotide.
        """
        self.seq = sequence
        self.tx_nuc_first = tx_nuc_first
        self.tx_nuc_last = tx_nuc_last

    def tabularize(self):
        """Tabularize the data as Pandas DataFrame format, stored as .df attribute

        Using dictionary as intermediate,
          key = index (0-based).
          value[0] is 1-based nucleotide position.
          value[1] is the nucleotide itself.
        """
        position_index_row = range(0, len(self.seq) + 2)
        position_nucleotide = range(1, len(self.seq) + 1)
        column_names = ["ref_pos", "nuc"]

        seq_dict = {
            idx: [pos, nuc] for (idx, pos, nuc) in zip(position_index_row, position_nucleotide, self.seq)
        }

        df = pd.DataFrame.from_dict(seq_dict, orient="index", columns=column_names)
        df["Tx_pos"] = 0
        df["Offset"] = 0
        df["Residue"] = 0

        self.df = df
        return self

    def position_transcript(self):
        """Annotate the transcript position on the reference sequence (1-based)
        """
        cursor = self.tx_nuc_first
        tx_nuc_counter = 1

        for i, row in self.df.iterrows():
            if row["ref_pos"] == cursor:
                if cursor <= self.tx_nuc_last:
                    self.df.at[i, "Tx_pos"] = tx_nuc_counter
                
                cursor += 1
                tx_nuc_counter += 1
        
        return self

    def position_offset(self):
        """Annotate the codon offset (1, 2, or 3) on the translated nucleotide position.
        """
        cursor = self.tx_nuc_first
        offset = 1

        for i, row in self.df.iterrows():
            if row["ref_pos"] == cursor:
                if cursor <= self.tx_nuc_last:
                    if offset < 3:
                        self.df.at[i, "Offset"] = offset
                        offset += 1
                    elif offset == 3:
                        self.df.at[i, "Offset"] = offset
                        offset = 1

                cursor += 1
        
        return self
    
    def position_residue(self):
        """Annotate the amino acid position on the 'Residue' column.
        """
        cursor = self.tx_nuc_first

        for i, row in self.df.iterrows():
            if row["ref_pos"] == cursor:
                if cursor <= self.tx_nuc_last:
                    residue_cursor = (row["Tx_pos"] / 3)
                    self.df.at[i, "Residue"] = math.ceil(residue_cursor)

                cursor += 1

        return self

    def empty_unassigned(self):
        """Change value of 0 to np.NaN
        """
        self.df.replace(0, np.nan, inplace=True)
        return self

    def translate_transcript(self):
        """Translate the transcript into amino acid and save into a separate DataFrame
        """
        transcript_string_series = self.df.dropna(subset="Tx_pos")["nuc"]
        transcript_string = "".join(transcript_string_series.to_list())

        residue = Seq(transcript_string).translate()

        position_index_row = range(0, len(residue) + 2)
        position_residue = range(1, len(residue) + 1)

        residue_dict = {
            idx: [pos, residue] for (idx, pos, residue) in zip(position_index_row, position_residue, residue)
        }

        column_names = ["position", "residue"]
        df = pd.DataFrame.from_dict(residue_dict, orient="index", columns=column_names)
        self.residue = df

        return self
