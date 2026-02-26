from Bio.Seq import Seq 
from Bio import SeqIO 
from Bio.SeqRecord import SeqRecord 
import re  


class SequencePreprocessor:

    def __init__(self, input_file):
        self.input_file = input_file 
        self.sequences = list(SeqIO.parse(input_file, "fasta"))
    

    def remove_gaps(self):
        gap_chars = ['-', '.', '~', '*']

        for record in self.sequences:
            seq_str = str(record.seq)
            for ch in gap_chars:
                seq_str = seq_str.replace(ch, '')
            record.seq = Seq(seq_str) 

        return self 


    def handle_ambiguous(self, action="remove", ambig_chars=None):
        if not ambig_chars:
            ambig_chars = ['B', 'Z', 'X', 'J']
        
        
        for record in self.sequences:
            str_seq = str(record.seq)
            for ch in ambig_chars:
                if ch in str_seq:
                    if action == "remove":
                        str_seq = str_seq.replace(ch, '')
                    elif action == "replace_with_gap":
                        str_seq = str_seq.replace(ch, '-')
            record.seq = Seq(str_seq)
        return self 
    
    def ensure_uniform_length(self, method='trim'):
        if method == 'trim':
            min_length = min(len(sequence.seq) for sequence in self.sequences)
            for sequence in self.sequences:
                if len(sequence.seq) > min_length:
                    sequence.seq = Seq(str(sequence.seq)[:min_length])
        
        elif method == 'pad':
            max_lenth = max(len(sequence.seq) for sequence in self.sequences)

            for sequence in self.sequences:
                exlen = max_lenth - len(sequence.seq)

                if len(sequence.seq) < max_lenth:
                    padding = '-' * (exlen)
                    sequence.seq = Seq(str(sequence.seq) + padding)
        
        return self 
    

    def filter_by_length(self, min_length=50, max_length=1000):
        """Remove sequences outside length range"""
        original_count = len(self.sequences)
        self.sequences = [r for r in self.sequences 
                         if min_length <= len(r.seq) <= max_length]
        print(f"Length filtering: {original_count} -> {len(self.sequences)} sequences")
        return self
    
    def convert_case(self, case='upper'):
        """Convert all sequences to upper or lower case"""
        for record in self.sequences:
            if case == 'upper':
                record.seq = Seq(str(record.seq).upper())
            else:
                record.seq = Seq(str(record.seq).lower())
        return self
    
    def save_processed(self, output_file):
        """Save preprocessed sequences to file"""
        SeqIO.write(self.sequences, output_file, "fasta")
        print(f"Saved {len(self.sequences)} processed sequences to {output_file}")
        return self.sequences

# Usage example
if __name__ == "__main__":
    preprocessor = SequencePreprocessor("haemoglobin_sequences.fasta")
    
    processed_sequences = (preprocessor
        .remove_gaps()
        .handle_ambiguous(action='remove')
        .filter_by_length(min_length=100, max_length=500)
        .convert_case('upper')
        .ensure_uniform_length(method='trim')
        .save_processed("hemoglobin_processed.fasta"))

                    
