from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
import time 
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 

class PairwiseAlignerAnalyzer:
    def __init__(self):
        self.aligner = Align.PairwiseAligner()
        self.results = []
    

    def load_sequences(self, fasta_file, seq1_id, seq2_id):
        sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        self.seq1 = sequences[seq1_id]
        self.seq2 = sequences[seq2_id]

        return self 


    def needleman_wunsch_global(self, matrix_name="BLOSUM62", open_gap=-10, extend_gap=-0.5):
        self.aligner.mode = "global"

        self.aligner.extend_gap_score = extend_gap
        self.aligner.open_gap_score = open_gap
        self.aligner.substitution_matrix = substitution_matrices.load(matrix_name)

        start_time = time.time()
        alignments = self.aligner.align(self.seq1.seq, self.seq2.seq)
        end_time = time.time()

        result = {
            'algorithm': 'Needleman-Wunsch (Global)',
            'matrix': matrix_name,
            'gap_open': open_gap,
            'gap_extend': extend_gap,
            'score': alignments.score,
            'time': end_time - start_time,
            'num_alignments': len(alignments),
            'seq1_id': self.seq1.id,
            'seq2_id': self.seq2.id,
            'seq1_len': len(self.seq1),
            'seq2_len': len(self.seq2)
        }

        best_aln = alignments[0]
        result['alignment'] = best_aln
        
        # Calculate identity and similarity
        identity, similarity = self.calculate_identity_similarity(best_aln)
        result['identity'] = identity
        result['similarity'] = similarity
        
        self.results.append(result)
        
        print(f"Score: {result['score']:.2f}")
        print(f"Identity: {identity:.2%}")
        print(f"Similarity: {similarity:.2%}")
        print(f"Time: {result['time']:.4f} seconds")
        print(best_aln)
        
        return result

    
    def smith_waterman_local(self, matrix_name, open_gap=-10, extend_gap=-0.5):
        self.aligner.mode = 'local'
        self.aligner.open_gap_score = open_gap
        self.aligner.extend_gap_score = extend_gap
        self.aligner.substitution_matrix = substitution_matrices.load(matrix_name)

        start_time = time.time()
        alignments = self.aligner.align(self.seq1.seq, self.seq2.seq)
        end_time = time.time()

        result = {
            'algorithm': 'Needleman-Wunsch (Global)',
            'matrix': matrix_name,
            'gap_open': open_gap,
            'gap_extend': extend_gap,
            'score': alignments.score,
            'time': end_time - start_time,
            'num_alignments': len(alignments),
            'seq1_id': self.seq1.id,
            'seq2_id': self.seq2.id,
            'seq1_len': len(self.seq1),
            'seq2_len': len(self.seq2)
        }

        best_aln = alignments[0]
        result['alignment'] = best_aln
        
        # Calculate identity and similarity
        identity, similarity = self.calculate_identity_similarity(best_aln)
        result['identity'] = identity
        result['similarity'] = similarity
        
        self.results.append(result)
        
        print(f"Score: {result['score']:.2f}")
        print(f"Identity: {identity:.2%}")
        print(f"Similarity: {similarity:.2%}")
        print(f"Time: {result['time']:.4f} seconds")
        print(best_aln)

        return result 

    
    def calculate_identity_similarity(self, alignment):
        """
        Calculate percent identity and similarity from alignment
        """
        aligned_seq1, aligned_seq2 = alignment[0], alignment[1]
        
        matches = 0
        similar = 0
        total = len(aligned_seq1)
        
        # Define conservative substitutions (for proteins)
        conservative = {
            'A': ['G', 'S'], 'R': ['K'], 'N': ['Q'], 'D': ['E'],
            'C': ['S'], 'Q': ['N'], 'E': ['D'], 'G': ['A'],
            'H': ['K', 'R'], 'I': ['L', 'V'], 'L': ['I', 'V'],
            'K': ['R', 'H'], 'M': ['L', 'I'], 'F': ['Y', 'W'],
            'P': ['A'], 'S': ['T'], 'T': ['S'], 'W': ['F', 'Y'],
            'Y': ['F', 'W'], 'V': ['I', 'L']
        }
        
        for a, b in zip(aligned_seq1, aligned_seq2):
            if a == b and a != '-':
                matches += 1
                similar += 1
            elif a != '-' and b != '-':
                # Check for conservative substitution
                if a in conservative and b in conservative[a]:
                    similar += 1
                elif b in conservative and a in conservative[b]:
                    similar += 1
        
        identity_pct = matches / total
        similarity_pct = similar / total
        
        return identity_pct, similarity_pct
    
    def experiment_with_matrices(self):
        """Compare different substitution matrices"""
        matrices = ['BLOSUM45', 'BLOSUM62', 'BLOSUM80', 'PAM250']
        results = []
        
        for matrix in matrices:
            print(f"\n=== Testing {matrix} ===")
            # Test with both algorithms
            self.needleman_wunsch_global(matrix_name=matrix)
            self.smith_waterman_local(matrix_name=matrix)
            
        return results
    
    def experiment_with_gap_penalties(self):
        """Compare different gap penalty combinations"""
        gap_combinations = [
            (-5, -1),    # Low penalties
            (-10, -0.5), # Standard
            (-15, -1),   # High open, medium extend
            (-20, -2),   # Very high penalties
            (-8, -4),    # Linear-like (open = extend)
        ]
        
        for gap_open, gap_extend in gap_combinations:
            print(f"\n=== Testing Gap Open={gap_open}, Gap Extend={gap_extend} ===")
            self.needleman_wunsch_global(gap_open=gap_open, gap_extend=gap_extend)
    
    def performance_analysis(self):
        """Analyze algorithm performance with different sequence lengths"""
        lengths = [50, 100, 200, 300, 400, 500]
        global_times = []
        local_times = []
        
        # Create sequences of increasing length
        base_seq = self.seq1.seq
        
        for length in lengths:
            # Truncate sequences
            seq_a = base_seq[:length]
            seq_b = base_seq[:int(length*0.8)]  # 80% length for variation
            
            # Test global
            start = time.time()
            self.aligner.mode = 'global'
            self.aligner.align(seq_a, seq_b)
            global_times.append(time.time() - start)
            
            # Test local
            start = time.time()
            self.aligner.mode = 'local'
            self.aligner.align(seq_a, seq_b)
            local_times.append(time.time() - start)
        
        # Plot results
        plt.figure(figsize=(10, 6))
        plt.plot(lengths, global_times, 'b-o', label='Needleman-Wunsch (Global)')
        plt.plot(lengths, local_times, 'r-o', label='Smith-Waterman (Local)')
        plt.xlabel('Sequence Length')
        plt.ylabel('Time (seconds)')
        plt.title('Algorithm Performance: Time vs Sequence Length')
        plt.legend()
        plt.grid(True)
        plt.savefig('pairwise_performance.png', dpi=300)
        plt.show()
        
        return lengths, global_times, local_times
    
    def generate_report(self):
        """Generate comprehensive report of all alignments"""
        df = pd.DataFrame(self.results)
        
        # Create summary statistics
        summary = df.groupby('algorithm').agg({
            'score': ['mean', 'std', 'min', 'max'],
            'time': ['mean', 'std'],
            'identity': ['mean'],
            'similarity': ['mean']
        })
        
        print("\n" + "="*60)
        print("PAIRWISE ALIGNMENT SUMMARY REPORT")
        print("="*60)
        print(summary)
        
        # Save to CSV
        df.to_csv('pairwise_alignment_results.csv', index=False)
        print("\nDetailed results saved to 'pairwise_alignment_results.csv'")
        
        return df

# Main execution
if __name__ == "__main__":
    # Initialize analyzer
    analyzer = PairwiseAlignerAnalyzer()
    
    # Load sequences (use human and mouse hemoglobin)
    analyzer.load_sequences("hemoglobin_processed.fasta", 
                           "NP_000509.1", "NP_032244.2")  # Human beta vs Mouse alpha
    
    # Run alignments
    print("\n=== STANDARD ALIGNMENTS ===")
    analyzer.needleman_wunsch_global()
    analyzer.smith_waterman_local()
    
    print("\n=== EXPERIMENT 1: Different Matrices ===")
    analyzer.experiment_with_matrices()
    
    print("\n=== EXPERIMENT 2: Different Gap Penalties ===")
    analyzer.experiment_with_gap_penalties()
    
    print("\n=== PERFORMANCE ANALYSIS ===")
    analyzer.performance_analysis()
    
    # Generate final report
    analyzer.generate_report()

        