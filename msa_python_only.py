"""
Multiple Sequence Alignment using only Biopython (no external tools needed)
"""
from Bio import AlignIO, SeqIO, Phylo
from Bio.Align import PairwiseAligner, MultipleSeqAlignment
from Bio.Align import substitution_matrices
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
import time

class ProgressiveMSA:
    """
    Simple progressive multiple sequence alignment implementation
    using Biopython's pairwise aligner
    """
    
    def __init__(self, sequences, matrix_name="BLOSUM62"):
        """
        Initialize with a list of SeqRecord objects
        """
        self.sequences = sequences
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.substitution_matrix = substitution_matrices.load(matrix_name)
        self.aligner.open_gap_score = -10
        self.aligner.extend_gap_score = -0.5
        print(f"Initialized Progressive MSA with {len(sequences)} sequences")
    
    def calculate_pairwise_distances(self):
        """
        Calculate pairwise distances between all sequences
        """
        n = len(self.sequences)
        distances = np.zeros((n, n))
        
        print("\nCalculating pairwise distances...")
        for i in range(n):
            for j in range(i+1, n):
                # Align sequences
                alignments = self.aligner.align(self.sequences[i].seq, 
                                               self.sequences[j].seq)
                score = alignments[0].score
                
                # Convert to distance (higher score = lower distance)
                max_len = max(len(self.sequences[i].seq), len(self.sequences[j].seq))
                distance = 1 - (score / (max_len * 5))  # Normalize
                
                distances[i, j] = distance
                distances[j, i] = distance
        
        return distances
    
    def guide_tree_upgma(self, distances):
        """
        Build a simple UPGMA guide tree from distance matrix
        Returns order of sequences for progressive alignment
        """
        n = len(self.sequences)
        clusters = [[i] for i in range(n)]
        
        while len(clusters) > 1:
            # Find minimum distance between clusters
            min_dist = float('inf')
            merge_i, merge_j = 0, 1
            
            for i in range(len(clusters)):
                for j in range(i+1, len(clusters)):
                    # Average distance between clusters
                    dist_sum = 0
                    count = 0
                    for seq_i in clusters[i]:
                        for seq_j in clusters[j]:
                            dist_sum += distances[seq_i, seq_j]
                            count += 1
                    avg_dist = dist_sum / count
                    
                    if avg_dist < min_dist:
                        min_dist = avg_dist
                        merge_i, merge_j = i, j
            
            # Merge clusters
            new_cluster = clusters[merge_i] + clusters[merge_j]
            clusters = [c for idx, c in enumerate(clusters) 
                       if idx not in [merge_i, merge_j]]
            clusters.append(new_cluster)
        
        return clusters[0]
    
    def align_profile_to_sequence(self, profile_alignment, sequence):
        """
        Align a profile (alignment) to a single sequence
        """
        # Use consensus sequence for the profile
        consensus = self.get_consensus(profile_alignment)
        
        # Align consensus to new sequence
        alignments = self.aligner.align(consensus, sequence.seq)
        best_alignment = alignments[0]
        
        # Extract aligned strings
        aligned_consensus = str(best_alignment[0])
        aligned_sequence = str(best_alignment[1])
        
        # Expand profile alignment with gaps
        new_profile = []
        cons_pos = 0
        
        for pos in range(len(aligned_consensus)):
            if aligned_consensus[pos] == '-':
                # Insert gap in all profile sequences
                for record in profile_alignment:
                    record.seq = Seq(str(record.seq) + '-')
            else:
                cons_pos += 1
        
        # Add the new sequence with proper gaps
        new_seq = SeqRecord(
            Seq(aligned_sequence),
            id=sequence.id,
            description=sequence.description
        )
        
        # Normalize all sequences to same length
        max_len = max(len(str(r.seq)) for r in profile_alignment)
        max_len = max(max_len, len(aligned_sequence))
        
        for record in profile_alignment:
            if len(str(record.seq)) < max_len:
                record.seq = Seq(str(record.seq) + '-' * (max_len - len(str(record.seq))))
        
        if len(aligned_sequence) < max_len:
            new_seq.seq = Seq(aligned_sequence + '-' * (max_len - len(aligned_sequence)))
        
        profile_alignment.append(new_seq)
        
        return profile_alignment
    
    def get_consensus(self, alignment):
        """
        Get consensus sequence from alignment
        """
        consensus = []
        aln_length = len(str(alignment[0].seq))
        
        for pos in range(aln_length):
            column = [str(rec.seq)[pos] for rec in alignment if str(rec.seq)[pos] != '-']
            if column:
                counts = Counter(column)
                consensus.append(counts.most_common(1)[0][0])
            else:
                consensus.append('-')
        
        return Seq(''.join(consensus))
    
    def progressive_alignment(self):
        """
        Perform progressive multiple sequence alignment
        """
        print("\n" + "="*50)
        print("Progressive Multiple Sequence Alignment")
        print("="*50)
        
        start_time = time.time()
        
        # Step 1: Calculate pairwise distances
        distances = self.calculate_pairwise_distances()
        
        # Step 2: Build guide tree
        print("Building guide tree (UPGMA)...")
        order = self.guide_tree_upgma(distances)
        
        # Step 3: Progressive alignment following guide tree
        print("Performing progressive alignment...")
        
        # Start with first two sequences
        first_pair = self.aligner.align(
            self.sequences[order[0]].seq,
            self.sequences[order[1]].seq
        )[0]
        
        # Create initial alignment
        alignment = [
            SeqRecord(Seq(str(first_pair[0])), 
                     id=self.sequences[order[0]].id,
                     description=self.sequences[order[0]].description),
            SeqRecord(Seq(str(first_pair[1])), 
                     id=self.sequences[order[1]].id,
                     description=self.sequences[order[1]].description)
        ]
        
        # Add remaining sequences
        for idx in order[2:]:
            print(f"  Adding sequence {idx+1}/{len(self.sequences)}: {self.sequences[idx].id}")
            # Align profile to new sequence - simple pairwise approach
            new_seq = self.sequences[idx]
            
            # Align new sequence to consensus
            consensus = self.get_consensus(alignment)
            pair_aln = self.aligner.align(consensus, new_seq.seq)[0]
            
            # Add gaps to existing alignment where needed
            aligned_consensus = str(pair_aln[0])
            aligned_new = str(pair_aln[1])
            
            # Rebuild alignment with proper gaps
            new_alignment = []
            gap_positions = [i for i, c in enumerate(aligned_consensus) 
                           if c == '-' and aligned_new[i] != '-']
            
            # Add gaps to all existing sequences
            for record in alignment:
                seq_str = str(record.seq)
                # Insert gaps at appropriate positions
                for gap_pos in sorted(gap_positions):
                    seq_str = seq_str[:gap_pos] + '-' + seq_str[gap_pos:]
                new_alignment.append(SeqRecord(
                    Seq(seq_str),
                    id=record.id,
                    description=record.description
                ))
            
            # Add new sequence
            new_alignment.append(SeqRecord(
                Seq(aligned_new),
                id=new_seq.id,
                description=new_seq.description
            ))
            
            alignment = new_alignment
        
        end_time = time.time()
        
        # Convert to MultipleSeqAlignment
        msa = MultipleSeqAlignment(alignment)
        
        print(f"\n✓ Alignment completed in {end_time - start_time:.2f} seconds")
        print(f"  Alignment length: {msa.get_alignment_length()}")
        print(f"  Number of sequences: {len(msa)}")
        
        return msa, end_time - start_time


class MSAAnalyzer:
    """
    Analyze multiple sequence alignments
    """
    
    def __init__(self, alignment, computation_time):
        self.alignment = alignment
        self.computation_time = computation_time
    
    def calculate_conservation(self):
        """
        Calculate conservation score for each column
        """
        n_seqs = len(self.alignment)
        aln_length = self.alignment.get_alignment_length()
        conservation = []
        
        for i in range(aln_length):
            column = self.alignment[:, i]
            # Remove gaps for calculation
            column_no_gaps = [c for c in column if c != '-']
            if len(column_no_gaps) == 0:
                conservation.append(0)
                continue
            
            # Calculate frequency of most common residue
            counts = Counter(column_no_gaps)
            most_common = counts.most_common(1)[0][1]
            conservation.append(most_common / len(column_no_gaps))
        
        return conservation
    
    def calculate_gap_statistics(self):
        """
        Calculate gap-related statistics
        """
        gap_stats = []
        
        for record in self.alignment:
            seq_str = str(record.seq)
            gap_count = seq_str.count('-')
            gap_percent = (gap_count / len(seq_str)) * 100
            
            gap_stats.append({
                'Sequence': record.id,
                'Total Gaps': gap_count,
                'Gap %': f"{gap_percent:.1f}%"
            })
        
        return pd.DataFrame(gap_stats)
    
    def calculate_pairwise_identity(self):
        """
        Calculate pairwise identity between all sequences
        """
        n = len(self.alignment)
        identities = []
        
        for i in range(n):
            for j in range(i+1, n):
                seq_i = str(self.alignment[i].seq)
                seq_j = str(self.alignment[j].seq)
                
                matches = sum(1 for a, b in zip(seq_i, seq_j) 
                            if a == b and a != '-')
                length = sum(1 for a, b in zip(seq_i, seq_j) 
                           if a != '-' or b != '-')
                
                if length > 0:
                    identity = (matches / length) * 100
                    identities.append({
                        'Seq1': self.alignment[i].id,
                        'Seq2': self.alignment[j].id,
                        'Identity %': f"{identity:.1f}%"
                    })
        
        return pd.DataFrame(identities)
    
    def visualize_alignment(self, filename='msa_visualization.png'):
        """
        Create comprehensive visualization of the alignment
        """
        conservation = self.calculate_conservation()
        
        fig, axes = plt.subplots(3, 1, figsize=(14, 10))
        
        # 1. Conservation profile
        axes[0].plot(conservation, color='blue', linewidth=1.5)
        axes[0].fill_between(range(len(conservation)), conservation, alpha=0.3)
        axes[0].set_xlabel('Alignment Position')
        axes[0].set_ylabel('Conservation Score')
        axes[0].set_title('Sequence Conservation Profile')
        axes[0].set_ylim(0, 1)
        axes[0].grid(True, alpha=0.3)
        
        # 2. Heat map of alignment
        aln_array = []
        aa_to_num = {
            'A': 1, 'R': 2, 'N': 3, 'D': 4, 'C': 5,
            'Q': 6, 'E': 7, 'G': 8, 'H': 9, 'I': 10,
            'L': 11, 'K': 12, 'M': 13, 'F': 14, 'P': 15,
            'S': 16, 'T': 17, 'W': 18, 'Y': 19, 'V': 20,
            '-': 0
        }
        
        for record in self.alignment:
            seq_nums = [aa_to_num.get(aa, 0) for aa in str(record.seq)]
            aln_array.append(seq_nums)
        
        im = axes[1].imshow(aln_array, aspect='auto', cmap='viridis', interpolation='nearest')
        axes[1].set_xlabel('Alignment Position')
        axes[1].set_ylabel('Sequence')
        axes[1].set_title('Alignment Heatmap (darker = gaps)')
        axes[1].set_yticks(range(len(self.alignment)))
        axes[1].set_yticklabels([rec.id for rec in self.alignment], fontsize=8)
        plt.colorbar(im, ax=axes[1], label='Amino Acid')
        
        # 3. Gap distribution
        gap_counts = [str(rec.seq).count('-') for rec in self.alignment]
        seq_ids = [rec.id for rec in self.alignment]
        
        axes[2].barh(seq_ids, gap_counts, color='coral')
        axes[2].set_xlabel('Number of Gaps')
        axes[2].set_ylabel('Sequence')
        axes[2].set_title('Gap Distribution per Sequence')
        axes[2].grid(True, alpha=0.3, axis='x')
        
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"\n✓ Visualization saved to '{filename}'")
        plt.show()
    
    def print_alignment_excerpt(self, start=0, width=60):
        """
        Print a portion of the alignment for inspection
        """
        print("\n" + "="*70)
        print("ALIGNMENT EXCERPT")
        print("="*70)
        
        aln_length = self.alignment.get_alignment_length()
        end = min(start + width, aln_length)
        
        for record in self.alignment:
            seq_str = str(record.seq)[start:end]
            print(f"{record.id:20s} {seq_str}")
        
        print(f"\nShowing positions {start+1}-{end} of {aln_length}")
    
    def generate_summary_report(self):
        """
        Generate comprehensive summary report
        """
        print("\n" + "="*70)
        print("MULTIPLE SEQUENCE ALIGNMENT SUMMARY REPORT")
        print("="*70)
        
        print(f"\nComputation Time: {self.computation_time:.2f} seconds")
        print(f"Number of Sequences: {len(self.alignment)}")
        print(f"Alignment Length: {self.alignment.get_alignment_length()}")
        
        # Conservation statistics
        conservation = self.calculate_conservation()
        print(f"\nConservation Statistics:")
        print(f"  Mean Conservation: {np.mean(conservation):.3f}")
        print(f"  Max Conservation: {np.max(conservation):.3f}")
        print(f"  Min Conservation: {np.min(conservation):.3f}")
        print(f"  Std Conservation: {np.std(conservation):.3f}")
        
        # Gap statistics
        print(f"\nGap Statistics:")
        gap_df = self.calculate_gap_statistics()
        print(gap_df.to_string(index=False))
        
        # Pairwise identity
        print(f"\nPairwise Identity Matrix:")
        identity_df = self.calculate_pairwise_identity()
        print(identity_df.to_string(index=False))
        
        # Print excerpt of alignment
        self.print_alignment_excerpt()
    
    def export_alignment(self, filename="alignment_output.aln", format="clustal"):
        """
        Export alignment to file
        """
        AlignIO.write(self.alignment, filename, format)
        print(f"\n✓ Alignment exported to '{filename}' (format: {format})")


# Main execution
if __name__ == "__main__":
    print("="*70)
    print("PYTHON-ONLY MULTIPLE SEQUENCE ALIGNMENT")
    print("(No external tools required)")
    print("="*70)
    
    # Load sequences
    input_file = "hemoglobin_processed.fasta"
    sequences = list(SeqIO.parse(input_file, "fasta"))
    print(f"\nLoaded {len(sequences)} sequences from '{input_file}'")
    
    # Perform progressive MSA
    progressive_msa = ProgressiveMSA(sequences)
    alignment, comp_time = progressive_msa.progressive_alignment()
    
    # Analyze results
    analyzer = MSAAnalyzer(alignment, comp_time)
    analyzer.generate_summary_report()
    
    # Visualize
    analyzer.visualize_alignment('python_msa_visualization.png')
    
    # Export
    analyzer.export_alignment('hemoglobin_python_msa.aln', 'clustal')
    analyzer.export_alignment('hemoglobin_python_msa.fasta', 'fasta')
    
    print("\n" + "="*70)
    print("✓ MSA Analysis Complete!")
    print("="*70)
    print("\nOutput files:")
    print("  - python_msa_visualization.png (visualization)")
    print("  - hemoglobin_python_msa.aln (Clustal format)")
    print("  - hemoglobin_python_msa.fasta (FASTA format)")
