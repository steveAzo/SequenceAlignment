"""
Hybrid MSA script that uses WSL tools from Windows Python
"""
from Bio import AlignIO, SeqIO, Phylo
from Bio.Align import MultipleSeqAlignment
import subprocess
import os
import time
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend for compatibility
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from io import StringIO
import platform

class WSLMSAAnalyzer:
    """
    MSA Analyzer that can use WSL tools from Windows
    """
    
    def __init__(self, input_fasta):
        self.input_fasta = input_fasta
        self.sequences = list(SeqIO.parse(input_fasta, "fasta"))
        self.results = {}
        self.use_wsl = self.check_wsl_available()
        
        print(f"Loaded {len(self.sequences)} sequences for MSA")
        print(f"Platform: {platform.system()}")
        print(f"Using WSL: {self.use_wsl}")
        
    def check_wsl_available(self):
        """
        Check if WSL is available on Windows
        """
        if platform.system() != "Windows":
            return False
            
        try:
            result = subprocess.run(["wsl", "--status"], 
                                  capture_output=True, 
                                  text=True,
                                  timeout=5)
            return result.returncode == 0
        except:
            return False
    
    def check_tool_installed(self, tool_name):
        """
        Check if a tool is installed in WSL
        """
        if not self.use_wsl:
            return False
            
        try:
            result = subprocess.run(
                ["wsl", "which", tool_name],
                capture_output=True,
                text=True,
                timeout=5
            )
            return result.returncode == 0 and result.stdout.strip() != ""
        except:
            return False
    
    def convert_windows_path_to_wsl(self, windows_path):
        """
        Convert Windows path to WSL path
        Example: C:\\Users\\Name\\file.txt -> /mnt/c/Users/Name/file.txt
        """
        # Get absolute path
        abs_path = os.path.abspath(windows_path)
        
        # Convert to WSL format
        # C:\Users\... -> /mnt/c/Users/...
        drive = abs_path[0].lower()
        path_without_drive = abs_path[2:].replace('\\', '/')
        wsl_path = f"/mnt/{drive}{path_without_drive}"
        
        return wsl_path
    
    def run_clustalw_wsl(self, output_base="clustalw_alignment"):
        """
        Run ClustalW via WSL
        """
        print(f"\n{'='*50}")
        print("Running ClustalW via WSL (Progressive Method)")
        print('='*50)
        
        if not self.use_wsl:
            print("✗ WSL not available")
            return None
        
        if not self.check_tool_installed("clustalw2"):
            print("✗ ClustalW not installed in WSL")
            print("  Install with: wsl sudo apt install clustalw")
            return None
        
        output_aln = f"{output_base}.aln"
        
        # Convert paths to WSL format
        wsl_input = self.convert_windows_path_to_wsl(self.input_fasta)
        wsl_output = self.convert_windows_path_to_wsl(output_aln)
        
        try:
            # Build WSL command
            cmd = [
                "wsl",
                "clustalw2",
                f"-INFILE={wsl_input}",
                f"-OUTFILE={wsl_output}",
                "-OUTPUT=CLUSTAL",
                "-OUTORDER=INPUT"
            ]
            
            print(f"Command: {' '.join(cmd)}")
            
            # Time the execution
            start_time = time.time()
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            end_time = time.time()
            
            if result.returncode != 0:
                print(f"ClustalW error: {result.stderr}")
                return None
            
            # Read alignment
            alignment = AlignIO.read(output_aln, "clustal")
            
            # Calculate conservation
            conservation = self.calculate_conservation(alignment)
            
            # Store results
            self.results['clustalw'] = {
                'time': end_time - start_time,
                'alignment': alignment,
                'conservation': conservation,
                'length': alignment.get_alignment_length(),
                'num_sequences': len(alignment)
            }
            
            print(f"✓ ClustalW completed in {end_time - start_time:.2f} seconds")
            print(f"  Alignment length: {alignment.get_alignment_length()}")
            print(f"  Mean conservation: {np.mean(conservation):.3f}")
            
            return alignment
            
        except subprocess.TimeoutExpired:
            print("✗ ClustalW timed out")
            return None
        except Exception as e:
            print(f"✗ ClustalW failed: {str(e)}")
            return None
    
    def run_muscle_wsl(self, output_base="muscle_alignment"):
        """
        Run MUSCLE via WSL
        """
        print(f"\n{'='*50}")
        print("Running MUSCLE via WSL (Iterative Method)")
        print('='*50)
        
        if not self.use_wsl:
            print("✗ WSL not available")
            return None
        
        if not self.check_tool_installed("muscle"):
            print("✗ MUSCLE not installed in WSL")
            print("  Install with: wsl sudo apt install muscle")
            return None
        
        output_aln = f"{output_base}.aln"
        
        # Convert paths to WSL format
        wsl_input = self.convert_windows_path_to_wsl(self.input_fasta)
        wsl_output = self.convert_windows_path_to_wsl(output_aln)
        
        try:
            # Build WSL command
            cmd = [
                "wsl",
                "muscle",
                "-in", wsl_input,
                "-out", wsl_output,
                "-maxiters", "16"
            ]
            
            print(f"Command: {' '.join(cmd)}")
            
            # Time the execution
            start_time = time.time()
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            end_time = time.time()
            
            if result.returncode != 0:
                print(f"MUSCLE error: {result.stderr}")
                return None
            
            # Read alignment
            alignment = AlignIO.read(output_aln, "fasta")
            
            # Calculate conservation
            conservation = self.calculate_conservation(alignment)
            
            # Store results
            self.results['muscle'] = {
                'time': end_time - start_time,
                'alignment': alignment,
                'conservation': conservation,
                'length': alignment.get_alignment_length(),
                'num_sequences': len(alignment)
            }
            
            print(f"✓ MUSCLE completed in {end_time - start_time:.2f} seconds")
            print(f"  Alignment length: {alignment.get_alignment_length()}")
            print(f"  Mean conservation: {np.mean(conservation):.3f}")
            
            return alignment
            
        except subprocess.TimeoutExpired:
            print("✗ MUSCLE timed out")
            return None
        except Exception as e:
            print(f"✗ MUSCLE failed: {str(e)}")
            return None
    
    def run_mafft_wsl(self, output_base="mafft_alignment"):
        """
        Run MAFFT via WSL
        """
        print(f"\n{'='*50}")
        print("Running MAFFT via WSL")
        print('='*50)
        
        if not self.use_wsl:
            print("✗ WSL not available")
            return None
        
        if not self.check_tool_installed("mafft"):
            print("✗ MAFFT not installed in WSL")
            print("  Install with: wsl sudo apt install mafft")
            return None
        
        output_aln = f"{output_base}.aln"
        
        # Convert paths to WSL format
        wsl_input = self.convert_windows_path_to_wsl(self.input_fasta)
        wsl_output = self.convert_windows_path_to_wsl(output_aln)
        
        try:
            # Build WSL command - MAFFT outputs to stdout
            cmd = [
                "wsl",
                "bash", "-c",
                f"mafft --auto --preservecase {wsl_input} > {wsl_output}"
            ]
            
            print(f"Command: {' '.join(cmd)}")
            
            # Time the execution
            start_time = time.time()
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            end_time = time.time()
            
            if result.returncode != 0:
                print(f"MAFFT error: {result.stderr}")
                return None
            
            # Read alignment
            alignment = AlignIO.read(output_aln, "fasta")
            
            # Calculate conservation
            conservation = self.calculate_conservation(alignment)
            
            # Store results
            self.results['mafft'] = {
                'time': end_time - start_time,
                'alignment': alignment,
                'conservation': conservation,
                'length': alignment.get_alignment_length(),
                'num_sequences': len(alignment)
            }
            
            print(f"✓ MAFFT completed in {end_time - start_time:.2f} seconds")
            print(f"  Alignment length: {alignment.get_alignment_length()}")
            print(f"  Mean conservation: {np.mean(conservation):.3f}")
            
            return alignment
            
        except subprocess.TimeoutExpired:
            print("✗ MAFFT timed out")
            return None
        except Exception as e:
            print(f"✗ MAFFT failed: {str(e)}")
            return None
    
    def calculate_conservation(self, alignment):
        """
        Calculate conservation score for each column
        """
        n_seqs = len(alignment)
        aln_length = alignment.get_alignment_length()
        conservation = []
        
        for i in range(aln_length):
            column = alignment[:, i]
            # Remove gaps for calculation
            column_no_gaps = [c for c in column if c != '-']
            if len(column_no_gaps) == 0:
                conservation.append(0)
                continue
            
            # Calculate frequency of most common residue
            from collections import Counter
            counts = Counter(column_no_gaps)
            most_common = counts.most_common(1)[0][1]
            conservation.append(most_common / len(column_no_gaps))
        
        return conservation
    
    def compare_tools(self):
        """
        Compare performance of different MSA tools
        """
        print("\n" + "="*60)
        print("MSA TOOL COMPARISON (via WSL)")
        print("="*60)
        
        # Run all tools
        self.run_clustalw_wsl()
        self.run_muscle_wsl()
        self.run_mafft_wsl()
        
        if not self.results:
            print("\n⚠ No tools completed successfully")
            print("\nTo install tools in WSL, run:")
            print("  wsl sudo apt update")
            print("  wsl sudo apt install -y clustalw muscle mafft")
            return None
        
        # Create comparison table
        comparison = []
        for tool, data in self.results.items():
            comparison.append({
                'Tool': tool.upper(),
                'Time (s)': f"{data['time']:.2f}",
                'Alignment Length': data['length'],
                'Sequences': data['num_sequences'],
                'Mean Conservation': f"{np.mean(data['conservation']):.3f}"
            })
        
        df = pd.DataFrame(comparison)
        print("\nPerformance Comparison:")
        print(df.to_string(index=False))
        
        # Create visualization
        self.create_comparison_plots()
        
        return df
    
    def create_comparison_plots(self):
        """
        Create comparison visualizations
        """
        if len(self.results) == 0:
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        tools = list(self.results.keys())
        colors = ['blue', 'green', 'orange'][:len(tools)]
        
        # Time comparison
        times = [self.results[t]['time'] for t in tools]
        axes[0, 0].bar(tools, times, color=colors)
        axes[0, 0].set_ylabel('Time (seconds)')
        axes[0, 0].set_title('MSA Tool Speed Comparison')
        axes[0, 0].set_xlabel('Tool')
        
        # Conservation profiles
        for i, tool in enumerate(tools):
            axes[0, 1].plot(self.results[tool]['conservation'], 
                          label=tool.upper(), alpha=0.7, color=colors[i])
        axes[0, 1].set_xlabel('Alignment Position')
        axes[0, 1].set_ylabel('Conservation Score')
        axes[0, 1].set_title('Conservation Profiles')
        axes[0, 1].legend()
        axes[0, 1].set_ylim(0, 1)
        
        # Alignment length comparison
        lengths = [self.results[t]['length'] for t in tools]
        axes[1, 0].bar(tools, lengths, color=colors)
        axes[1, 0].set_ylabel('Alignment Length')
        axes[1, 0].set_title('Alignment Length Comparison')
        
        # Mean conservation
        mean_cons = [np.mean(self.results[t]['conservation']) for t in tools]
        axes[1, 1].bar(tools, mean_cons, color=colors)
        axes[1, 1].set_ylabel('Mean Conservation')
        axes[1, 1].set_title('Average Conservation Score')
        axes[1, 1].set_ylim(0, 1)
        
        plt.tight_layout()
        plt.savefig('wsl_msa_comparison.png', dpi=300)
        print(f"\n✓ Comparison plot saved to 'wsl_msa_comparison.png'")

# Main execution
if __name__ == "__main__":
    print("="*70)
    print("WSL-HYBRID MSA ANALYSIS")
    print("Calls WSL tools from Windows Python")
    print("="*70)
    
    # Initialize analyzer
    msa = WSLMSAAnalyzer("hemoglobin_processed.fasta")
    
    if not msa.use_wsl:
        print("\n⚠ WSL not detected!")
        print("\nPlease ensure WSL is installed:")
        print("  1. Open PowerShell as Administrator")
        print("  2. Run: wsl --install")
        print("  3. Restart your computer")
        print("  4. Set up your Linux distribution")
        print("\nOr use the Python-only version: python msa_python_only.py")
    else:
        # Compare tools
        comparison = msa.compare_tools()
        
        if comparison is not None:
            print("\n✓ MSA analysis complete!")
            print("  Results saved to 'wsl_msa_comparison.png'")
        else:
            print("\nFor installation instructions, see: SETUP_MSA_TOOLS.md")
