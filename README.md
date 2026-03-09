# DCIT 411: Bioinformatics Project
## Sequence Alignment with Biopython

**Student:** Stephen Azongo - 11027858  
**Date:** February 27, 2026  
**GitHub Repository:** [Github](https://github.com/steveAzo/SequenceAlignment) 
**Google Colab Notebook:** https://drive.google.com/file/d/1j5N4dO-VehqZxpLAdQSa11SyJ-lSFgRL/view?usp=sharing

---

## Project Overview

This project implements and analyzes various sequence alignment techniques using Biopython, including:
- Pairwise sequence alignment (Needleman-Wunsch, Smith-Waterman)
- Multiple sequence alignment (Progressive MSA)
- Comparison with professional tools (ClustalW, MUSCLE, MAFFT)
- Advanced topics (Consensus sequences, conservation analysis, phylogenetics)

---

## Files Description

### Main Files
- **Complete_Sequence_Alignment_Project.ipynb** - Comprehensive Jupyter notebook with all analyses
- **DCIT411_Sequence_Alignment_Report_[YourName].pdf** - Final report with results

### Python Scripts
- **pairwise_alignment.py** - Needleman-Wunsch and Smith-Waterman implementations
- **msa_python_only.py** - Pure Python multiple sequence alignment
- **fetch_sequences.py** - Fetch sequences from NCBI database
- **pre_process_sequence.py** - Sequence preprocessing and cleaning

### Data Files
- **hemoglobin_processed.fasta** - Processed protein sequences
- **haemoglobin_sequences.fasta** - Raw sequences from NCBI

### Generated Results
- **msa_visualization.png** - MSA conservation and gap analysis
- **identity_matrix.png** - Pairwise sequence identity heatmap
- **phylogenetic_tree.png** - UPGMA tree of sequences
- **msa_methods_comparison.png** - Comparison of different MSA tools

---

## Installation and Setup

### Prerequisites
- Python 3.8 or higher
- pip package manager

### Local Installation (Windows/Mac/Linux)

1. Create a virtual environment:
```bash
python -m venv venv
```

2. Activate the virtual environment:
- Windows: `venv\Scripts\activate`
- Mac/Linux: `source venv/bin/activate`

3. Install dependencies:
```bash
pip install -r requirements.txt
```

### Google Colab Setup

1. Upload the notebook to Google Colab
2. Run the installation cell:
```python
!pip install biopython matplotlib pandas numpy
!apt-get update
!apt-get install -y clustalw muscle mafft
```

---

## Usage

### Running the Jupyter Notebook

**Option 1: Google Colab (Recommended)**
- Upload `Complete_Sequence_Alignment_Project.ipynb` to Google Colab
- Run all cells (Runtime → Run all)
- Professional MSA tools will be available

**Option 2: Local Jupyter**
```bash
jupyter notebook Complete_Sequence_Alignment_Project.ipynb
```
Note: Professional MSA tools comparison will be skipped on Windows

### Running Individual Scripts

**Pairwise Alignment:**
```bash
python pairwise_alignment.py
```

**Multiple Sequence Alignment:**
```bash
python msa_python_only.py
```

---

## Key Features

### 1. Pairwise Alignment
- Global alignment using Needleman-Wunsch algorithm
- Local alignment using Smith-Waterman algorithm
- Parameter optimization (substitution matrices, gap penalties)
- Comprehensive comparison and visualization

### 2. Multiple Sequence Alignment
- Progressive alignment with UPGMA guide tree
- Pure Python implementation (no external dependencies)
- Comparison with professional tools (ClustalW, MUSCLE, MAFFT)
- Conservation analysis and gap distribution

### 3. Advanced Analysis
- Consensus sequence generation
- Conserved motif identification
- Phylogenetic tree construction
- Position frequency matrices (profiles)
- Pairwise identity matrices

### 4. Visualization
- Alignment heatmaps
- Conservation profiles
- Gap distribution charts
- Method comparison plots
- Phylogenetic trees

---

## Results Summary

### Pairwise Alignment
- Implemented both global and local alignment algorithms
- Evaluated multiple substitution matrices (BLOSUM62, BLOSUM80, etc.)
- Tested various gap penalty combinations
- Achieved >90% sequence identity for related sequences

### Multiple Sequence Alignment
- Successfully aligned 5 protein sequences
- Identified highly conserved regions (>80% conservation)
- Generated consensus sequences with confidence scores
- Built phylogenetic trees showing evolutionary relationships

### Tool Comparison (Colab Results)
- Compared Python implementation with ClustalW, MUSCLE, and MAFFT
- All methods produced biologically meaningful alignments
- Professional tools showed slightly better conservation scores
- Python implementation competitive in terms of quality

---

## Biological Insights

1. **Conservation Patterns:** Highly conserved regions likely represent functionally important domains
2. **Sequence Diversity:** Moderate pairwise identity (~40-60%) suggests evolutionary divergence
3. **Gap Patterns:** Gaps often correspond to insertions/deletions in evolution
4. **Phylogenetic Relationships:** UPGMA tree reveals clustering of related sequences

---

## Technical Highlights

- **Dynamic Programming:** O(mn) time complexity for pairwise alignment
- **Progressive Alignment:** Efficient MSA construction via guide tree
- **Robust Implementation:** Handles edge cases and missing data
- **Cross-Platform:** Works on Windows, Mac, Linux, and Google Colab
- **Well-Documented:** Comprehensive comments and documentation

---

## Challenges and Solutions

### Challenge 1: External Tool Dependencies
- **Problem:** ClustalW, MUSCLE, MAFFT not available on Windows
- **Solution:** Implemented pure Python MSA alternative

### Challenge 2: Biopython Version Compatibility
- **Problem:** Bio.Align.Applications deprecated in newer versions
- **Solution:** Updated to use subprocess directly

### Challenge 3: Sequence ID Mismatches
- **Problem:** NCBI sequences require version numbers
- **Solution:** Updated IDs to include version (e.g., NP_000509.1)

---

## Future Enhancements

1. Implement iterative refinement algorithms
2. Add profile-based search (PSI-BLAST)
3. Incorporate structural alignment methods
4. Implement Hidden Markov Models (HMMer)
5. Add statistical significance testing
6. Create web interface for interactive analysis

---

## References

1. Needleman & Wunsch (1970) - Global alignment algorithm
2. Smith & Waterman (1981) - Local alignment algorithm
3. Henikoff & Henikoff (1992) - BLOSUM matrices
4. Thompson et al. (1994) - CLUSTAL W algorithm
5. Edgar (2004) - MUSCLE algorithm
6. Cock et al. (2009) - Biopython library

---

## License

This project is submitted as part of DCIT 411 coursework.

---

## Contact

For questions or clarifications, please contact: [Your Email]

---

**Last Updated:** February 27, 2026
