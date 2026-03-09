# Setting Up Multiple Sequence Alignment Tools

## Option 1: Using WSL (Windows Subsystem for Linux) - RECOMMENDED

Since you have WSL installed, this is the best approach for using professional MSA tools.

### Step 1: Install the tools in WSL

Open your WSL terminal and run:

```bash
# Update package list
sudo apt update

# Install ClustalW
sudo apt install -y clustalw

# Install MUSCLE
sudo apt install -y muscle

# Install MAFFT
sudo apt install -y mafft
```

### Step 2: Verify installations

```bash
# Check if tools are installed
clustalw2 -help
muscle -version
mafft --version
```

### Step 3: Run the analysis from WSL

Option A - Copy files to WSL and run there:
```bash
# Navigate to your WSL home directory
cd ~

# Create a directory for your project
mkdir dna-alignment
cd dna-alignment

# Copy files from Windows to WSL
# (Windows drives are mounted at /mnt/c, /mnt/d, etc.)
cp /mnt/c/Users/User/Desktop/DNA-sequnce-alignment/*.py .
cp /mnt/c/Users/User/Desktop/DNA-sequnce-alignment/*.fasta .

# Install Python packages in WSL
pip install biopython matplotlib pandas numpy

# Run the analysis
python msa_analysis.py
```

Option B - Run directly from Windows filesystem:
```bash
# Navigate to your Windows folder from WSL
cd /mnt/c/Users/User/Desktop/DNA-sequnce-alignment

# Make sure Python packages are installed
pip install biopython matplotlib pandas numpy

# Run the script
python msa_analysis.py
```

---

## Option 2: Installing on Windows (More Complex)

### ClustalW for Windows
1. Download from: http://www.clustal.org/download/current/
2. Extract the executable (clustalw2.exe)
3. Add to PATH or place in script directory

### MUSCLE for Windows
1. Download from: https://drive5.com/muscle/downloads.htm
2. Rename to `muscle.exe`
3. Add to PATH or place in script directory

### MAFFT for Windows
1. Download from: https://mafft.cbrc.jp/alignment/software/windows.html
2. Install and add to PATH

### Add to Windows PATH
```powershell
# Add to current session
$env:PATH += ";C:\path\to\your\tools"

# Add permanently (requires admin)
[Environment]::SetEnvironmentVariable("PATH", "$env:PATH;C:\path\to\your\tools", "Machine")
```

---

## Option 3: Use WSL from Windows Python

I can create a hybrid script that calls WSL tools from Windows PowerShell!

---

## Recommended Approach

**For learning and occasional use:** Use the Python-only implementation (`msa_python_only.py`)

**For serious research:** Set up WSL properly and use the professional tools (`msa_analysis.py`)

**For convenience:** Use the WSL hybrid approach (see `msa_wsl_hybrid.py`)

## Troubleshooting

### "Command not found" in WSL
- Make sure you ran `sudo apt update` first
- Try `which clustalw2` to see if it's installed
- Some distros might have different package names

### Permission denied
- Use `chmod +x filename` to make files executable
- Run with `sudo` if needed for installation

### Python packages missing in WSL
```bash
# Install pip if needed
sudo apt install python3-pip

# Install required packages
pip3 install biopython matplotlib pandas numpy
```

### Display issues with plots in WSL
```bash
# Install X server for Windows (VcXsrv or Xming)
# Or save plots without displaying:
# Add this to your Python script:
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend
```
