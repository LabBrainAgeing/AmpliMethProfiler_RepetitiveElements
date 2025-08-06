AmpliMethProfiler for Repetitive Elements
AmpliMethProfiler is a bioinformatic tool developed by Scala et al. for analyzing average DNA methylation and epihaplotypes from deep targeted bisulfite sequencing data.

📌 For installation and usage of the base version, please refer to the original repository:
🔗 https://sourceforge.net/projects/amplimethprofiler/

🔬 About This Fork
Our group optimized AmpliMethProfiler to analyze DNA methylation in Repetitive Elements (RE) such as:

LINE-1

Alu

Ribosomal DNA repeats

This was achieved using an Illumina-based, targeted-deep bisulfite sequencing pipeline.

⚙️ Key Modifications
The base version uses blastn for read alignment. However, its default low-complexity masking (DUST, by Morgulis et al., 10.1089/cmb.2006.13.1028) causes bisulfite-converted unmethylated reads to be filtered out, biasing the DNA methylation profile towards hypermethylation.

We resolved this by integrating a --dust argument into the AmpliMethProfiler command line, allowing users to disable the DUST filtering:

bash
Copia
Modifica
--dust no
This argument uses the _dust_ property from Bio.Blast.Applications.NcbiblastnCommandline.

🚀 Installation & Usage
1. Clone the Repository
bash
Copia
Modifica
git clone https://github.com/LabBrainAgeing/AmpliMethProfiler_RepetitiveElements.git
cd AmpliMethProfiler_RepetitiveElements/AmpliMethProfiler
2. Set Up the Environment
Option A: Follow the original setup guide
🔗 https://sourceforge.net/projects/amplimethprofiler/

Option B: Use our Conda environment

bash
Copia
Modifica
# Download and create the environment
conda env create -f $PWD/AmpliMethProfilerEnv.yml

# Activate it
conda activate AmpliMethProfilerEnv
3. Run AmpliMethProfiler with DUST Disabled
Use the --dust no flag when running the tool.

📄 Example command script:
🔗 Amplimeth_command.sh

🧠 Citation
Please cite the original paper:

Scala, G. et al. "AmpliMethProfiler: a pipeline for the analysis of DNA methylation profiles of targeted bisulfite sequencing of amplicons." BMC Bioinformatics, 2016.
10.1186/s12859-016-1380-3
