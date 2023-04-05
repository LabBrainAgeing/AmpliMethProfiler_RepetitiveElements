=============================================
PACKAGE NAME ---> AmpliMethProfiler 
VERSION ----> 1.2
LICENCE---> GPL
CREATED BY-> Giovanni Scala <giovanni.scala[at]helsinki.fi>
==============================================

INDEX:

-> Requirements
-> About this package
-> Installation
-> Files and their Proper Directories
-> Playing this plugin
-> Known Bugs & Issues
-> Version History
-> Incompatibilities & Save game warnings
-> Credits & Usage


==============================================
REQUIREMENTS:
==============================================
Before running AmpliMethProfiler scripts make sure following python modules are installed:

python >= 2.7
biopython >= 1.65
blast >= 2.2.25
biom-format >= 2.1.5
qiime >= 1.9


==============================================
ABOUT THIS PACKAGE:
===============================================
This package contains python scripts and modules to run AmpliMethProfiler analysis.
Two executable script are provided with this version of AmpliMethProfiler.
  -preprocessFasta.py:
    -> takes as input the the output of an amplicon bisulphite sequencing
       in fasta format along with a series of filtering and demultiplexing parameters
       and generates one filtered fasta file for each sequenced region
  -methProfiles.py:
    -> takes as input a list of regions, the corresponding reference sequences and reads in fasta format and
       a series of filtering parameters. Generates alignment, methylation profiles and summary statistics for
       each input region.
    -qiime_analysis.py:
         -> takes as input a file in BIOM format containing methylation profile abundances for each sample,
            a tab separated file containing samples to analyse along with sample features
            (see the complete docs for info abut the fields of this file), the id of the analysed region, the sample ID,
            the local qiime environment script path, and an output folder name.
             Generates for the input biom three kind of analyses: summary files, Alpha diversity plots and Beta diversity plots.
    -ampliMethProfiler.py:
           -> performs all of the above analyses sequentially. The last one (qiime_analysis) is optional since it requires
              local installation of qiime software.


===============================================
INSTALLATION:
===============================================

PREREQUISITES:
===============================================
AmpliMethProfiler depends on blast alignment tool and biopyton module for the extraction and generation of methylation profiles. 
biom-format and qiime suite are needed to perform epihaplotype analysis through the module qiime_analysis.py.

Users can install each dependency:
  1) independently or 
  2) building an anaconda environment where all dependencies are satisfied (suggested).

1) MANUAL INSTALLATION 
  -to install biopython users can use the PyPi repositories using the command "pip install biopython".
  -to install blast suite users can download the installer from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ 
  -to install biom-format users can use PyPi repositories by first giving the command "pip install numpy" and then "pip install biom-format"
  -to install qiime suite follow the instructions at http://qiime.org/install/alternative.html or MacOS X users can use the MacQuiime bundle deployment from Jeff Werner  lab (http://www.wernerlab.org/software/macqiime/macqiime-installation).

2) (SUGGESTED) SETTING UP AmpliMethProfiler ENVIRONMENT FOR YOUR ANACONDA INSTALLATION 
This is the easiest way to prepare an environment where to run AmpliMethProfiler scripts, it requires a local installation of Anaconda tool (https://www.continuum.io/downloads).
You’ll primarily interact with Anaconda using the conda command.

Create your AmpliMethProfiler environment:
> conda create -n AmpliMethProfilerEnv python=2.7 qiime matplotlib=1.4.3 mock nose biom-format blast biopython -c bioconda
Using QIIME after installation with Anaconda:

Anytime you want to use AmpliMethProfiler after installation with Anaconda, you’ll need to reactivate your AmpliMethProfilerEnv environment using this command:
 > source activate AmpliMethProfilerEnv
When using this configuration, the path to blast executable to give as input to AmpliMethProfiler scripts is found in the bin directory of the created environment (e.g. /home/john/anaconda/envs/AmpliMethProfilerEnv/bin).

To exit the virtual environment, simply run the deactivate command:
        > source deactivate

If you decide later that you don’t want the environment or its packages anymore, deactivate the environment and then run this command:
       > conda remove --name AmpliMethProfilerEnv --all


AmpliMethProfiler SCRIPTS:
===============================================
Extract the three python files from the archive and give them executable privileges.


===============================================
FILES USED AND THEIR PROPER DIRECTORIES:
===============================================
With this distribution are provided five python files.
One master script for the execution of the whole pipeline, three scripts (preprocessFasta.py, methProfiles.py, qiime_analysis.py) for the execution of the single analysis steps and a last python module (methylUtils.py) containing necessary functions to run the scripts.


===============================================
USING AmpliMethProfiler:
===============================================
If installad with Anaconda, remember to reactivate your AmpliMethProfilerEnv environment.
AmpliMethProfiler scripts can be launched:
 -by providing their path to your python interpreter (suggested)
 -placing them into a directory listed in your PATH environment variable and then running them as command line programs.

For detailed instructions on the options and input/output formats of the AmpliMethProfiler scripts we remand the user to the AmpliMethProfiler_pipeline_help documentation contained in this package.


===============================================
KNOWN ISSUES OR BUGS:
===============================================
Fixed bug on bisulphite efficiency evaluation.


===============================================
VERSION HISTORY
===============================================
AmpliMethProfiler 1.0 BETA
AmpliMethProfiler 1.1 BETA
AmpliMethProfiler 1.1
AmpliMethProfiler 1.2



