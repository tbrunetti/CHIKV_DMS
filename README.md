# CHIKV_DMS
Deep mutational scanning repository for Thomas (Tem) Morrison Lab  

This repository contains code to run analysis for deep mutational scanning projects related to analysis for CHIKV.  

## Table of Contents
---------------------
- [CHIKV\_DMS](#chikv_dms)
  - [Table of Contents](#table-of-contents)
    - [Introduction and Overview](#introduction-and-overview)
    - [Software Requirements and Dependencies](#software-requirements-and-dependencies)
    - [Installation](#installation)
    - [General Usage](#general-usage)

* script: [plot_mutational_frequency_and_qc_stats.py](https://github.com/tbrunetti/CHIKV_DMS/wiki/DMS_plots)  
  


### Introduction and Overview
-----------------------------
The repository contains code for analyzing deep mutational scanning projects, particularly related to CHIKV.  


### Software Requirements and Dependencies
------------------------------------------  

For script `plot_mutational_frequency_and_qc_stats.py`:  

* [python >= v3.9](https://www.python.org/downloads/)  
* [pandas](https://pandas.pydata.org/docs/getting_started/install.html#installing-from-pypi)  
* [numpy](https://numpy.org/install/)  
* [matplotlib](https://matplotlib.org/stable/users/installing/index.html#installation)
* [logomaker](https://pypi.org/project/logomaker/) 


### Installation
-----------------  
There is no installation required, please just clone the repository:  
```
git clone https://github.com/tbrunetti/CHIKV_DMS.git
```  

Then cd into the `code` directory where all scripts will be present:  
```
cd code  
```

### General Usage  
-----------------  

For using the `logo_plot_standalone.py` script, general use is as follows:  

```
python3 logo_plot_standalone.py --input wtDNA_filtered_df1_dedup.csv --sampleName wtDNA --annotConfig ../ref/annotations_config.csv --codonStartPos 9 
```
<br/>  
For all possible arguments available, you can run the following:  

```
python3 logo_plot_standalone.py --help
```
,which will show the following options and their defaults:  

```
usage: logo_plot_standalone.py [-h] --input INPUT [--sampleName SAMPLENAME] [--annotConfig ANNOTCONFIG] [--codonStartPos CODONSTARTPOS] [--codonEndPos CODONENDPOS]
                               [--aaSpacing AASPACING] [--minAnnotLabel MINANNOTLABEL]

Generates logo plot for predefined input matrix

options:
  -h, --help            show this help message and exit  
  --input INPUT         Path to input csv matrix containing data to plot (default: None)  
  --sampleName SAMPLENAME  
                        string indicating the name to give to sample (default: sample_1)  
  --annotConfig ANNOTCONFIG  
                        Path to csv containing annotations. Example file located in ref folder of github repo (default: None)  
  --codonStartPos CODONSTARTPOS  
                        The position of which codon position you want to start at (must be present in your csv matrix provided to --input; default is to plot every position in
                        your matrix) (default: None)  
  --codonEndPos CODONENDPOS  
                        The position of which codon position you want to end at (must be present in your csv matrix provided to --input; default is to plot every position in your
                        matrix) (default: None)  
  --aaSpacing AASPACING  
                        the number of amino acids to show per line on the logo plot (default: 65)  
  --minAnnotLabel MINANNOTLABEL  
                        the minimum length of consecutive amino acids under an annotation bar; anything smaller (non-inclusive) than this value will not have text written in the
                        bar, to help prevent text from overflowing into margins (default: 7)  

```

 

