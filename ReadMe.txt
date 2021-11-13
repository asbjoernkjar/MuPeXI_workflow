###### Workflow to run mupexi on cluster ######

This workflow calls and ranks neoantigens for all mutations in supplied vcf files.
VCF FILES HAS TO BE HG38 --> USE LIFTOVER IF HG19
Only mutations with PASS filter are used (manually change to PASS for inclusion)

Results are merged into two files (see below)
File with all neoantigens
Files with neoantigens filtered based on RankEL or BA (for easier handling)

netMHCpan4.1b is used for ranking neoantigens
MuPeXI updated to run with netMHCpan4.1b is used to extract neoantigens

### MuPexi ###
This workflow uses MuPeXI updated for netMHCpan4.1b
This can be found at https://github.com/asbjoernkjar/MuPeXI.git

### netMHCpan4.1b ###
Can be found at DTU website:
https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1

Both should be placed in main folder of workflow
They will be setup automatically 


###### RUNNING THE PIPELINE #######

As any snakemake workflow (tested for 6.8.0)
Snakemake –profile config/slurm -j10

Note  change account in config/slurm/cluster.yaml

###### Setting up the Pipeline #######

Controlled by file: config/config.yaml

### VEP ###

Vep 104.3 is downloaded from anaconda (much more convinient that manual install)
To run vep via mupexi a bash file is run, that just calls the anaconda version:
    #!/usr/bin/bash
    vep $@"

This file should then be referenced in "PATH_VEP" in the config file

The program still needs a VEP cache
The VEP cache can be downloaded via: 
$vep_install -a cf -s homo_sapiens -y GRCh38 -c /OUT/PATH/ --CACHE_VERSION 104 --CONVERT

To do this either activate the enviroment created by snakemake 
Or create a new environment and install ensembl-vep=104.3 from bioconda

### cDNA/pep/cosmic ###

Mupexi uses reference cDNA/pep files + cosmic annotation file
See mupexi user manual for thorough instructions 
CDNA/pep can be downloaded from Ensemble & census from cosmic website 

NOTE: THE VERSION SHOULD MATCH VEP VERSION (104)


### Sample_File ###

The samples are specified by csv file. It should have the following columns:
Sample_ID --> Name of samples, will be used as output names
VCF_file_name --> Name of the bam file for that sample. 
	PATH_VCFS + bam_file_name is used to access files
Race --> race of individual {Caucasian, Black, Asian, Unknown}
	If there is no Race column Unknown is used by default 

NOTE --> the pipeline drops any row with a NA value in any col!!!

### HLA_CALLS ###

CSV file with HLA calls in correct 4-digit format: (HLA-A03:01)
The file should have 8 columns:
Sample_ID --> Name of samples (same as in Sample_file)
Caller --> Name of caller used --> config specifies which caller to use
HLA_A1, HLA_A2, HLA_B1, HLA_B2, HLA_C1, HLA_C2 --> with HLA calls

### SPLIT VCFS ###

MuPeXI can take a very long time to run with large vcf files. (>9H for 5000 mutations)
This time seems to increase non-linearly
This workflow can split VCF files into smaller subfiles and run on these.
NOTE this will filter out any non PASS mutations
SPLIT_SIZE describes size of vcf files
TOLERANCE describes the smallest VCF file allowed
MuPeXI WILL CRASH if no protein altering mutations are found
Having high tolerance should fix this problem

(Rare cases can be “solved” by adding an empty output file with the right name
to the results folder (needs mupexi file header!))

set "SPLIT_VCFS: false", to run without splitting

If SPLIT_VCFS = true, the first run will create a list of split files to be created
This opens all vcf files on the local node. You might want to do this on 
a cluster node with a  dryrun (-n)
NODE$ snakemake -n




### feature notes ###

mupexi will ONLY run for mutations with filter "PASS". 
if you want to run for all samples you need to change 
the filter to “PASS”

The workflow splits the vcf according to number if mutations.
This is done since a large amount of mutations takes a very
long time to run with mupexi (>9H for 5000 mutations)
set "SPLIT_VCFS: false", to remove this behaviour

