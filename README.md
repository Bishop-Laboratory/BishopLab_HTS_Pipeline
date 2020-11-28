# BishopLab HTS Pipeline
This pipeline is used to automate general HTS analysis tasks in the Bishop Lab. 

Usage info:

```
###############################################################################################################################
  ______ _     _                 _           _         _   _ _____ _____    ______ _            _ _
  | ___ (_)   | |               | |         | |       | | | |_   _/  ___|   | ___ (_)          | (_)
  | |_/ /_ ___| |__   ___  _ __ | |     __ _| |__     | |_| | | | \ `--.    | |_/ /_ _ __   ___| |_ _ __   ___
  | ___ \ / __| '_ \ / _ \| '_ \| |    / _` | '_ \    |  _  | | |  `--. \   |  __/| | '_ \ / _ \ | | '_ \ / _ \
  | |_/ / \__ \ | | | (_) | |_) | |___| (_| | |_) |   | | | | | | /\__/ /   | |   | | |_) |  __/ | | | | |  __/
  \____/|_|___/_| |_|\___/| .__/\_____/\__,_|_.__/    \_| |_/ \_/ \____/    \_|   |_| .__/ \___|_|_|_| |_|\___|
                          | |                                                       | |
                          |_|                                                       |_|
###############################################################################################################################

BishopLab_HTS_Pipeline [-m mode] [-a SRP_accession] [-p project_directory] [-g genome_directory] [-P num_threads]

General options:

  -m|--mode              mode      Choose analysis mode ('genomeBuild', 'RNASeq', 'HiCSeq', or 'ChIPSeq')
  -a|--accession         SRA       SRA accession to download data from. [e.g. SRP045672].
  -l|--accessionList     file      File with newline-separated list of SRA accessions.
  -c|--python27Env       env       If in ChIPSeq mode or running TEcount, name of conda env with python 2.7.
  -p|--projectDir        dir       Project directory. [default = 'RNASeq_SRA_Pipeline_Project']
  -g|--genomeDir         dir       Genome directory. [Not required if genome built with '-m genomeBuild']
  -P|--num_threads       int       Specify number of threads. [default = 1]
  --synapseID            synID     Parent of desired synapse file directory.
  --synapseIDList        file      File with newline-separated list of synpaseIDs (samples or sample folders).
  --synUserName          uName     Synapse username for downloading data (if --synapseID is specified).
  --runInfoFile          file      SRA run info table file (optional).
  --noFastp                        Do not use fastp to perform adapter trimming, filtering, and fastq QC.
  --keepFastq                      Do not delete fastq files once processing is finished.
  --downloadOnly                   Will only download fastqs, merge tech replicates + do fastp (--mode not required).
  --returnBamsOnly                 Returns bams, coverage tracks, and normalized signal tracks.
  --noMerge                        Do not attempt to merge technical replicates from SRA.
  --help                           Display usage info

genomeBuild --mode options:

  -s|--species           str       Select 'human', 'mouse', or 'both' [default = 'both'].
  -t|--toBuild           str       's' [Salmon] 'S' [STAR] 'f' [STAR-Fusion] 'b' [BWA] Default: 'sSfb'

RNASeq --mode options:

  -r|--testsRNASeq       string    's' [salmon] 't' [TECount] 'S' [Splicing-STAR] 'f' [STAR-Fusion] Default: 'stSf'
  --starFusionEnv        env       Conda ENV with STAR-Fusion and STAR 2.7.0 [Due to conflict with 2.7.1].
  --force_spliceSE                 Not recommended: If 'S' in tests, single-end reads are run with splicing STAR.
  --force_fusionSE                 Not recommended: If 'f' in tests, single-end reads are run with STAR-fusion.

ChIPSeq --mode options:

  -G|--groupFile         file      TSV with columns 'SampleName', 'Group', 'Replicate', and 'Condition'. Can include 'PeakType'.

```

