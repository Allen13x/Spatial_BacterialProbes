# Spatial Visium 10x script for analysis with Bacterial Probes

Script and tools used for the study:

_Probe-based spatial host-pathogen genes expression to study bacterial pathogenesis and the regulation of bacterial virulence factors in tissue_

Data available at:
GSE242471

## Step 1: Desgin Probes

Starting from the sequence of the interested target gene (Gene.fasta), paste the sequence into the Kmer.py script
This will generate a list for all the 50mer present in the sequence into the file Gene50mer

Use the script_1.sh to generate a fasta file containing all 50mers compliant with the 10x specifics

```
./scrpti_1.sh Gene
```

Use the obtained Gene_50mer.fasta file on ncbi blast to check for specificity and sensibility of the probes and save the good probes ID into a Gene_goodhit.txt file

Use the script_2.sh to obtain the final list of possible probes to use for the experiment

```
./script_2.sh Gene
```
The ouptut will be put into the OUTPUT folder as Gene.tab


## Step 2: Spatial Visium10x preprocessing 

Follow the 10x instructions to build you probe dataset and reference genome accordingly your new probes
Follow Cellranger and 10x tutorial fro obtaining the preprocessed dataset


## Step 3: Spatial data analysis using Seurat

Look at the SpatialProbes.R script






