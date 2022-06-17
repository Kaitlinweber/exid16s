# ExId16S

**Ex**traction and **Id**entification of **16S** rDNA (ExId16S) 

## Description

ExId16S is a package in which an extraction and bacterial identification with 16S rDNA data is performed on assembled Whole Genome Sequence data. The tool [Barrnap](https://github.com/tseemann/barrnap) is used for the extraction of 16S rDNA, this tool localizes the ribosomal DNA sequences using Hidden Markov Models (HMM). From the output with Barrnap, the 16S rDNA sequence is used for bacterial identification with [Kraken 2](https://github.com/DerrickWood/kraken2input), which identifies bacteria using an alignment free algorithm. For bacterial identification a 16S rDNA database is required, this database must be compatible with Kraken 2.

<p align="center">
    <img src= "https://user-images.githubusercontent.com/64156013/174275271-26740583-4cb5-48ea-b590-e92656fe7eb3.png", alt="flowchart ExId16S">
</p>


## Installation

1. Clone the repository:

```
git clone https://github.com/Kaitlinweber/exid16s.git
```

2. Enter directory and install conda enviroment:

```
cd exid16s
conda env install -f envs/exid16s.yaml
```

3. Install python package 

```
pip install .
```

## Required parameters 
* ```-i, --input``` Pathway to all the input FASTA files. If you give a pathway to a directory, it is important you select all the files with (/*)
* ```-o, --output``` Directory pathway, in which the output subdirectories (barrnap_result, FASTA_16S_sequence and Kraken2_kreports) will be created. If the output directory does not excists, it will be created.
* ```-db, --database``` Directory for Kraken 2 compatible database. This database needs to contain at least 3 files: hash.k2d, opts.k2d, taxo.k2d. 


## Usage - Basic command to run the package

```
exid16s -i [path/to/input/dir/*] -o [path/to/output/dir] -db [path/to/database/dir]
```


## Output explanation 
* **barrnap_result:** In this directory the barrnap multi-FASTA files will be stored. This are the files with the 5S 23S and 16S sequences.
* **FASTA_16S_sequence:** Directory with the extracted 16S sequences from Barrnap output, are redirected to FASTA or multi-FASTA files based on the number of copies. 
* **Kraken2_kreports:** This directory contains one kreport (Kraken 2 identification output) per sample, and a combined summary file with the genus and species information of each kreport. 

