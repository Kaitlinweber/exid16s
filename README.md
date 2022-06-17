# ExId16S

**Ex**traction and **Id**entification of **16S** rDNA (ExId16S) 

## Description

ExId16S is a package in which an extraction and identification with 16S rDNA data is performed on assembled Whole Genome Sequence data. The tool Barrnap is used for the extraction of 16S rDNA, this tool localizes the ribosomal DNA sequences using Hidden Markov Models (HMM).
From the output with Barrnap, the 16S rDNA sequence is used for as Kraken2 input. Kraken2 identifies bacteria using an alignment free algorithm. 

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


