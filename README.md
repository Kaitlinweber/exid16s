# exid16s

**Ex**traction and **Id**entification of **16S** rDNA (ExId16S) 

## General information 

ExId16S is a package in which an extraction and identification with 16S rDNA data is performed on assembled Whole Genome Sequence data. The tool Barrnap is used for the extraction of 16S rDNA, this tool localizes the ribosomal DNA sequences using Hidden Markov Models (HMM).
From the output with Barrnap, the 16S rDNA sequence is used for as Kraken2 input. Kraken2 identifies bacteria using an alignment free algorithm. 

![exid16s (1)](https://user-images.githubusercontent.com/64156013/163177370-bdbd06fa-fa33-435f-8234-b3192b60cd41.jpg)

## Installation

1. Clone the repository:

```
git clone https://github.com/Kaitlinweber/exid16s.git
```

2. 
