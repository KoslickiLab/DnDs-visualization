# Hierarchical Edge Bundling of Genomic dN/dS

Hierarchical edge bundling is an easy way to display quantitative pairwise relationships among thousands of genomes [1]. Here, we use the Hierarchical Edge Bundling R Package to present genomic dN/dS estimations between archaeal and bacterial genomes from the Genome Taxonomy Database, where red and blue represent positive and negative selection, respectively. Please find the reproducible script to generate the Hierarchical Edge Bundling Figure presented in manuscript for "Leveraging FracMinHash Containment for Genomic dN/dS". 

![Figure](./Hierarchical_Edge_Bundling_tree/GTDB/output_100samples.png)

# Table of contents

- [Environment setupe](#Environment-Setup)
- [Dataset Information](#Dataset-Information)
- [Figure Generation](#Figure-Generation)
- [Acknowledgements](#Acknowledgements)
- [References](#References)
- [Our paper](#Our-paper)

# Environment Setup

Create either a conda or pip environment to reproduce this figure.

## Conda environment

```
conda env create -f environment.yml
conda activate hierarchical_env
```

## Pip environment

```
python3 -m venv hierarchical_pip_env
source hierarchical_pip_env/bin/activate
pip install -r requirements.txt
```


# Dataset Information

This file is to help you out on what the files mean and there purpose. If you have any questions, feel free to reach out!

**Filename:** marinus_taxonomy.csv

**Purpose:** This file contains taxonomic information that can be used to produle a phylogenetic tree. We should be using the genus information as the cut off.
i


| Column Name | Description |
|---|---|
| NCBI_fasta_name | the identifier also known as an accession number used by ncbi.com |
| custom_fasta_name | a name Judith gave to the sequence based on published data |
| phylum | taxonomic classifier |
| class | taxonomix classifier |
| order | taxonomix classifier |
| family | taxonomix classifier |
| genus | taxonomix classifier that we would potentially use for tree generation |
| species | taxonomix classifier |
| strain | taxonomix classifier |


**Filename:** dnds_constant_15.csv 

**Purpose:** This file contains the dN/dS values. dN/dS above 1 is postive selection (generally shown in red) and a dN/dS below 1 is negative selection (generally shown in blue). The most important columns here are A, B, and dNdS_ratio_constant.


| Column Name | Description |
|---|---|
| A | genome A |
| B | genome B |
| containment_nt | containment estimated among nucleotide sequences of genomes A and B |
| ksize | the size of the subsets used for sketching |
| containment_protein | containment estimated among protein sequences of genomes A and B |
| sequence_comparison | comparison between genome name A and B |
| dN_constant | estimated nonsynonymous mutations with constant |
| dS_constant | estimated synonymous mutations with constant |
| dN | estimated nonsynonymous mutations |
| dS | estimated synonymous mutations |
| dNdS_ratio | the estimated dN/dS without constant |
| dNdS_ratio_constant | the estimated dN/dS with constant |

# Figure Generation

Please follow and run jupyter notebook instructions here: [DnDs-visualization/Hierarchical_Edge_Bundling_tree/GTDB/test-code_GTDB copy_jzr_modify.ipynb](https://github.com/KoslickiLab/DnDs-visualization/blob/main/Hierarchical_Edge_Bundling_tree/GTDB/test-code_GTDB%20copy_jzr_modify.ipynb)

<!-- ![Figure](./Hierarchical_Edge_Bundling_tree/GTDB/output_100samples.png) -->

# Acknowledgements

Thank you Jinglin Feng for your contributions to this figure!

# References

[1] D. Holten (2006). Hierarchical Edge Bundles: Visualization of Adjacency Relations in Hierarchical Data. IEEE Transactions on Visualization and Computer Graphics, 12(5)

# Our paper

Leveraging FracMinHash Containment for Genomic dN /dS. Judith S. Rodriguez,
Mahmudur Rahman Hera, and David Koslicki. In preparation.

