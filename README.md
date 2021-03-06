MitoIMP
====

## Description
MitoIMP: A computational framework for imputation of missing data in low-coverage human mitochondrial genome

## Requirement
・python 2.7 (pypy2)

・[MAFFT](https://mafft.cbrc.jp/alignment/software/): Multiple alignment program for amino acid or nucleotide sequences

  　Please refer to [the official site](https://mafft.cbrc.jp/alignment/software/) for details of the installation of MAFFT.

## Installation

Clone this repository into your local machine  
`git clone https://github.com/omics-tools/mitoimp.git`  

## Usage

**Basic Usage**

`mitoimp.py -i input.fasta [-k 5] [-f 0.7] [-t 4]`

**Example1**

An imputed sequence and a summary table are output to the same directory as the input file.

`mitoimp.py -i ./sample_data/A_cov65.fasta -k 5 -f 0.7 -t 4`

**Example2**

If you want to use a sequence oriented to the rCRS position by other alignment software, please set -no_aln flag.

`mitoimp.py -i ./sample_data/Z_cov85.fasta -k 5 -f 0.7 -t 4 -no_aln`

**optional arguments:**

| Flag | Description | File Format, Parameter etc. |
|:-----------|:------------|:------------|
| **-i**       | query sequence (**required**) | Single-FASTA format |
| **-p**       | in-house (customized) panel sequences | Multi-FASTA format |
| **-w**       | window-size  (default: 16569)           | 1 〜 16569  |
| **-k**       | k-number  (default: 5)     |1 〜 max of panel sequences  |
| **-f**       | the threshold frequency to determine a genotype  (default: 0.7)  | 0.5 〜 1.0 |
| **-t**       | multiprocessing numbers (default: the max of available CPU-threads) | 1 〜 (max: -1) |
| **-no_aln**  | set a switch to non-alignment mode  (default: Disable)  |  |
| **-v**       | show program's version number and exit  | |
| **-h**       | show this help message and exit         | |

## Version

1.0.1 (beta)

## Licence

[MIT] https://github.com/omics-tools/mitoimp/blob/master/LICENSE

## Citation

Ishiya, K., Mizuno, F., Wang, L., & Ueda, S. (2019). [MitoIMP: A Computational Framework for Imputation of Missing Data in Low-Coverage Human Mitochondrial Genome. Bioinformatics and biology insights.](https://journals.sagepub.com/doi/10.1177/1177932219873884)
