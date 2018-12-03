mitoimp
====

## Description
An imputation tool for human mitochondrial genome based on shared allele pairwise distances.

## Requirement
・python 2.7 (pypy2)

・[MAFFT](https://mafft.cbrc.jp/alignment/software/): Multiple alignment program for amino acid or nucleotide sequences

## Usage

**Basic Usage**

`mitoimp.py -i input.fasta [-k 5] [-f 0.7] [-t 4]`

**optional arguments:**

| Flag | Description | File Format, Parameter etc. |
|:-----------|:------------|:------------|
| **-i**       | query sequence (**required**) | Single-FASTA format |
| **-p**       | in-house (customized) panel sequences | Multi-FASTA format |
| **-w**       | window-size  (default: 16569)           | 1 〜 16569  |
| **-k**       | k-number  (default: 5)     |1 〜 max of panel sequences  |
| **-f**       | the threshold frequency to determine a genotype  (default: 0.7)  | 0.5 〜 1.0 |
| **-t**       | multiprocessing numbers (detault: the max of available CPU-threads) | 1 〜 (max: -1) |
| **-no_aln**  | set a switch to non-alignment mode  (default: Disable)  |  |
| **-v**       | show program's version number and exit  | |
| **-h**       | show this help message and exit         | |

## Version

1.0.0 (beta)

## Licence

[MIT] https://github.com/omics-tools/mitoimp/blob/master/LICENSE

## Author

Koji Ishiya