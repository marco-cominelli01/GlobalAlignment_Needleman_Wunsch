# Parallel implementation of the Needleman-Wunsch algorithm for Global Sequence Alignment

This repository contains all the necessary files and documentation for the project.

## Full Documentation

The complete README, including detailed descriptions, images, and other resources, is available in PDF format. 
You can view or download it by clicking the link below:

[View the complete README in PDF format](./README.pdf)

## Project Overview

This program computes in parallel the optimal global alignment of two sequences (DNA, RNA or protein sequences).

### Usage
GlobalAlignment_Needleman_Wunsch.py [-h] (-s SEQUENCE SEQUENCE | -i INPUT INPUT) [-t {d,r,p}] [-g GAP] [-mm MISMATCH] [-m MATCH] [-b BLOSUM] 
[-o OUTPUT] [-q] [-c CORES] [-v {1,2}]

### Options
  -h, --help           							                    show this help message and exit\
  -s SEQUENCE SEQUENCE, --sequence SEQUENCE SEQUENCE		the two sequences to be aligned\
  -f FASTA FASTA, --fasta FASTA FASTA					          the two FASTA files containing the sequences to be aligned\                
  -t {d,r,p}, --type {d,r,p}							              type of sequences to align: 'd' for DNA sequences, 'r' for RNA sequences and 'p' for protein\ 
                                                        sequences. Default is 'd'\
  -g GAP, --gap GAP     							                  negative GAP penalty. Default is -4\
  -mm MISMATCH, --mismatch MISMATCH				              negative MISMATCH penalty. Default is -5\
  -m MATCH, --match MATCH						                    positive MATCH score. Default is 5\
  -b BLOSUM, --blosum BLOSUM					                  use the specified BLOSUM matrix for MATCHES/MISMATCHES (compatible only with protein sequences)\
  -o OUTPUT, --output OUTPUT						                save the alignment(s) in the specified output file\
  -q, --quiet           							                  don't display the output\
  -c CORES, --cores CORES						                    number of cores to use. Default is 3 (if available)\
  -v {1,2}, --verbose {1,2}						                  increase verbosity. Default is 1\
  -a, --approximation							                      don't show all possible alignments but a reduced number of them\

For more detailed information, please refer to the PDF document linked above.
