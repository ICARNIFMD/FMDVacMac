
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FMDVacMac

<!-- badges: start -->
<!-- badges: end -->

The goal of FMDVacMac is to …

## Installation

You can install the development version of FMDVacMac from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("ICARNIFMD/FMDVacMac")
```

## Example

This is a basic example which shows you how to solve a common problem:

Users can either paste nucleotide sequences (for fewer isolates) or upload a FASTA file (.fasta or .txt) for larger datasets. Input must be in standard FASTA format, where each sequence starts with a definition line beginning with a carat (>) followed by a unique identifier. Sequences should contain only standard nucleotides: A, T, G, C. An example dataset of sequences is as follows:
>ICFMD182/2022
ACAACCTCCACAGGTGAGTCGGCTAATCCCGTGACTGCCACCGTTGAAAACTACGGAGGCGAGACACAGG
TCCAGAGACGCCAGCACACGGACGTCTCTTTCATATTGGACAGATTTGTAAAAGTGACACCAAAAGACCA
AATTAATGTACTGGACCTGATGCAAACCCCTGCTCACACTTTGGTGGGTGCCCTTCTCCGCACCGCCACC
TACTACTTCGCAGATTTAGAGGTGGCAGTGAAACACGAGGGGAACCTCACCTGGGTCCCGAACGGGGCGC
CCGAAAAAGCCTTGGACAACACCACCAATCCAACGGCTTACCACAAGGCACCACTCACCCGACTTGCACT
GCCGTACACGGCACCCCACCGTGTCATGGCTACTGTTTACAACGGGAACTGCAAGTACGGCGAGAGCCAC
GCAACTAGCGTGAGAGGTGACCTGCAAGTGTTGGCCCAGAAGGCAGCAAGGACGCTGCCTACCTCCTTCA
ACTACGGTGCCATCAAAGCTACTCGGGTGACTGAACTGCTTTATCGCATGAAGAGGGCTGAAACATACTG
CCCTCGGCCTCTTTTGGCCATCCACCCGAGTGAAGCTAGACACAAACAAAAGATTGTGGCACCTGTGAAG
CAACTTCTG
Machine Learning Method Selection:
FMDVacMac is a machine learning framework for predicting FMD vaccine matching scores (r1-values) from VP1 sequences. Users can select the best-performing model for prediction.
Sl. No.	Machine learning technique	Abbreviation
1	Support Vector Machine	         SVM
2	Random Forest	                 RFF
3	XGBoost	                         XGB

After the data is processed, FMDVacMac provides result in a table

Output: The output is provided in a tabular format which shows the predicted vaccine matching score (r1-value) of the Foot-and-Mouth Disease (FMD) virus isolate for two options: (i) when single method is selected; (ii) when all the methods are selected.
Case I: Single method (SVM/RF/XGB) selection. Suppose, SVM machine learning option was selected for predictive analysis. Then the output in following format is obtained.
Sl. No.	Isolate_ID	r ̂_1-value	Remarks
1	ICFMD182/2022	0.45	Protection
2	K154/12	        0.67	Protection
3	IND 168/2004	0.37	Protection
 
Interpretation of predicted r1-value:  As per OIE-recommended guidelines, the interpretation of the predicted r1-values (for VNT) is as follows:
The r ̂_1-value equal to or above 0.3 indicate a close antigenic relationship between the vaccine strain and the field isolate, i.e., it is likely that the vaccine strain will confer cross-protection against the field strain. Whereas, the r ̂_1-value less than 0.3 indicate a lack of such cross-protection. 

34.6 samb2 35.9 sams3 27.6 … … samg7 53.2

