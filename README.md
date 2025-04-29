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
devtools::install_github("ICARNIFMD/FMDVacMac")
```
## Example

This is a basic example which shows you how to solve a common problem:

The user provides two inputs: the field isolate and the vaccine strain, both in nucleotide sequence format (.fasta or .txt).Each sequence must start with a FASTA definition line (e.g., >Sequence_ID), followed by the nucleotide sequence in standard single-letter code (A, T, G, C).The field isolate is the FMD virus collected from infected animals during outbreaks, representing the circulating strain. Accurate identification is important because FMDV is highly variable, affecting vaccine effectiveness. The vaccine strain is a reference strain used in vaccines, selected for its ability to provide broad protection. However, due to the virus’s diversity and mutation rate, vaccine strains may not always offer effective protection against all field isolates.

field isolate:
>ICFMD182/2022
ACAACCTCCACAGGTGAGTCGGCTAATCCCGTGACTGCCACCGTTGAAAACTACGGAGGCGAGACACAGG
TCCAGAGACGCCAGCACACGGACGTCTCTTTCATATTGGACAGATTTGTAAAAGTGACACCAAAAGACCA
>K154/12
ACCACTGCGACGGGAGAGTCAGCAGACCCTGTTACCACCACCGTTGAGAACTACGGCGGTGAAACACAGG
TCCAAAGACGGCACCACACAAGTGTCGAGTTCGTCATGGACAGGTTTGTGAAACTGGAAGCTCCCAGCCC
>IND 168/2004
ACCACCACAACCGGTGAGTCGGCGGACCCGGTGACAACCACGGTTGAGAACTACGGAGGAGAAACTCAGA
CGGCCAGACGGCTCCACACTGACGTCGCCTTCGTTCTCGACAGGTTTGTGAAACTCACTGCACCCAAGAA

vaccine strains:
>VP1_R2/1975
ACCACCTCCCCGGGTGAGTCAGCTGACCCCGTGACCGCCACTGTTGAAAACTACGGCGGTGAGACACAGG
CCAACACACGGACGTCTCATTCATTTTGGACAGATTTGTAAAAGTGACGCCAAAAGACCAAATTAATGTA


Machine Learning Method Selection:
FMDVacMac is a machine learning framework for predicting FMD vaccine matching scores (r1-values) from VP1 nucleotide sequences. Users can select the best-performing model for prediction.
Sl. No.	Machine learning technique	Abbreviation
1	Support Vector Machine	         SVM
2	Random Forest	                 RFF
3	XGBoost	                         XGB

After the data is processed, FMDVacMac provides result in a table

Output: The output is provided in a tabular format, displaying the predicted vaccine matching score (r1-value) for the Foot-and-Mouth Disease (FMD) virus isolate.
The user can select the prediction method (SVM, RF, or XGB). For example, if the SVM (Support Vector Machine) option is selected for the predictive analysis, the output will be presented in the following format:
Sl. No.	Isolate_ID	r ̂_1-value	Remarks
1	ICFMD182/2022	0.45	      Protection
2	K154/12	        0.67	      Protection
3	IND 168/2004	0.37	      Protection
 
Interpretation of predicted r1-value:  As per OIE-recommended guidelines, the interpretation of the predicted r1-values (for VNT) is as follows:
The r ̂_1-value equal to or above 0.3 indicate a close antigenic relationship between the vaccine strain and the field isolate, i.e., it is likely that the vaccine strain will confer cross-protection against the field isolate. Whereas, the r ̂_1-value less than 0.3 indicate a lack of such cross-protection. 

34.6 samb2 35.9 sams3 27.6 … … samg7 53.2
