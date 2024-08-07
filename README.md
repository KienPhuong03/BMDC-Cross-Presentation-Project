# BMDC-Cross-Presentation-Project
This repository contains the workflow my mentor and I made for our project as well as the proteomics data.
Everything is grouped by experiments, with data and analyses being in the same file. A powerpoint presentation will also be provided
as a summary of the analyses made.

Usage:
The workflow has 2 parts:
1. Unlabeled MHC Filter:
The parameters are:
--path: the path to the folder containing the data files
--i: name of the input file
--o: name of output file
--type: "Mascot" for unlabeled (global proteomics), for MHC use "Sequest"
2. Workflow file:
--path: the path to filtered csv file
--data: name of input csv files
--sample_names: corresponding list of names of samples
--label: list of labeled amino acids
--sep: analyze each amino acid separately
--not-sep: analyze all labels in aggregate
--o: output file name

