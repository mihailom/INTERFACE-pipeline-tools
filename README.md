# INTERFACE-pipeline-tools
Codes to analyze paired reads from INTERFACE assays (in .bed) format to output the weighted average for each asRNA, normalized between 0 and 1 within unique target RNAs. Note that the INTERFACE assay quantifies the hybridization efficacy of a library of user-defined asRNAs (targeting distinct RNA regions in vivo) via a transcriptional elongation mechanism. Codes written by mathematics intern Daniel Herrera under supervision of grad student Mia Mihailovic
1. First, use [matrix_creator_ends.py](https://github.com/mihailom/INTERFACE-pipeline-tools/blob/master/matrix_creator_ends.py) to process each .bed file into a matrix containing frequency of ends for each nucleotide position in the range (1,max)
2. Then, use [final_frequency_analysis.py](https://github.com/mihailom/INTERFACE-pipeline-tools/blob/master/final_frequency_analysis.py) with your synthetic INTERFACE genome for statistical analysis on each set of end frequency matrices corresponding to distinct replicates within a single folder




