Two subdirectories:

`Models/` : all python and R codes to run the SVZ, SVRZ, SIVZ, and SIVRZ models simulations, the parameter optimization to target concentration and histograms of measured abundances from Carlson *et al.* 2022 (https://www.nature.com/articles/s41564-022-01088-x).    
Generates pdfs to create Figure 2 to 6, Figure S5 to S17 and .txt file with data from Table 1 and 2

`trait_data/` : all R codes to generate life history trait models.  
Generates pdfs to create Figure S1 to S4 and some data used in `Models/`

In each subdirectory is a `README.md` describing precisely how to run the scripts, the pdfs it generates and where the panel of each figure in the paper can be found.

The scripts from the two subdirectory can be run independently from each other (the scripts in `Models/` take some data generated in `trait_data/` as input but they are already pre-generated so each subdirectory is independent)

R version 4.3.2 and Python 3.12.9 were used
