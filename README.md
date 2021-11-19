# ripp-design
Jupyter notebooks and associated data for designing enzyme-modified peptides

# To install:
All of the .py files in the packages folder along with the .py files in the ripp-analysis (github.com/VoigtLab/ripp-analysis) packages folder need to be downloaded and accessible either in the path or in the working directory (they need to be importable in python).


## In particular, files needed are:

- modification_rules.py
- ripp_design.py
- enzymeanalysis.py
- enzymeplots.py
- lcms.py
- lcmsanalysis.py

## The following additional packages (versions) were used:

- Code was written, run, and tested with Python 3.7.10

- matplotlib (tested with version 3.0.2)
- pandas (tested with version 0.23.4)
- numpy (tested with version 1.19.4)
- seaborn (tested with version 0.9.0)
- regex (tested with version 2.5.77)
- biopython (tested with version 1.79)
- scipy (tested with version 1.5.4)


# To Use:
Ipython notebooks are the best way to see how the software is used and serve as a "How To".

## Folder "analysis"
Includes ipython notebooks and related metadata for analyzing LC-MS data of peptide variants to calculate enzyme-peptide constraints on modification.

extracts.xlsx is an excel spreadsheet that contains all metadata for the extracts, including the peptide and modifying enzyme present in the extract, the expected mass, the expected mass shift after modification, etc. This file is required to parse the raw LC-MS data.

### Notebooks are split up by goal:

- Import extracts to dataframe.ipynb -- This must be run first, analyzes raw LC-MS data to export a processed dataframe that is used by the other notebooks. With the full dataset, this takes about 3 hours on our server running 20 threads. The final result is a pickled pandas dataframe. The folder 'extract_dataframes' must be present in the working directory in order to run this notebook. The 'extracts.xlsx' spreadsheet must also be present.                                              
- Leader varaint analysis (Figure 1, SI Notes).ipynb -- This notebook details the process used to generate plots shown in Figure 1 and Supplementary Notes. It pulls data from 'dataset.pickle'
- Core variant analysis (Figure 2, SI Notes).ipynb -- This notebook details the process used to generate plots shown in Figure 2 and Supplementary Notes. It pulls data from 'dataset.pickle'
- SI Figures 2-5, 8-12 (Mod Validation).ipynb -- This notebook details the process used to generate plots shown in Supplementary Figures 2-5 and 8-12. It pulls data from 'dataset.pickle'
- Raw Chromatograms (SI Figure 6).ipynb -- This notebook details the process used to generate plots shown in Supplementary Figure 6. It pulls data from 'dataset.pickle'
                                                
## Folder "design"
Includes ipython notebook for designing new RiPPs based on peptide constraints. The example detailed in the notebook is the same as what is shown in Figure 3 of the manuscript.

## Folder "tgn-sample-dataset"
Contains all the raw data, extract.xlsx file, and notebook for analyzing and exporting plots for just the enzyme TgnB. It is meant to be an example dataset to interact with the software. To use - open the ipython notebook in the folder, ensure that all packages are importable and dependencies installed, and execute the code.
  
## Folder "lcms"
Contains all the raw data acquired/used for this study, split between folders "qqq" and "qtof". These folders are required if running the "Import extracts to dataframe" notebook, and the paths must be updated in that notebook
