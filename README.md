# FetalThyroidAtlas
This repository contains the code needed to reproduce the custom analyses used in the publication *A developmental cell atlas of the human thyroid gland*. We recommend reading the Methods section of the paper before using this code.

bioRxiv paper: [https://www.biorxiv.org/content/10.1101/2024.08.22.609152v1](https://www.biorxiv.org/content/10.1101/2024.08.22.609152v1)

**Note to users**  
  
In some of the scripts, absolute paths are used when loading data. These have been kept intentionally for legacy and reproducibility reasons, as they reflect the original analysis environment used for the paper. Please feel free to modify these paths to point to the input files generated on your own system when running the code.  
  
We also provide a set of helper functions that are called across all `.py` files, along with the colour definitions used for the different clusters throughout the analysis. These shared utilities are collected in the *snippet_and_more* folder.  
  
We hope this structure makes it easier to reproduce, adapt, and extend the analysis.
