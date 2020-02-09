# FBMN_Evaluator

Feature-Based Molecular Networking (FBMN) is a powerful unsupervised approach for organizing
LCMS data bases on retention time, MS1, and MS2 over multiple samples.

https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/

Louis Felix Nothias, et. al.
bioRxiv 812404; doi: https://doi.org/10.1101/812404

A key challenge in generating FBMN's is determining appropriate processing settings when MS1
features are generated, using MZmine or other tools.  If merge settings are too permissive,
different molecular species will be inappropriately merged, whereas if merge settings are
too restrictive, the same molecular species will be inappropriately split.  The first case
artificially hides molecular diversity in the sample, whereas the second overestimates it.
The second case is especially a challenge for downstream visualization (e.g. 'ili).

Here a simple command line Python script is presented for comparing the relative quality of
seveal FBMN's is presented.  This script can also be run from a Jupyter notebook.  

FBMN's can be evaluated as a visual plot of retention timeand m/z features assigned via MS/MS.  

As well, figures such as the total number of FBMN features, and FBMN features within a 
certain retention time and m/z tolerance are reported.

Usage:

1. Run jobs via FMBN workflow.
2. After job is complete, click on "done" then "Download Cytoscape Data"
3. Unzip the network, and in the "quantification_table" directory there is your network!
4. Save one or multiple networks to a known directory
5a. Go to appropriate command line interface: terminal (Mac/Linux) or Windows Power Shell
5b. Open appropriate Jupyter Notebook and run from cell.
6. Run command as below, substituting uppercase text:

python FBMN_Evaluator.py --path PATH_TO_NETWORK_DIR --mztol MZ_TOLERANCE_IN_PPM --rttol
RETENTION_TIME_TOL_IN_SEC

Output:

Path to plots, number of FBMN features, number of FBMN features within mz/rt tolerance.
