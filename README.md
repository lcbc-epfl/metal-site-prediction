# Metal3D and Metal1D: Accurate prediction of transition metal ion location via deep learning

[![DOI](https://zenodo.org/badge/456988168.svg)](https://zenodo.org/badge/latestdoi/456988168)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/lcbc-epfl/metal-site-prediction/blob/main/Metal3D/ColabMetal.ipynb)

If using this work please cite:

>Accurate prediction of transition metal ion location via deep learning
>S.L. Dürr, A. Levy, U. Rothlisberger
>bioRxiv 2022.XXXXXXXX; doi: https://doi.org/XXXXXXXXXXXX


# How to run predictions

No installation, no account required use [Metal3D on Huggingface Spaces](https://hf.space/simonduerr/metal3d)

If you prefer a notebook based environment [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/lcbc-epfl/metal-site-prediction/blob/main/Metal3D/ColabMetal.ipynb)

Command line usage is described below

# How to install and run locally

For local installation run the following commands to setup the environment. 

```
conda env create -f environment.yml
conda activate metalprediction
cd Metal3D
```
You need to have VMD installed to view predictions directly from the commandline program (connect with ssh -X if working on a remote machine), download VMD from [uiuc.edu](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD). 
Alternatively you can use `--writecube --cubefile nameofcube.cube --softexit` and view the predicted maps in UCSF Chimera or any other viewer that supports cube file. 

**Typical commands would be:** 

Analyze all ASP, CYS, ASN, GLN, GLU and HIS residues, write a pdb file with the found probes and write the maximum probabilty to a text file. Will open VMD viewer.

`./metal3d.py --pdb PDB.pdb --metalbinding --writeprobes --probefile metalsites.pdb --maxp `

Analyze only specific residues in the pdb file and write a cubefile to disk without openening VMD. 
`./metal3d.py --pdb PDB.pdb --id 91 94 116 --writecube --cubefile test.cube --softext`

Display all possible options
`./metal3d.py --help`

# Data

The PDB codes used for training, validation and testing are available in `data`. 
The PDB codes used for the selectivity analysis including the residue ids of the coordinating residues are available in `data` as `selectivity_analysis_sites.csv`. 


# License

All code is licensed under MIT license, the weights of the network are licensed under CC BY 4.0.
