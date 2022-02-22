---
author-meta:
- "Simon L. D\xFCrr"
- Andrea Levy
- Ursula Rothlisberger
bibliography: []
date-meta: '2022-02-14'
header-includes: "<!--\nManubot generated metadata rendered from header-includes-template.html.\nSuggest improvements at https://github.com/manubot/manubot/blob/master/manubot/process/header-includes-template.html\n-->\n<meta name=\"dc.format\" content=\"text/html\" />\n<meta name=\"dc.title\" content=\"Accurate prediction of metalsites using deep learning\" />\n<meta name=\"citation_title\" content=\"Accurate prediction of metalsites using deep learning\" />\n<meta property=\"og:title\" content=\"Accurate prediction of metalsites using deep learning\" />\n<meta property=\"twitter:title\" content=\"Accurate prediction of metalsites using deep learning\" />\n<meta name=\"dc.date\" content=\"2022-02-14\" />\n<meta name=\"citation_publication_date\" content=\"2022-02-14\" />\n<meta name=\"dc.language\" content=\"en-US\" />\n<meta name=\"citation_language\" content=\"en-US\" />\n<meta name=\"dc.relation.ispartof\" content=\"Manubot\" />\n<meta name=\"dc.publisher\" content=\"Manubot\" />\n<meta name=\"citation_journal_title\" content=\"Manubot\" />\n<meta name=\"citation_technical_report_institution\" content=\"Manubot\" />\n<meta name=\"citation_author\" content=\"Simon L. D\xFCrr\" />\n<meta name=\"citation_author_institution\" content=\"Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland\" />\n<meta name=\"citation_author_orcid\" content=\"0000-0002-4304-8106\" />\n<meta name=\"twitter:creator\" content=\"@simonduerr\" />\n<meta name=\"citation_author\" content=\"Andrea Levy\" />\n<meta name=\"citation_author_institution\" content=\"Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland\" />\n<meta name=\"citation_author_orcid\" content=\"0000-0003-1255-859X\" />\n<meta name=\"citation_author\" content=\"Ursula Rothlisberger\" />\n<meta name=\"citation_author_institution\" content=\"Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland\" />\n<meta name=\"citation_author_orcid\" content=\"0000\u20100002\u20101704\u20108591\" />\n<meta property=\"og:type\" content=\"article\" />\n<meta property=\"twitter:card\" content=\"summary_large_image\" />\n<link rel=\"icon\" type=\"image/png\" sizes=\"192x192\" href=\"https://manubot.org/favicon-192x192.png\" />\n<link rel=\"mask-icon\" href=\"https://manubot.org/safari-pinned-tab.svg\" color=\"#ad1457\" />\n<meta name=\"theme-color\" content=\"#ad1457\" />\n<!-- end Manubot generated metadata -->"
keywords:
- metal
- protein engineering
- deep learning
lang: en-US
manubot-clear-requests-cache: false
manubot-output-bibliography: output/references.json
manubot-output-citekeys: output/citations.tsv
manubot-requests-cache-path: ci/cache/requests-cache
title: Accurate prediction of metalsites using deep learning
...





<!-- 
<small><em>
This manuscript
was automatically generated
on February 14, 2022.
</em></small> -->

## Authors

<!-- 

+ **Simon L. Dürr** * <br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [0000-0002-4304-8106](https://orcid.org/0000-0002-4304-8106)
    · ![Twitter icon](images/twitter.svg){.inline_icon}
    [simonduerr](https://twitter.com/simonduerr)<br>
  <small>
     Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland
  </small>

+ **Andrea Levy** * <br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [0000-0003-1255-859X](https://orcid.org/0000-0003-1255-859X)<br>
  <small>
     Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland
  </small>

+ **Ursula Rothlisberger** † <br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [0000‐0002‐1704‐8591](https://orcid.org/0000‐0002‐1704‐8591)<br>
  <small>
     Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland
  </small>
 -->



[![ORCID icon](images/orcid.svg){height="11px" width="11px"}](https://orcid.org/0000-0002-4304-8106)
Simon L. Dürr<sup>1,*</sup>,
[![ORCID icon](images/orcid.svg){height="11px" width="11px"}](https://orcid.org/0000-0003-1255-859X)
Andrea Levy<sup>1,*</sup>,
[![ORCID icon](images/orcid.svg){height="11px" width="11px"}](https://orcid.org/0000‐0002‐1704‐8591)
Ursula Rothlisberger<sup>1,†</sup>


<!--<sup>☯</sup> --- These authors contributed equally. <br> -->
<sup>†</sup> --- To whom correspondence should be addressed: ursula.roethlisberger@epfl.ch
<small>

###### Affiliations

1. Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland




### Abstract

## Introduction

## Methods

## Results

### Metal1D


### Metal3D

The model was trained using a train/test/val split based on sequence identity. We used the MMSeqs2 clustered PDB at 30% sequence identity and used the highest resolution structure from each cluster that contained a zinc, did not contain DNA and had resolution <2.5A. In case no structure was found the cluster was discarded. 

The training examples were sampled from the chosen structure by choosing a balanced number of boxes from each protein that contained or did not contain a zinc. Each box was randomly rotated such that the model is insensitive to rotation. 

Metal3D predicts a per residue score that can then be averaged over all residues or used individually (e.g for protein design). 

Hyperparameter tuning using Ray tune


Evaluation of quality of predictions per box using discretized Jaccard Score (similarity of two sets). We noticed that at the edges often spurious density is predicted . 


![Discretized Jaccard index using different edge cutoffs and different probability cutoff](images/jaccard.png){#fig:jaccard}



#### Selectivity for other metals

<!-- Probably no figure needed  -->

- Trained exclusively on zinc, predictions are similar for metals with different binding modes e.g Copper
- For sodium (binds via backbone O) not predicted
- Some sites for Mg are predicted with high but not super high probability

Figure 4 panels Copper protein, False positive Calcium Zinc and Zinc->Calcium site, Dimetallic zinc site

Discussion:
Lack in selectivity could be related to smoothing the gaussian quite a bit when training (anything >0.05) is a hit. 

Resolution of grid might be an issue 
Might be improved by improving the grid resolution to 0.5 A

#### Alpha Fold 

AlphaFold often predicts sidechains in metal ion binding sites in the holo conformation. Services like AlphaFill use homology to transplant metals from similar PDB structures to the AlphaFold structure. Metal3D does not use homology can even deal with metal sites that are only partially in the holo conformation. 

![Metal3D, AlphaFill and 3RZV zinc positions. Metal3D places the metal with high accuracy even if coordination is not perfectly predicted. Probability map colored at p=0.6](images/3rzv_alphafold_metal3d_alphafill.png)

#### Hidden/transient metalsites

Metallodrugs are in important class of drugs that rely on binding inhibitors to a protein (or DNA). Metal3D can be used to screen the hidden metalloproteome by finding transient metal ion binding sites. 

The site where Rapta binds is detected with p=0.3 but in a high resolution structure without (1KX4) there is a salt bridge with a lysine that might occlude metal detection. One could weight by the rotamer/do MD simulation e.g similar as for cryptic pockets. 

![RAPTA (Ru-arene-phosphaadamantane) inhibitor bound to Nucleosome core particle. Metal3D identifies the two Ruthenium sites with low probabilities.](images/5xf6_rapta_vis.png)

#### HCA2: case study

HCA2 is the first enzyme where a catalytic zinc was discovered and is therefore one of the best studied metalloenzymes to date with a rich throve of mutational data available. 
On the wildtype enzyme crystal structure (2CBA) Metal3D perfectly recapitulates the location of the active site metal when using a high probability cutoff (p>0.4). The sites predicted with lower probability all look like reasonable transient binding sites at the surface of the protein. 

![Probability evolution in HCA2 for different probability cutoffs A) p=0.1 B) p=0.2 C) p=0.3 D) p=0.4 E) p= 0.75](images/2CBA_probabilities.png)

We used in silico generated mutants matching mutants in the first and second shell of the active site zinc and probed the effect on the predicted metal probability. For mutants that decrease zinc binding also a drop in probability can be observed that correlates well with the experimentally measured $K_d$. We used a consistent set of $K_d$ values from the literature.
![Probability evolution in HCA2 for different probability cutoffs A) p=0.1 B) p=0.2 C) p=0.3 D) p=0.4 E) p= 0.75](images/kd_vs_p_nolog_newmethod.jpg)

### Comparison of Metal1D, Metal3D, MIB and BioMetAll

Many metal ion predictors exists that can be subdivided in two categories: binding site predictors and binding location predictors. The former label only the residues binding the ion, the latter also predict a location of the ion. 

In addition to Metal1D and Metal3D we also compared two recent predictors BioMetAll and MIB. MIB uses a fragment method to identify homologus binding sites to the motifs it finds in a given structure and will extract the location of the metal from the homologous structures in its database. The main performance regulator of MIB is the tscore cutoff which is a parameter for the template similiarity with higher values requiring higher similiarity. 
BioMetAll was calibrated on the PDB and places probes on a regular grid at all sites where they find the criteria to be fulfilled. For each collection of probes also a center of the probes is given which we used to assess performance as there is no individual ranking of the probes given by the program. The main parameter for BioMetAll is the cluster cutoff which indicates how many probes in reference to the largest cluster a specific cluster has. We used a cutoff of 0.5 requiring all chosen clusters to have at least 50% of the probes of the most popolous. 

For both tools the recommended settings match the accuracy of Metal3D p=0.75 with a lot more false positives. 
Metal1D offers high detection capabilites but also with a high number of false positives. 
While MIB also offers high precision, BioMetAll (using the cluster center) is not very precise with a MAD for correctly identified sites of 2.8 +- XX. Metal1D which identifies more sites than BioMetAll is slightly more precise than BioMetAll. MIB detects less sites but does so with high precision because it can use homologues sites to correctly place the metal ligand. BioMetAll also often provides probes that correctly identify the metal but as there is no ranking of the probes any probe could be closest to the actual location. 

![Comparison of Metal1D, Metal3D, BioMetAll and MIB on the testset used to train Metal1D and Metal3D. Predicted sites are counted if within 5 A of true metal location. False positive probes are clustered and counted once per cluster. ](images/metal3d_biometall_comparison.jpg)

![Mean absolute deviation of predicted location of zincs using Metal1D, Metal3D, BioMetAll and MIB on the testset used to train Metal1D and Metal3D for all correctly identified sites. ](images/mad_violin.jpg)



## Discussion

3D CNN model accuracy. Recent work EquiDock uses no sidechains at all, gets ligand RMSD where only 25% are under a 2 threshold, https://arxiv.org/pdf/2202.05146.pdf Mean RMSD 8.3 A, Centroid 42.4

## Conclusion

## Supplement

### Metal site detection using Metal1D


### Metal site detection using Metal3D

Figure MAD testset


### Comparison 

Figure MAD only good zincs comparison
Figure  only good zincs TP/FN/FP 