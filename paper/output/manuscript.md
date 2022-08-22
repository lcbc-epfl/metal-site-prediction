---
title: Accurate prediction of transition metal ion location via deep learning
keywords:
- metal
- protein engineering
- deep learning
lang: en-US
date-meta: '2022-08-22'
author-meta:
- Simon L. Dürr
- Andrea Levy
- Ursula Rothlisberger
header-includes: |-
  <!--
  Manubot generated metadata rendered from header-includes-template.html.
  Suggest improvements at https://github.com/manubot/manubot/blob/main/manubot/process/header-includes-template.html
  -->
  <meta name="dc.format" content="text/html" />
  <meta name="dc.title" content="Accurate prediction of transition metal ion location via deep learning" />
  <meta name="citation_title" content="Accurate prediction of transition metal ion location via deep learning" />
  <meta property="og:title" content="Accurate prediction of transition metal ion location via deep learning" />
  <meta property="twitter:title" content="Accurate prediction of transition metal ion location via deep learning" />
  <meta name="dc.date" content="2022-08-22" />
  <meta name="citation_publication_date" content="2022-08-22" />
  <meta name="dc.language" content="en-US" />
  <meta name="citation_language" content="en-US" />
  <meta name="dc.relation.ispartof" content="Manubot" />
  <meta name="dc.publisher" content="Manubot" />
  <meta name="citation_journal_title" content="Manubot" />
  <meta name="citation_technical_report_institution" content="Manubot" />
  <meta name="citation_author" content="Simon L. Dürr" />
  <meta name="citation_author_institution" content="Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland" />
  <meta name="citation_author_orcid" content="0000-0002-4304-8106" />
  <meta name="twitter:creator" content="@simonduerr" />
  <meta name="citation_author" content="Andrea Levy" />
  <meta name="citation_author_institution" content="Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland" />
  <meta name="citation_author_orcid" content="0000-0003-1255-859X" />
  <meta name="citation_author" content="Ursula Rothlisberger" />
  <meta name="citation_author_institution" content="Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland" />
  <meta name="citation_author_orcid" content="0000-0002-1704-8591" />
  <meta property="og:type" content="article" />
  <meta property="twitter:card" content="summary_large_image" />
  <link rel="icon" type="image/png" sizes="192x192" href="https://manubot.org/favicon-192x192.png" />
  <link rel="mask-icon" href="https://manubot.org/safari-pinned-tab.svg" color="#ad1457" />
  <meta name="theme-color" content="#ad1457" />
  <!-- end Manubot generated metadata -->
bibliography: []
manubot-output-bibliography: output/references.json
manubot-output-citekeys: output/citations.tsv
manubot-requests-cache-path: ci/cache/requests-cache
manubot-clear-requests-cache: false
...





<!--  -->
<div style="text-align:center">
<small><em>
This manuscript
was automatically generated
on August 22, 2022. <br>View [non-interactive version on BioRxiv](doi).
</em></small> 
</div>


<!--  -->
<!-- 
+ **Simon L. Dürr**  <br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [0000-0002-4304-8106](https://orcid.org/0000-0002-4304-8106)
    · ![Twitter icon](images/twitter.svg){.inline_icon}
    [simonduerr](https://twitter.com/simonduerr)<br>
  <small>
     Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland
  </small>

+ **Andrea Levy**  <br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [0000-0003-1255-859X](https://orcid.org/0000-0003-1255-859X)<br>
  <small>
     Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland
  </small>

+ **Ursula Rothlisberger** ✉ <br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [0000-0002-1704-8591](https://orcid.org/0000-0002-1704-8591)<br>
  <small>
     Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland
  </small>
 -->



[![ORCID icon](images/orcid.svg){height="11px" width="11px"}](https://orcid.org/0000-0002-4304-8106)
Simon L. Dürr <sup>1</sup>,
[![ORCID icon](images/orcid.svg){height="11px" width="11px"}](https://orcid.org/0000-0003-1255-859X)
Andrea Levy <sup>1</sup>,
[![ORCID icon](images/orcid.svg){height="11px" width="11px"}](https://orcid.org/0000-0002-1704-8591)
Ursula Rothlisberger <sup>1,✉</sup>


<!--<sup>☯</sup> --- These authors contributed equally. <br> -->
<sup>✉</sup> --- To whom correspondence should be addressed: ursula.roethlisberger@epfl.ch
<small>

###### Affiliations

1. Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland




## Abstract

Metal ions are essential cofactors for many proteins. In fact, currently, about half of the structurally characterized proteins contain a metal ion. Metal ions play a crucial role for many applications such as enzyme design or design of protein-protein interactions because they are biologically abundant, tether to the protein using strong interactions, and have favorable catalytic properties e.g. as Lewis acid. Computational design of metalloproteins is however hampered by the complex electronic structure of many biologically relevant metals such as zinc that can often not be accurately described using a classical force field. In this work, we develop two tools - Metal3D (based on 3D convolutional neural networks) and Metal1D (solely based on geometric criteria) to improve the identification and localization of zinc and other metal ions in experimental and computationally predicted protein structures. Comparison with other currently available tools shows that Metal3D is the most accurate metal ion location predictor to date outperforming geometric predictors including Metal1D by a wide margin using a single structure as input. Metal3D outputs a confidence metric for each predicted site and works on proteins with few homologes in the protein data bank. The predicted metal ion locations for Metal3D are within 0.70 ± 0.64 Å of the experimental locations with half of the sites below 0.5 Å. Metal3D predicts a global metal density that can be used for annotation of structures predicted using e.g. AlphaFold2 and a per residue metal density that can be used in protein design workflows for the location of suitable metal binding sites and rotamer sampling to create novel metalloproteins. Metal3D is available as easy to use webapp, notebook or commandline interface. 


## Introduction

Metalloproteins are ubiquitous in nature and are present in all major enzyme families [@doi:10.1021/cr400458x;@doi:10.1093/protein/gzw026].The metals predominantly found in biological systems are the first and second row alkali and earth alkali metals and the first row transition metals such as zinc and copper.  Zinc is the most common transition metal (present in ~10% of deposited structures) and can fulfill both a structural (e.g. in zinc finger proteins) or a catalytic role in up to trinuclear active sites. Zn<sup>2+</sup> is an excellent Lewis acid and is most often found in tetrahedral, pentavalent, or octahedral coordination. About 10 % of all reactions catalyzed by enzymes use zinc as cofactor[@pmid:18604568].

Metalloproteins are well studied because metal cofactors are essential for the function of many proteins and loss of this function is an important cause of diseases [@doi:10.1038/s42256-019-0119-z]. Industrial applications for metalloproteins capitalize on the favorable catalytic properties of the metal ion where the protein environment dictates (stereo)-selectivity [@doi:10.1126/science.aau3744, @doi:10.1021/bi201881p, @doi:10.1039/C5SC01065A; @doi:10.1038/nature17968; @doi:10.1038/s41570-021-00339-5]. To crystallize proteins, metal salts are also often added to the crystallization buffer as they can help in the formation of protein crystals overcoming the enthalpic cost of association of protein surfaces. Metal ion binding sites can be used to engineer protein-protein interactions (PPI) [@doi:10.1038/nchem.1290; @doi:10.1021/ja208015j; @doi:10.1021/ar900273t] and the hypothesis has been put forward that one origin of macromolecular complexity is the superficial binding of metal ions in early single domain proteins [@doi:10.1021/ar900273t].

While simple metal ion binding sites can be rapidly engineered because initial coordination on a protein surface can for example be achieved by creating an i, i+4 di-histidine site on an alpha-helix [@doi:10.1002/anie.202009226] or by placing cysteines in spatial proximity [@doi:10.1038/nchem.1201], the engineering of complex metal ion binding sites e.g. in the protein interior is considerably more difficult [@doi:10.1093/protein/gzw026;@doi:10.1021/ja208015j] as such sites are often supported by a network of hydrogen bonds. A complication for computational design of metalloproteins is the unavailability of good (non-bonded) force fields for zinc and other transition metals that accurately reproduce (e.g. tetrahedral) coordination with the correct coordination distances which renders design using e.g. Rosetta very difficult [@doi:10.1093/protein/gzw026; @doi:10.1021/jacs.0c01329]. In fact, the latest parametrization of the Rosetta energy function (ref2015) [@doi:10.1021/acs.jctc.7b00125] did not refit the parameters for the metal ions which originally are from CHARMM27 with empirically derived Lazaridis-Karplus solvation terms. To adequately treat metal sites in proteins quantum mechanical treatments such as in hybrid quantum mechanics/molecular mechanics (QM/MM) simulations [@doi:10.1021/cr500628b;@doi:10.1021/acs.jcim.1c01109] is needed whose computational cost is prohibitive for regular protein design tasks. QM/MM simulations can however be used to verify coordination chemistry for select candidate proteins [@doi:10.1021/jacs.7b10660]. On the other hand, neural network potentials have been developed for zinc however those require the experimental zinc location as input [@doi:10.3389/fchem.2021.692200].

Many tools exist to predict whether a protein contains metals (e.g. ZincFinder [@doi:10.1186/1471-2105-8-39]), which residues in the protein bind a metal (e.g. IonCom [@doi:10.1093/bioinformatics/btw396], MIB [@doi:10.1021/acs.jcim.6b00407]) and where the metal is bound (AlphaFill [@doi:10.1101/2021.11.26.470110], FindsiteMetal [@doi:10.1002/prot.22913], BioMetAll [@doi:10.1021/acs.jcim.0c00827],MIB [@doi:10.1021/acs.jcim.6b00407] ). The input for these predictors is based on sequence and/or structure information. Sequence-based predictors use pattern recognition to identify the amino acids which might bind a metal [@doi:10.1002/minf.201800169].
Structure-based methods use homology to known structures (MIB, Findsite-metal, AlphaFill) or distance features (BioMetAll) to infer the location of metals. Some tools like Findsite-metal or ZincFinder employ machine learning based approaches such as support vector machines.

Structure based deep learning approaches have been used in the field of protein research for a variety of applications such as protein structure prediction [@doi:10.1038/s41586-021-03819-2, @doi:10.1126/science.abj8754], prediction of identity of masked residues [@doi:10.1186/s12859-017-1702-0; @doi:10.1021/acssynbio.0c00345; @doi:10.1038/s41467-022-28313-9 ], functional site prediction [@doi:10.1093/bioinformatics/bty813; @doi:10.1038/s41467-021-24070-3], for ranking of docking poses [@doi:10.1038/s41467-021-27396-0; @doi:10.1038/s41592-019-0666-6], prediction of the location of ligands [@doi:10.1093/bioinformatics/btx350; @doi:10.1093/bioinformatics/bty583; @doi:10.48550/arXiv.2202.05146 ; @doi:10.1038/s41592-019-0666-6; @doi:10.1021/acs.jcim.2c00306], and prediction of effects of mutations for stability and disease [@doi:10.1038/s42256-019-0119-z, @doi:10.1371/journal.pcbi.1008291].
Current state of the art predictors for metal location are MIB [@doi:10.1371/journal.pone.0039252;@doi:10.1021/acs.jcim.6b00407], which combines structural and sequence information in the “Fragment Transformation Method” to search for homologous sites in its database, and BioMetAll [@doi:10.1021/acs.jcim.0c00827], a geometrical predictor based on backbone preorganization. Both methods have significant drawbacks: MIB excludes metal sites with less than 2 coordination partners from its analysis and is limited by the availability of templates in its database. BioMetAll does not use templates but provides many possible locations for putative binding sites on a regular grid. The individual probes in BioMetAll do not have a confidence metric therefore only allowing to rank sites by the number of probes found, which results in a large uncertainty in the position. Both tools suffer from many false positives.  In this work, we present two metal ion location predictors that do not suffer from these drawbacks. The deep learning based Metal3D predictor operates on a voxelized representation of the protein environment and predicts a per residue metal density that can be averaged to get a smooth metal probability density over the whole protein. The distance based predictor Metal1D predicts the location of metals using distances mined from the protein data bank (PDB) directly predicting coordinates of the putative metal binding site. These tools pave the way to perform in silico design of metal ion binding sites without relying on predefined geometrical rules or expensive quantum mechanical calculations. 


## Results

A dataset of experimental high resolution crystal structures (2085 structures/252324 voxelized environments) containing zinc sites was used for training of the geometric predictor Metal1D and the deep learning predictor Metal3D  (Figure @fig:method). For training, we used the crystal environment including crystal contacts. For predictions, the biological assembly was used. 

![**Workflow of Metal3D and Metal1D** **A** Training of Metal3D and Metal1D is based on experimental Zn<sup>2+</sup> sites. Metal1D extracts coordination environments from LINK records, Metal3D is a fully convolutional 3DCNN trained to predict the metal density from voxelized protein environments. **B** In inference mode Metal3D predicts the location of a metal ion by computing per residue metal densities and then averaging them to obtain a global metal density for the input proteins. The ions can then be placed using the weighted average of voxels above a cutoff. For Metal1D all residues in the protein are scanned for compatibility with the probability map. Metals are placed at the geometric center of residues with high scores according to the probability map. A final ranking of sites is obtained using the probability map. 
](images/Metal3D_Metal1D_method.png){#fig:method}

### Metal3D

Metal3D takes a protein structure and a set of residues as input, voxelizes the environment around each of the residues and predicts the per residue metal density. The predicted per residue densities (within a 16 x 16 x 16 Å<sup>3</sup> volume)  can then be averaged to yield a zinc density for the whole protein. At high probability cutoffs the predicted metal densities are spherical (Figure @fig:2cbaprobabilities E), at low probability cutoffs the predicted densities are non-regular (Figure @fig:2cbaprobabilities A).

![**Metal3D probability density** Probability evolution in HCA2 (PDB 2CBA) for different probability cutoffs A) p=0.1 B) p=0.2 C) p=0.3 D) p=0.4 E) p=0.75.](images/2CBA_probabilities.png){#fig:2cbaprobabilities}


We evaluated the quality of the metal densities generated by the model with the discretized Jaccard similarity (Figure @fig:jaccard) for all environments in the test set. We noticed that at the edges of the residue-centered output densities often spurious density is predicted wherefore we evaluated the similarity of the test set metal density and the predicted metal probability density taking into account a smaller box with zeroed outer edges. Figure  @fig:jaccard shows that the similarity of the boxes does not depend much on the probability cutoff chosen with higher cutoffs yielding slightly higher discretized Jaccard similarity values (0.02 - 0.04 difference between p=0.5 and p=0.9). Reducing the size of the analyzed boxes (i.e trimming of the edges) increases the Jaccard similarity from ≈ 0.64 to 0.88 showing that the metal density in the center of the box is more accurate than the density at the edges.

Metal3D is available as self-contained notebook on [GoogleColab](https://colab.research.google.com/github/lcbc-epfl/metal-site-prediction/blob/main/Metal3D/ColabMetal.ipynb) and on [Huggingface Spaces](https://hf.space/simonduerr/metal3d).

<script type="module" src="https://gradio.s3-us-west-2.amazonaws.com/3.0.18/gradio.js"></script>
<gradio-app space="simonduerr/metal3d"></gradio-app>


### Metal1D
The statistical analysis for the geometric predictor uses the `LINK` records present in deposited PDB structures. A probability map for all zinc coordination motifs was extracted from all training structures (Figure @fig:method A). The mean coordination distance in the training set was found to be  2.2 ± 0.2 Å, and the default search radius for the predictions was therefore set to 5.5 Å (Table @tbl:metal1dreferencepoints). In total 208 different environments with more than 5 different proteins (at 30 % sequence identity) were identified. 
Metal1D is available as self-contained notebook on [Google Colab](https://colab.research.google.com/github/lcbc-epfl/metal-site-prediction/blob/main/Metal1D/ColabMetal1D.ipynb).



### Comparison of Metal1D, Metal3D, MIB and BioMetAll

Existing metal ion predictors can be subdivided into two categories: binding site predictors and binding location predictors. The former identify only the residues binding the ion, the latter predict the coordinates of the metal ion itself. Both Metal1D and Metal3D can predict the coordinates of putative binding sites. We therefore assessed their performance by comparing to recent binding location predictors with available code/webserver: BioMetAll [@doi:10.1021/acs.jcim.0c00827]  and MIB [@doi:10.1021/acs.jcim.6b00407].  The main tuning parameter of MIB is the template similarity `t`,  with higher values requiring higher similarity of the templates available for the search in structurally homologous metalloproteins. BioMetAll on the other hand was calibrated on available protein structures and places probes on a regular grid at all sites where the criteria for metal binding are fulfilled. The main adjustable parameter for BioMetAll is the cluster cutoff `c`, which indicates how many probes in reference to the largest cluster a specific cluster has. We used the recommended cutoff of 0.5 requiring all chosen clusters to have at least 50% of the probes of the most populous cluster and used the cluster center to compute distances.

We first investigated the potential of all tools to detect the location of a zinc ion binding site in a binary fashion (zinc site or no zinc site). We defined a correctly identified binding site (true positive, TP) as a prediction within 5 Å of an experimental zinc site. In case a tool predicted no metal within the 5 Å radius, we counted this site as false negative (FN). False positive (FP) predictions, i.e sites where a metal was placed spuriously, were clustered in a 5 Å radius and counted once per cluster. All tools were assessed against the held out test biological assemblies for Metal3D and Metal1D. When the performance of MIB (`t=1.25`) and BioMetAll is compared against Metal3D with probability cutoff p=0.75 we find that Metal3D identifies more sites (85) than MIB (78) or BioMetAll (75) with a much lower number of false positives (Figure @fig:comparison).
MIB predicts 180 false positive sites, BioMetAll 134 sites whereas Metal3D only predicts 9 false positive sites at the p=0.75 cutoff. Metal1D (`t=0.5`) offers similar detection capabilities (78 sites detected) with a lower number of false positives (47) compared to MIB and BioMetAll. 
We removed 56 sites from the list of zinc sites in the test set (189 total) that had less than 2 unique protein ligands within 2.8 Å of the experimental zinc location. The amount of correct predictions in this reduced set is almost unchanged for all tools (Figure @fig:comparison) indicating that most tools correctly predict sites if they have 2 or more protein ligands. For Metal3D at p=0.75 and p=0.9 as well as MIB with t=1.9 all sites that are correctly predicted to contain a metal are sites with more than 2 protein ligands. The number of false negatives is reduced for all tools by about 50 sites indicating that most tools do not predict these crystallographic artifacts that might depend on additional coordinating residues from an adjacent molecule in the crystal. Of all tools, Metal3D has the least false positives (1 FP at p=0.9) and the highest number of detected sites (110 at p=0.25). The single false positive at p=0.9 does not contain a zinc ion but is a calcium binding site with three aspartates and one backbone carbonyl ligand (Figure @fig:4jjjFP).

![**Identification of metal sites** Comparison of Metal1D, Metal3D, BioMetAll and MIB on the test set held out from training of Metal1D and Metal3D. Predicted sites are counted as true positives (TP) if they are within 5 Å of a true metal location and as false negatives (FN) otherwise. False positive (FP) probes are clustered and counted once per cluster. *For MIB we used 2 structures less because the server did not accept these structures.](images/all_2coord_metal_1d_3d_min_biometall_comparison.jpg){#fig:comparison}


After assessment of how many sites the tools predict, another crucial metric is the spatial precision of the predictions. For the correctly identified sites (TP) we measured the mean absolute distance (MAD) between experimental and predicted position (Figure @fig:selectivity-mad). The MAD for Metal3D at p=0.9 is 0.70 ± 0.64 Å and 0.73 ± 0.65 Å at p=0.25 indicating that low confidence predictions are still accurately placed inside the protein. The median MAD of predictions for Metal3D at p=0.9 is 0.52 Å indicating that for half of the predictions the model predicts at or better than the grid resolution of 0.5 Å.

BioMetAll is not very precise with a MAD for correctly identified sites of 2.80 ± 1.30 Å. BioMetAll predicts many possible locations per cluster with some of them much closer to the experimental metal binding site than the cluster center. However, it does not provide any ranking of the probes within a cluster and therefore the cluster center was used for the distance calculation. Metal1D (MAD 2.06 ± 1.33 Å) which identifies more sites than BioMetAll is also more precise than BioMetAll. MIB t=1.9 detects sites with high precision (MAD 0.77 ± 1.09 Å) but it relies on the existence of homologous sites to align the found sites.

![**Precision of predicted sites and selectivity for Zn<sup>2+</sup>** **A** Mean absolute deviation (MAD) of predicted zinc ion locations using Metal1D, Metal3D, BioMetAll and MIB on the test set used to train Metal1D and Metal3D for all correctly identified (TP) sites. n is the number of sites predicted by the tool. *For MIB we used 2 structures less because the server did not accept these structures. For each tool the whisker plot indicates the median (white dot) and the first quartiles (black box). **B** Recall for the zinc test set and 25 randomly drawn structures for other transition, alkali and earth-alkali metal ions for Metal3D using p=0.5 as cutoff. ](images/Figure3_selectivity_mad.png){#fig:selectivity-mad}

### Selectivity for other metals

Both Metal3D and Metal1D were exclusively trained on zinc and we assessed their performance on sodium (Na+, PDB code NA), potassium (K+, PDB code K), calcium (Ca2+, PDB code CA), magnesium (Mg2+, PDB code MG), and various transition metals (Fe2+, Fe3+, Co2+, Cu2+, Cu+, Mn2+ with corresponding PDB codes FE2, FE, CO, CU, CU1, MN, respectively) from 25 randomly drawn structures from the clustered PDB at 30% identity. Only sites with at least 3 unique protein ligands were used for the analysis to exclude crystallographic artifacts and use only highly defined sites which should exhibit most selectivity towards a specific metal. Figure @fig:selectivity-mad B shows that recall for Metal3D is high for all transition metals, meaning that the model correctly finds most sites even though it was only trained on zinc. For the alkali and earth alkali metals recall is much lower as the model only finds some sites. The mean probability for found zinc structures (ZN p=0.95 ± 0.10) in the test set is higher than for the other transition metals (Figure @fig:selectivity-probability-metal3d) and significantly higher than the probability for alkali metals (NA p=0.61 ± 0.10, K p=0.79 ± 0.16) while the probability for the earth alkali metals is slightly higher with MG (p=0.77 ± 0.16)  similar to CA (p=0.73 ± 0.16). The MAD for each found metal site is again lowest for zinc (0.56 ± 0.59 Å). The MAD for the found sodium (n=2) and potassium (n=5) sites are as low as for the other transition metals. The only metal with significantly higher MAD (1.45 ± 0.93 Å) is CU1 (Figure @fig:selectivity-distance-metal3d).

The only two structures where a sodium is detected by Metal3D  (2OKQ[@doi:10.2210/pdb2OKQ/pdb], 6KFN[@doi:10.2210/pdb6KFN/pdb]) have at least 2 side chain coordinating ligand atoms and only one backbone (2OKQ) or no backbone ligand atom (6KFN). Canonical sodium binding sites e.g. such as in PDB 4I0W [@doi:10.2210/pdb4I0W/pdb] with two coordinating backbone carbonyl oxygen atoms and one asparagine side chain have probabilities around 5 % and are basically indistinguishable from background noise of the model.
For Metal1D overall recall is lower with similar differences in the detection of main group metals versus transition metals (Figures @fig:selectivity-metal1d and  @fig:selectivity-distance-metal1d).


### Applications

After having evaluated the accuracy of Metal3D on held out test structures we also investigated possible uses in downstream applications such as protein function annotation and protein design. 

#### Alpha Fold 

AlphaFold2 often predicts side chains in metal ion binding sites in the holo conformation [@doi:10.1038/s41586-021-03819-2]. Tools like AlphaFill [@doi:10.1101/2021.11.26.470110]use structural homology to transplant metals from similar PDB structures to the predicted structure. Metal3D does not require explicit homology based on sequence or structural alignment like AlphaFill so it is potentially suited to annotate the dark proteome that is now accessible from the AlphaFold database with metal binding sites. Metal3D identifies both the catalytic site (1) and the zinc finger (2) for the example (PDB 3RZV[@doi:10.2210/pdb3RZV/pdb] , Figure @fig:alphafold A) used in ref [@doi:10.1101/2021.11.26.470110] with high probability (p=0.99) even though one of the sites in the AlphaFold model is slightly disordered with one of the binding residues in the solvent facing conformation (D309). The distances between predicted and modeled metal locations for Metal3D are 0.22 Å and 0.37 Å, for AlphaFill they are 0.21 Å and 0.41 Å. 

AlphaFill uses a 25% sequence identity cutoff which can be problematic for certain proteins with no structurally characterized homologues. For human palmitoyltransferase ZDHHC23 (Uniprot Q8IYP9) a high confidence AlphaFold2 prediction exists but AlphaFill cannot place the zinc ions because the sequence identity is 24% to the closest PDB structure (PDB 6BMS [@doi:10.1126/science.aao6326]), i.e below the 25% cutoff. For the identical site in another human palmitoyltransferase ZDHHC15 (Uniprot Q96MV8) AlphaFill is able to place the metal because of higher sequence identity to 6BMS (64%) (Figure @fig:alphafold B). For ZDHHC23 Metal3D is able to place the metal with high confidence (MAD 0.75 Å for site 1 and 0.48 Å for site 2, p>0.99 ) based on the single input structure alone. 





![**Annotation of AlphaFold2 structures** **A** Predicted metal binding sites (a and b) from Metal3D, respectively AlphaFill compared to the experimentally found zinc positions for Uniprot O95630. Metal3D places the metal with high accuracy even if sidechains are not perfectly predicted by AlphaFold for site 1 **B** palmitoyltransferase ZDHHC23 (Uniprot Q8IYP9)  and ZDHHC15 (Uniprot Q96MV8). AlphaFill can only place the metal for ZDHHC15 because sequence identity for ZDHHC23 is only 24 %. Probability isosurfaces from Metal3D for both structures at p=0.6, colored in gray.](images/3rzv_6bsm_AF2_alphafill_metal3d.png){#fig:alphafold}



#### Metalloprotein engineering

Human carbonic anhydrase II (HCA2) is a well studied metalloenzyme with a rich amount of mutational data available. For  the crystal structure of the wildtype enzyme (PDB 2CBA [@doi:10.2210/pdb2cba/pdb]),  Metal3D recapitulates the location of the active site metal with a RMSD to the true metal location of 0.21 Å with a probability of p=0.99. At lower probability cutoffs (p<0.4) the probability map indicates further putative metal ion binding sites with interactions mediated by surface residues (e.g. H36, D110, p=0.22) (Figure @fig:2cbaprobabilities). 


<iframe src="interactive/index.html" title="2CBA probabilities" width="100%" height="900px" frameborder="0"></iframe> 


To investigate the capabilities for protein engineering we used mutational data for first and second shell mutants of the active site residues in HCA2 with corresponding K<sub>d</sub> values from a colorimetric assay [@pmid:3887984]. For most mutants no crystal structures are available so we used the structure builder in the EVOLVE package <!--REF-->  to choose the most favorable rotamer for each single point mutation based on the EVOLVE-ddg energy function <!--REF-->  with explicit zinc present (modeled using a dummy atom  approach [@pmid:11106157]). The analysis was run for every single mutant and the resulting probability maps from Metal3D were analyzed. For the analysis we used the maximum predicted probability as a surrogate to estimate relative changes in K<sub>d</sub>. For mutants that decrease zinc binding drastically we observe a drop in the maximum probability predicted by Metal3D (Figure @fig:hca-kd).The lowest probability mutants are H119N and H119Q with p=0.23 and 0.38. The mutant with the largest loss in zinc affinity H94A has a zinc binding probability of p=0.6. Conservative changes to the primary coordination motif (e.g. H &rarr; C) reduce the predicted probability by 10 - 30 %. For second shell mutants the influence of the mutations is less drastic with only minor changes in the predicted probabilities. 

![**Protein design application** Experimentally measured K<sub>d</sub> values [@doi:10.1021/bi00255a003;@doi:10.1021/ja00079a046;@doi:10.1021/bi00089a005;@doi:10.1073/pnas.92.11.5017;@doi:10.1021/bi9526692] for 1<sup>st</sup> and 2<sup>nd</sup> shell active site mutants of HCA2 and predicted max probability for zinc using Metal3D.](images/kd_vs_p_nolog_newmethod_newmodel_0.5.png){#fig:hca-kd}



## Discussion

Metal3D predicts the probability distribution of zinc ions in protein crystal structures based on a neuronal network model trained on natural protein environments. The model performs a segmentation task to determine if a specific point in the input space contains a zinc ion or not. Metal3D predicts zinc ion sites with high accuracy making use of high resolution crystal structures (<2.5 Å). The use of high resolution structures is necessary because at resolutions greater than the average zinc ligand coordination distance (2.2 Å) the uncertainty of the zinc location noticeably increases [@doi:10.1021/ic401072d] which would likely hamper the accuracy of the site prediction.


In contrast to currently available tools, for Metal3D, it is not necessary to filter the training examples for certain coordination requirements (i.e only sites with at least 2 protein ligands). The model thus sees the whole diversity of zinc ion sites present in the PDB.  Such a model is advantageous since metalloprotein design workflows require models to score the full continuum of zinc sites starting from a suboptimal binding site only populated at high metal concentration to a highly organized zinc site in an enzyme with nanomolar metal affinity. The predicted probability can be used as a confidence metric or as an optimization target where mutations are made to increase probability of zinc binding.

### Site quality
The fraction of artifactual zinc binding sites in the PDB is estimated to be about 1/3 [@doi:10.1021/ic401072d] similar to our test set used with 70% (133) well coordinated zinc sites with at least 2 distinct protein ligands. To reduce the amount of artifactual sites in the training set we presented the model with as many complete sites as possible by using crystal symmetry to add adjacent coordinating protein chains. The frequency of artifacts in the training set is therefore much lower than 30 %. The sites which still remain incomplete or that are wrongly modeled and not excluded through the resolution cutoffs and filtering procedures likely present only a small fraction of the training set and their signal is drowned out by the numerical superiority of the correctly modeled sites. If the model is used on artifactual sites or partially disordered ones it can still predict the metal location with high spatial accuracy but often indicates a lower confidence for the prediction (Figures @fig:2cbaprobabilities and @fig:alphafold). 

Metal ion locators that rely on homology such as MIB perform worse on partial binding sites because reducing the quality of the available templates by including 1- or 2-coordinate sites would yield many false positives (similar to including less homologous structures for the template search). The deep learning based Metal3D can likely circumvent this because it does not require any engineered features to predict the location of the metal and learns directly from a full representation of the environment surrounding the binding site. This allows looking at low confidence sites in the context of a given environment. 

### Influence of non-protein ligands

Exogenous ligands play an important role for metals in biology as all empty coordination sites of metals are filled with water molecules in case there is no other exogenous ligand with higher affinity present (e.g. a thiol).  Like other predictors, both Metal1D and Metal3D do not consider water molecules or other ligands in the input as the quality of ligand molecules in the PDB varies [@pmid:3736419;@doi:10.1021/acs.jcim.2c00306].  In addition, other potential sources of input such as AlphaFold do not provide explicit waters wherefore models should not rely on water as an input source. It is also not possible to use in silico water predictions because common water placement algorithms to place deep waters [@doi:10.1002/prot.25081; @doi:10.1371/journal.pone.0172743; @doi:10.1021/acs.jcim.2c00306] either rely on metal ions being present in the input or ignore them completely.  Moreover, in protein design algorithms, water is usually only implicitly modeled (e.g. in Rosetta). 

For Metal3D, the input channel that encodes the total heavy atom density also encodes an implicit water density where all empty space can be interpreted as the solvent. For Metal1D, the contribution of water molecules is considered in an implicit way when the score is assigned to a site by considering coordinations including water compatible with the one observed (e.g. a HIS<sub>3</sub>Wat site is equivalent to a His<sub>3</sub> site for the scoring).


### Choice of architecture
This work is the first to report a modern deep learning based model destined for identification of metal ligands in proteins. Similar approaches have been used in the more general field of protein-ligand docking where a variety of architectures and representations have been used. 3D CNN based approaches such as LigVoxel [@doi:10.1093/bioinformatics/bty583] and DeepSite [@doi:10.1093/bioinformatics/btx350)] commonly use a resolution of 1 Å and similar input features as our model to predict the ligand density. However, predicting the density of a multi-atomic ligand is more complex than predicting the density of mononuclear metal ions. We therefore did not deem it necessary to include a conditioning on how many metal ions are present in the box and rather chose to reflect this in the training data where the model needs to learn that only about half of the environments it sees contain one or more metals. This choice is validated by the fact that the output probability densities at sufficiently high probability cutoffs are spherical with their radius approximately matching the van der Waals radius of zinc.

Mesh convolutional neural networks trained on a protein surface representation [@doi:10.1038/s41592-019-0666-6] also have been used to predict the location and identity of protein ligands but this approach can only label the regions of the surface that bind the metal ion and is conceptually not able to return the exact location of the metal. Some metal ion binding sites are also heavily buried inside proteins as they mediate structural stability rendering them inaccessible to a surface based approach. The most recent approaches such as EquiBind [@doi:10.48550/arXiv.2202.05146] use equivariant neural networks such as En-Transformer [@doi:10.48550/arXiv.2102.09844] to predict binding keypoints (defined as 1/2 distance between the Cɑ of the binding residue and a ligand atom). Explicit side chains are still too expensive for such models and these models assume a fixed known stoichiometry of the protein and ligand. Metal3D can also deal with proteins that do not bind a metal and does not assume that the amount of ions is known. The lack of explicit side chain information renders equivariant models unsuitable for the design of complex metal ion binding sites supported by an intricate network of hydrogen bonds that need to be positioned with sub-angstrom accuracy. Our model in contrast is less data- and compute-efficient than approaches representing the protein as graph due to the need to voxelize the input and provide different rotations of the input environment in training but the overall processing time for our model is still low taking typically 25 seconds for a 250 residue protein on a multicore GPU workstation (20 CPUs, GTX2070).
Sequence based models [@doi:10.1101/2021.12.22.473759;@doi:10.1038/s41598-018-34533-1] can only use coevolution signals to infer residues in spatial proximity that can bind a metal. This might be difficult when it comes to ranking similar amino acids such as aspartate and glutamate or even ranking different rotamers where sub-angstrom level precision is needed to identify the mutant with the highest affinity for zinc. 

### Selectivity

In terms of selectivity both of our methods have a clear preference for transition metals over main group metals after having been trained exclusively on zinc binding sites. The only sites that Metal3D identifies for sodium in the test set are the ones that have side chain ligands. Many sodium and potassium sites are using backbone carbonyl coordination exclusively, which is not common for zinc and those sites are therefore not detected. Both of our methods could be rapidly adapted to predict not only location but also the identity of the metal. In the framework of Metal3D even a semantic metal prediction would be possible where the same model predicts different output channels for each metal it was trained on. To achieve perfect selectivity using such a model will be difficult because sometimes non-native metals are used for crystallization experiments. In this work we chose to work exclusively with zinc because it is the most redox stable transition metal and because many training examples are available. 


### Application for protein design
Protein design using 3DCNNs trained on residue identity has been successfully demonstrated and we anticipate that our model could be seamlessly integrated into such a workflow[@doi:10.1038/s41467-022-28313-9] to enable fully deep learning based design of metalloproteins. We are currently also investigating the combination of Metal3D combined with a classic energy-based genetic algorithm-based optimization to make design of metalloproteins [@doi:10.1021/jacs.7b10660] easier without having to explicitly model the metal to compute the stability of the protein. As the model computes a probability density per residue it can be readily integrated into established software like Rosetta relying on rotamer sampling.


The HCA2 application demonstrates the utility of Metal3D for protein engineering  (Figure @fig:2cbaprobabilities). The thermodynamics of metal ion binding to proteins are complicated [@doi:10.1021/ic301645j] and there are currently no high-throughput based experimental approaches that could generate a dataset large enough to train a model directly on predicting K<sub>d</sub>. The data we use were obtained from a colorimetric assay with very high affinity of zinc in the picomolar range [@doi:10.1021/bi00255a003;@doi:10.1021/ja00079a046;@doi:10.1021/bi00089a005;@doi:10.1073/pnas.92.11.5017;@doi:10.1021/bi9526692]. More recent studies using ITC [@doi:10.1021/ic301645j]  instead of the colorimetric assay indicate lower K<sub>d</sub> values  in the nanomolar range for wild type HCA2. We can therefore only use the colorimetric data to estimate how well the model can recapitulate relative changes in the K<sub>d</sub> for different mutations in the first and second shell of a prototypical metalloprotein. 


Metal3D allows moving away from using rational approaches such as the i, i + 4 di His motifs used for the assembly and stabilization of metalloproteins to a fully automated approach where potential metal binding configurations can be scored computationally. [@doi:10.1126/science.8346440; @doi:10.1126/science.1648261; @doi:10.1038/nsb723]

### Metal1D vs. Metal3D
Metal1D is inferior to Metal3D for the prediction of metal ion binding sites because it produces more false positives while at the same time detecting fewer metal sites. Also the positioning of sites is somewhat imprecise. This demonstrates the inherent limitation of using solely distance based features for prediction of metal location. BioMetAll which is the tool most similar to Metal1D also suffers from many false positive predictions. In contrast, Metal1D is more data-efficient than Metal3D and provides predictions faster. 





 <!-- KD values discussion (new ITC data vs. old data by Kiefer, Fierke etc): 
 
 Old method: Enzyme-
bound zinc (E—Zn) was quantitated using the colorimetric
4-(2-pyridylazo)resorcinol (PAR) method of Hunt et al.
(1984) and measuring the absorbance at 500 nm. 4-(2-Pyridylazo)resorcinol (PAR) is a dibasic acid that forms the protonated complexes with most metal ions. It serves as a metallochromic indicator and is suitable as a chromogenic agent for the quantitative determination of over 50 elements.

removing unbound zinc by chromatography on a
PD-10 column, and measuring the protein concentration and
bound zinc concentration in the eluant using the PAR assay
(Hunt et al., 1984).
The concentration of free zinc in the
dialysis buffer was calculated from the Tris—zinc stability
constants (Dawson et al., 1986). The dissociation constant
was calculated using KaleidaGraph program with eq:
[E-ZN]/[E]tot = C/(1+Kd/[Zn]free)

Newer method: [@doi:10.1021/ic301645j]
 -->


## Conclusion

We present two metal ion location predictors: Metal3D based on 3D convolutional neural networks and Metal1D based on distances and amino acid propensity maps. Metal3D is the first tool with sub-angstrom level precision to predict the location of metal ions in proteins that does not rely on searching for structurally homologous proteins in a database. We therefore anticipate different applications such as protein-function annotation of structures predicted using AlphaFold2 [@doi:10.1038/s41586-021-03828-1], integration in protein design software and detection of cryptic metal binding sites that can be used to engineer PPIs. Such cryptic metal ion binding sites in common drug targets could also be used to engineer novel metallodrugs. Many of these applications will allow us to explore the still vastly untapped potential of proteins as large multi-dentate metal ligands with programmable surfaces.



<!-- Zinc is generally available in the intracellular volume and extracellularly with the total concentration of zinc in cells is estimated to be about about 200 μM [@doi:10.1021/cr800556u]. The real available concentration of free zinc is picomolar[@doi:10.1021/cr800556u] wherefore it is crucial that the design tool produces well designed sites with high spatial precision of the metal ion.  -->




## Materials and Methods

### Dataset
The input PDB files for training were obtained from the RCSB [@doi:10.1093/nar/28.1.235] protein databank (downloaded 5th March 2021). We use a clustering of the structures at 30% sequence identity using mmseqs2 [@doi:10.1038/nbt.3988] to largely remove sequence and structural redundancy in the input dataset. 
For each cluster, we check whether a zinc is contained in one of the structures, whether the resolution of these structures is better than 2.5 Å, if the experimental method is x-ray crystallography and whether the structure does not contain nucleic acids. If there are multiple structures fulfilling these criteria, the highest resolution structure is used. All structures larger than 3000 residues are discarded. We always use the first biological assembly to sample the training environments. The structures were stripped of all exogenous ligands except for zinc . If there are multiple models with e.g. alternative residue conformations for a given structure, the first one is used. For each biological assembly we used the symmetry of the asymmetric unit to generate a protein structure that contains all neighboring copies of the protein in the crystal such that metal sites at crystal contacts are fully coordinated.

The train/val/test split was performed based on sequence identity using `easy-search` in mmseqs2. All proteins that had no (partial) sequence overlap with any other protein in the dataset were put into the test/val set (85 proteins) which we further split into a test set of 59 structures and a validation set of 26 structures. The training set contained 2085 structures. (Supplemental Data 1). 

For the analysis, we always used the biological assembly and not the symmetry augmented structure. For the specificity analysis with respect to other transition metals, clusters from the PDB were randomly sampled to extract 25 biological assemblies per metal. 

By default all zinc sites in the test and validation set were used for the analysis. Since some of the sites might be affected by the crystallization conditions, we also created a subset of all sites that contained at least 2 amino acid ligands to largely exclude crystallization artifacts.  To analyze metal ion selectivity, we selected sites with at least 3 unique protein ligands to only use biologically significant sites with a high degree of metal preorganization as such sites should exhibit more selectivity for specific metals compared to sites with only 2 unique protein-ligands. 

### Metal 1D

Metal1D uses a probability map derived from `LINK` records in protein structures (Figure @fig:method). The `LINK` section of a PDB file specifies the connectivity between zinc (or any other ligand) and the amino acids of the protein, and each `LINK` record specifies one linkage. This is an extension of the approach by Barber-Zucker *et al.*, [@doi:10.1038/s41598-017-16777-5] in which `LINK` records were used to investigate the propensity of transition metals to bind different amino acids. 


Using the training set we generated a probability map for the propensity of different coordination environments to bind a zinc (e.g CCCC, CCHH etc.). For each zinc ion the coordination is extracted from  the`LINK` records excluding records involving only single amino acids (weak binding sites). Also, `LINK` records containing water molecules are excluded because of the difficulties in placing water molecules a posteriori in 3D structures when metal ions are present and because data quality of modelled water molecules varies. The probability map contains the counts of coordination environments found. 
Making a prediction using Metal1D consists of two steps (Figure @fig:method): Identification of possible metal coordinating residues in the structure via the scoring of each amino acid, and scoring of the likelihood of coordination for putative sites predicted by placing a metal between the identified coordinating residues. 


The protein structure is analyzed using the BioPandas python library [@doi:10.21105/joss.00279]. To identify coordinating residues, a per residue score is assigned by performing a geometrical search from a reference point, defined as the coordinate of the most probable metal binding atom, within a search radius considered as roughly twice the typical distance between the metal ion and the binding atom of amino acids in proteins (2.2 ± 0.2 Å as determined from `LINK` records). The search radius used was 5.5 Å in order to be able to take into account also deviations from the ideal coordination. In the case of amino acids which present more than one putative coordinating atom, such as e.g. histidine, the mid-point between the donor atoms is used as reference point and the search radius is enlarged accordingly. The atoms used as reference points for each amino acid and the increase in the search radius are reported in Supplemental Table @tbl:metal1dreferencepoints  .
The score is assigned to each amino acid considering all the other reference points of other amino acids within the search radius, and summing the probabilities in the probability map for coordinations compatible with the one observed. In the ideal case, a score of 1 corresponds to an amino acid surrounded by all possible coordinating amino acids observed in the probability map. In practice, scores result between 0 and < 1. Once all amino acids in the chain are scored, the metal location predictions are made grouping the highest-scored amino acids in clusters (defined as the ones within the chosen threshold with respect to the highest-scored one) based on distance. This is done using  `scipy.spatial.distance_matrix` and grouping together highest-scored amino acids closer than twice the search radius. For each cluster, a site prediction is made as a weighted average between the coordinates of the reference point of each amino acid, using as weighting factor the amino acid score. For isolated amino acids with a high score (e.g. a single histidine) the same score is assigned to the closest reference point from another amino acid, to be able to compute the position of the metal as before. Possible artifacts resulting from this fictitious score are resolved in the final step of the prediction.

After the metal has been placed the likelihood of the putative sites can be assessed by performing a geometrical search centered on the predicted metal coordinates (within 60% of the search radius, i.e 3.3 Å) and a final score is now assigned to the site. The final score is assigned in the same way as the amino acid scores based on the probability map, and has the advantage of being able to sort the predicted metal sites based on their frequency in the training set. A cutoff parameter is used to exclude sites with a probability lower than a certain threshold with respect to the highest-scored one. This final scoring also mitigates the errors which can be introduced by calculating the coordinates of the site simply as a weighted average excluding or assigning a low probability to the site ending in unfavorable positions in space.

### Metal 3D

#### Voxelization
We used the moleculekit python library [@doi:10.1021/acs.jctc.6b00049; @doi:10.1093/bioinformatics/bty583] to voxelize the input structures into 3D grids. 8 different input channels are used: aromatic, hydrophobic, positive ionizable, negative ionizable, hbond donor, hbond acceptor, occupancy, and metal ion binding site chain (Supplemental Table @tbl:voxelchannels).  The channels are assigned using AutoDockVina atom names and a boolean mask. For each atom matching one of the categories a pair correlation function centered on the atom is used to assign the voxel value [@doi:10.1093/bioinformatics/bty583]. For the target tensor only the zinc ions were used for the voxelization. The target tensor was discretized setting any voxel above 0.05 to 1 (true location of zinc), all other to 0 (no zinc). We used a box size of 16 Å centered on the Cɑ atom of a residue, rotating each environment randomly for training before voxelization. The voxel grid used a 0.5 Å resolution for the input and target tensors. Any alternative side chain conformations modeled were discarded keeping only the highest occupancy. For the voxelization only heavy atoms were used. For all structures selected for the respective sets we partitioned the residues of the protein into residues within 12 Å of a zinc ion and those further away (based on the distance to the Cɑ atom). A single zinc site will therefore be present many times in the dataset but each time translated and rotated in the box. A balanced set of examples was used sampling equal numbers of residues that are close to a zinc and residues randomly drawn from the non-zinc binding residues. The sampling of residues is based on the biological assembly of the protein, the voxelization is based on the full 3D structure including neighboring asymmetric units in the crystal structure. The environments are precomputed and stored using lxf compression in HDF5 files for concurrent access during training. In total, 252324 environments were voxelized for the training set, 6550 for the test set, 3067 for the validation set. The voxelization was implemented using ray [@doi:10.48550/arXiv.1712.05889].

#### Model training
We used PyTorch 1.10 [@doi:10.48550/arXiv.1912.01703] to train the model.  All layers of the network are convolutional layers with filter size 1.5 Å except for the fifth layer where a 8 Å filter is used to capture long range interactions. We use zero padding to keep the size of the boxes constant. Models were trained on a workstation with NVIDIA GTX3090 GPU and 32 CPU cores. Binary Cross Entropy [@doi:10.1007/s10479-005-5724-z] loss is used to train the model. The rectified linear unit (ReLU) non-linearity is used except for the last layer which uses a sigmoid function that  yields the probability for zinc per voxel. A dropout layer (p = 0.1) was used between the 5th and 6th layers. The network was trained using AdaDelta employing a stepped learning rate (lr=0.5, ɣ=0.9), a batch size of 150, and 12 epochs to train.


#### Hyperparameter tuning
We used the ray[tune] library [@doi:10.48550/arXiv.1712.05889] to perform a hyperparameter search choosing 20 different combinations between the following parameters with the best combination of parameters in bold:

- filtersize: **3**,4 (in units of 0.5 Å)
- dropout : **0.1**, 0.2, 0.4, 0.5
- learning rate : **0.5**, 1.0, 2.0
- gamma: 0.5, 0.7, 0.8, **0.9**
- largest dimension 80, 100, **120**

#### Grid Averaging 
The model takes as input a `(8,32,32,32)` tensor and outputs a `(1,32,32,32)` tensor containing the probability density for zinc centered on the Cɑ atom of the input residue. 
Predictions for a complete protein were obtained by voxelizing select residues of the protein (default all cysteines, histidines, aspartates, glutamates) and averaging the boxes using a global grid  (Figure @fig:method B). 98 % of the metal sites in the training data have at least one of those residues closeby wherefore this significant decrease in computational cost seems appropriate for most uses.  The global grid is obtained by computing the bounding box of all points and using a regular spaced (0.5 Å) grid. For each grid point in the global grid the predicted probability maps within 0.25 Å of the grid point are averaged. The search is sped up using the KD-Tree implementation in scipy.[@doi:10.1038/s41592-019-0686-2] 

##### Metal ion placement
The global probability density is used to perform clustering of voxels above a certain probability threshold (default p=0.15, cutoff 7 Å) using AgglomerativeClustering implemented in scikit-learn [@doi:10.48550/arXiv.1201.0490]. For each cluster the weighted average of the voxels in the cluster is computed using the probabilities for each point as the weight. This results in one metal placed per cluster.

##### Visualization 

We make available a command line program and interactive notebook allowing the user to visualize the results. The averaged probability map is stored as a `cube` file. The most likely metal coordinates for use in subsequent processing are stored in a `pdb` file. The command line program uses VMD [@pmid:8744570] to visualize the input protein and the predicted density, for the jupyter notebook 3Dmol.js/py3Dmol [@doi:10.1093/bioinformatics/btu829] is used. 

[@2cbastructurepaper]:  doi:10.1016/0022-2836(92)90531-N



#### Evaluation 

##### Comparison 
In order to standardize the evaluation between different tools, we always used the same test set used for the training of Metal1D and Metal3D. In order to compute standard metrics such as precision and recall, we chose to assess the performance of all assessed tools (Metal1D, Metal3D, BioMetAll, MIB) in a binary fashion. Any prediction within 5 Å of an experimental metal site is counted as true positive (TP). Multiple predictions by the same tool for the same site are counted as 1 TP. Any experimental site that has no predicted metal within 5 Å is counted as false negative (FN). A false positive (FP) prediction is a prediction that is not within 5 Å of a zinc site and also not within 5 Å of any other false positive prediction. If two or more false positive predictions are within 5 Å, they are counted as a single false positive prediction for the same site. In practice we first evaluate the true positive and false negative predictions and remove those from the set of predicted positions. The remaining predictions are all false positives and are clustered using AgglomerativeClustering with a radius of 5 Å. The number of false positives is determined from the number of clusters. Using the binary metric we assessed how good the models are at discovering sites and how much these predictions can be trusted.


In order to assess the quality of the predictions, we additionally compute for all the true positive predictions the mean of the Euclidean distance between the true and predicted site (mean absolute deviation MAD). For Metal1D, MIB, and BioMetAll, MAD was computed for all predictions above the threshold within 5 Å of a true zinc site where  $\sum\text{predicted sites} \geq \sum\text{ TP}$. This was done as some tools predict the same site for different residue combinations and we wanted to assess the general performance for all predicted sites above a certain cutoff and not just for the best predicted site above the cutoff. For Metal3D the weighted average of all voxels above the cutoff was used.


Precision was calculated as 

$$
\text{Precision} =\frac{\#\text{ correct metal sites}}{\#\text{ correct metal sites} + \#\text{ false positive clustered}} = \frac{\text{TP}}{\text{TP}+\text{FP}}
$$

Recall was calculated as 

$$
\text{Recall} =\frac{\#\text{ correct metal sites}}{\#\text{ correct metal sites} + \#\text{ not found metal sites}}  = \frac{\text{TP}}{\text{TP}+\text{FN}}
$$

##### Model assessment Metal3D
To evaluate the trained models we monitored loss and how accurately the model predicts the metal density of the test set. We used a discretized version of the Jaccard index setting each voxel either as 0 (no metal) or 1 (zinc present). We tested multiple different decision boundaries (0.5, 0.6, 0.75, 0.9) and also compared a slightly smaller centered box to remove any spurious density at the box edges, where the model has only incomplete information to make predictions.

The Jaccard index is computed as 
$$
J=\frac{\#\left|V_{p} \cap V_{exp}\right|}{\#\left|V_{p} \cup V_{exp}\right|},
$$
where $Vp$ is the array of voxels with predicted probability above the decision boundary and $V_{exp}$ is the array of voxels with the true metal locations also discretized at the same probability threshold.

##### HCA2 mutants 
The data for human carbonic anhydrase 2 (HCA2) mutants was extracted from refs [@doi:10.1021/bi00255a003; @doi:10.1021/ja00079a046;@doi:10.1021/bi00089a005; @doi:10.1073/pnas.92.11.5017; @doi:10.1021/bi9526692] and the crystal structure 2CBA [@2cbastructurepaper; @doi:10.2210/pdb2CBA/pdb] was used. The zinc was modeled using the zinc cationic dummy model forcefield [@pmid:11106157] and we verified that energy minimization produced the correct coordination environment. The Richardson rotamer library [@pmid:10861930] was used with the EVOLVE-ddG energy function to compute the most stable rotamer for a given mutation with the zinc present. The lowest-energy mutant was used for the prediction of the location of metals using Metal3D.

## Additional information 

### Acknowledgement
Supported by Swiss National Science Foundation Grant Number 200020-185092 with computational resources from the Swiss National Computing Centre CSCS.
We thank the developers of [3Dmol.js](https://3dmol.csb.pitt.edu/) and [manubot](https://manubot.org/) that made the interactive version of this manuscript possible.

### Data availibility
Code and training data are available on [lcbc-epfl/metal-site-prediction](https://github.com/lcbc-epfl/metal-site-prediction) and on Zenodo.

### Conflict of interest
None declared

### Author contributions
S.L.D and A.L designed research, S.L.D, A.L, U.R conceptualized research, S.L.D and A.L developed methodology and software, S.L.D and A.L wrote first draft, S.L.D, A.L, U.R revised and edited draft, U.R supervised research and acquired funding.



## Supplement

### Metal1D 

| Amino acid    | Residue name | Label(s) | Search radius increase (Å) |
|:--------------|:------------:|:--------:|:--------------------------:|
| Alanine       | ALA          | O        | 0                          |
| Arginine      | ARG          | NH1, NH2 | 1.2                        |
| Asparagine    | ASN          | OD1      | 0                          |
| Aspartic acid | ASP          | OD1, OD2 | 1.105                      |
| Cysteine      | CYS          | SG       | 0                          |
| Glutamic acid | GLU          | OE1, OE2 | 1.105                      |
| Glutamine     | GLN          | OE1      | 0                          |
| Glycine       | GLY          | O        | 0                          |
| Histidine     | HIS          | ND1, ND2 | 1.08                       |
| Isoleucine    | ILE          | O        | 0                          |
| Leucine       | LEU          | O        | 0                          |
| Lysine        | LYS          | NZ       | 0                          |
| Methionine    | MET          | SD       | 0                          |
| Phenylalanine | PHE          | O        | 0                          |
| Proline       | PRO          | O        | 0                          |
| Serine        | SER          | OG       | 0                          |
| Threonine     | THR          | OG1      | 0                          |
| Tryptophan    | TRP          | O        | 0                          |
| Tyrosine      | TYR          | OH       | 0                          |
| Valine        | VAL          | OH       | 0                          |

Table: Atoms used as reference points for each amino acid in Metal1D. In the case of amino acids with more than one possible ligand atom, the search radius is enlarged, the increase is computed from the midpoint between all ligating atoms. Typical values computed for structure data files downloaded from the PDB are reported. {#tbl:metal1dreferencepoints tag="S1"}



| Channel name    | Selected atoms                                                                                      |
|:----------------|:----------------------------------------------------------------------------------------------------|
| aromatic        | HIS TRP TYR PHE sidechain without CB                                                                |  
| hydrophobic     | element C                                                                                           |
| occupancy       | all protein heavy atoms                                                                             |
| hbond donor     | (ASN GLN TRP MSE SER THR MET CYS and name ND2 NE2 NE1 SG SE OG OG1) and name N                      |
| hbond acceptor  | (resname ASP GLU HIS SER THR MSE CYS MET and name ND2 NE2 OE1 OE2 OD1 OD2 OG OG1 SE SG) or name O   |
| metalbinding    | (name ND1 NE2 SG OE1 OE2 OD2) or (protein and name O N)                                             |
| positive charge | resname LYS ARG HIS and name NZ NH1 NH2 ND1 NE2 NE                                                  |
| negative charge | resname ASP GLU and name OD1 OD2 OE1 OE2                                                            |

Table: Atom selections used for voxelization of proteins using moleculekit {#tbl:voxelchannels tag="S2"}

![Discretized Jaccard indices using different cutoffs for edge trimming and different probability cutoffs (p(metal)) showing that Metal3D predictions well reproduce the target environments in the test set.](images/jaccard_0.5.jpg){#fig:jaccard tag="S1"}


### Precision/Recall

![Precision recall curve for Metal1D and Metal3D with the probability cutoffs (Metal3D) or thresholds (Metal1D) used in the analysis.](images/precisionrecall_0.5.jpg){#fig:precisionrecall tag="S2"}


### Metal selectivity 

#### Metal1D

![Recall for zinc testset and a 25 randomly drawn structures for other transition, alkali and earth-alkali metals for Metal1D](images/metal1D_metal_selectivity.jpg){#fig:selectivity-metal1d tag="S3"}

![Distance distribution Metal1D](images/model_0.5metal1D_distances_violin.jpg){#fig:selectivity-distance-metal1d tag="S4"}

#### Metal3D


![MAD for Metal3D for all sites with 3+ unique protein ligands in the test set and for the selected structures for the other metals. For each ion the whisker plot indicates the median probability (white dot) and the first quartiles (black box).](images/model_0.5metal3D_distances_violin_0.5.jpg){#fig:selectivity-distance-metal3d tag="S5"}

![Probability distribution for Metal3D on sites with 3+ unique protein ligands in the test set and for the selected structures for the other metals.  For each ion the whisker plot indicates the median probability (white dot) and the first quartiles (black box).](images/probability_violin.jpg){#fig:selectivity-probability-metal3d tag="S6"}

### Comparison 

![MAD only 2 residue coordinated zincs. For each tool the whisker plot indicates the median (white dot) and the first quartiles (black box).](images/mad_violin_0.5_2+.jpg){#fig:madonlyGoodZnmetal3d tag="S7"}

##### False positive Metal3D p=0.9
![False positive for Metal3D at p=0.9 in PDB 4JJJ. A calcium site is misclassified as zinc site.](images/4JJJ_FalsePositive_p=0.9_annotated.png){#fig:4jjjFP tag="S8"}



##### Mean absolute deviations

The measured MAD for all true positive zinc sites for each tool in the database. 

| tool            |  MAD [Å]             |   median  [Å]|
|:----------------|----------------------:|---------:|
| BioMetAll c=0.5 |          2.72 ± 1.33  |  2.86    |
| MIB t=1.25      |          1.13 ± 1.24  |  0.60    |
| MIB t=1.9       |          0.77 ± 1.09  | **0.44**    |
| Metal1D t=0.5   |          2.07 ± 1.33  |  2.06    |
| Metal1D t=0.75  |          2.19 ± 1.26  |  2.12    |
| Metal3D p=0.25  |          0.74 ± 0.66  |  0.61    |
| Metal3D p=0.5   |          0.73 ± 0.66  |  0.54    |
| Metal3D p=0.75  |          0.71 ± 0.64  |  0.51    |
| Metal3D p=0.9   |          **0.70 ± 0.64**  |  0.52    |

Table: MAD and median on zinc test set for all zinc sites (n=189). {#tbl:compmadmediantest tag="S3"}


| tool            |   MAD [Å]         | median [Å]|
|:----------------|----------------------:|-------:|
| BioMetAll c=0.5 |          2.68 ± 1.33  |  2.84  |
| MIB t=1.25      |          1.09 ± 1.21  |  0.60  |
| MIB t=1.9       |          0.77 ± 1.09  |  **0.44**  |
| Metal1D t=0.5   |          1.97 ± 1.29  |  1.99  |
| Metal1D t=0.75  |          2.06 ± 1.24  |  2.09  |
| Metal3D p=0.25  |          0.69 ± 0.58  |  0.56  |
| Metal3D p=0.5   |          0.69 ± 0.59  |  0.54  |
| Metal3D p=0.75  |          0.71 ± 0.64  |  0.51  |
| Metal3D p=0.9   |          **0.70 ± 0.64**  |  0.52  |

Table: MAD and median on zinc test set for zinc sites with at least 2 protein ligands (n=133). {#tbl:compmadmediantestgoodonly tag="S4"}



## References {.page_break_before}

<!-- Explicitly insert bibliography here -->
<div id="refs"></div>
