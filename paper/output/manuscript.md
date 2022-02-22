---
author-meta:
- "Simon L. D\xFCrr"
- Andrea Levy
- Ursula Rothlisberger
bibliography: []
date-meta: '2022-02-22'
header-includes: "<!--\nManubot generated metadata rendered from header-includes-template.html.\nSuggest improvements at https://github.com/manubot/manubot/blob/master/manubot/process/header-includes-template.html\n-->\n<meta name=\"dc.format\" content=\"text/html\" />\n<meta name=\"dc.title\" content=\"Accurate prediction of metalsites using deep learning\" />\n<meta name=\"citation_title\" content=\"Accurate prediction of metalsites using deep learning\" />\n<meta property=\"og:title\" content=\"Accurate prediction of metalsites using deep learning\" />\n<meta property=\"twitter:title\" content=\"Accurate prediction of metalsites using deep learning\" />\n<meta name=\"dc.date\" content=\"2022-02-22\" />\n<meta name=\"citation_publication_date\" content=\"2022-02-22\" />\n<meta name=\"dc.language\" content=\"en-US\" />\n<meta name=\"citation_language\" content=\"en-US\" />\n<meta name=\"dc.relation.ispartof\" content=\"Manubot\" />\n<meta name=\"dc.publisher\" content=\"Manubot\" />\n<meta name=\"citation_journal_title\" content=\"Manubot\" />\n<meta name=\"citation_technical_report_institution\" content=\"Manubot\" />\n<meta name=\"citation_author\" content=\"Simon L. D\xFCrr\" />\n<meta name=\"citation_author_institution\" content=\"Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland\" />\n<meta name=\"citation_author_orcid\" content=\"0000-0002-4304-8106\" />\n<meta name=\"twitter:creator\" content=\"@simonduerr\" />\n<meta name=\"citation_author\" content=\"Andrea Levy\" />\n<meta name=\"citation_author_institution\" content=\"Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland\" />\n<meta name=\"citation_author_orcid\" content=\"0000-0003-1255-859X\" />\n<meta name=\"citation_author\" content=\"Ursula Rothlisberger\" />\n<meta name=\"citation_author_institution\" content=\"Laboratory of Computational Chemistry and Biochemistry, Institute of Chemical Sciences and Engineering, Swiss Federal Institute of Technology (EPFL) CH-1015 Lausanne, Switzerland\" />\n<meta name=\"citation_author_orcid\" content=\"0000\u20100002\u20101704\u20108591\" />\n<meta property=\"og:type\" content=\"article\" />\n<meta property=\"twitter:card\" content=\"summary_large_image\" />\n<link rel=\"icon\" type=\"image/png\" sizes=\"192x192\" href=\"https://manubot.org/favicon-192x192.png\" />\n<link rel=\"mask-icon\" href=\"https://manubot.org/safari-pinned-tab.svg\" color=\"#ad1457\" />\n<meta name=\"theme-color\" content=\"#ad1457\" />\n<!-- end Manubot generated metadata -->"
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
on February 22, 2022.
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

Metalloproteins protein many biological functions [Percora review]? 

Understanding where metals bind in biology is related to health [@doi:10.1038/s42256-019-0119-z], biocatalyis [Kuhlman, Hilvert] and PPIs [Tezcan]. 

Starting with pioneering studies in the 1990 s,[@doi:10.1126/science.8346440;@doi:10.1146/annurev.biochem.68.1.779]there have been notable successes inthe de novo design of functional metalloproteins, which arepredominantly based on four-helix bundle anda-helicalcoiled-coiled motifs with readily parametrizable structure [@doi:10.1002/anie.202009226]. 


### Deep learning on proteins
Torng/Shroff 3DCNN bio stuff 
Ananad DeepRank
Correia surface studies Nat Methods 


### Enzymes

### Interfaces
It has been hypothesized that some modern metalloproteins may have emerged through the metal-nucleated oligomerization of small peptides or protein domains, followed by the evolution of the resulting assemblies into stable, functional architecture.


Metal-Templated Interface Redesign (MeTIR) ). These strategies, inspired by both the proposed evolutionary roles of metals and their prevalence in natural PPIs, take advantage of the favorable properties of metal coordination (bonding strength, directionality, and reversibility) to guide protein self-assembly with minimal design and engineering

In order to circumvent the complexity of constructing extensive noncovalent interfaces, which are typically involved in natural PPIs

### Existing approaches

Computational predictors of metal-binding sites built on sequence analyses are mostly based on scanning the sequence of a target protein to identify those regions where amino acid patterns match a metal-binding site fingerprint. For zinc binding sites often two histidine spaced by one residue which allows to readily detect the motif are easy to detect. These predictors yield the identites of the coordinating residues. 

Structural detectors often used distance features to identify sites based on statistical mining in the protein databank.  Predictors trained like this can identify highly preorganized motifs (e.g 4x Cys in close spatial proximity) but are often not very good at identifying weakly preorganized motifs. Current state of the art predictors (MIB, BioMetAll) use fragments/homology to predict the location of the metal or backbone preorganization predicting an approximate position of the metal that is less sensitive to the exact side chain geometry thus affording higher sensitivity to detect metal sites(BioMetAll). MIB [@doi:10.1021/acs.jcim.6b00407] uses the fragmentation transformation method to search for homologus sites in its database


In our work we develop two new predictors primarily intended for zinc binding sites - Metal1D and Metal3D that are more accurate and sensitive than existing approaches in predicting metal ion binding sites. We evaluate their capability with respect to 



## Materials and Methods

### Dataset

### Metal 1D
![Workflow of Metal1D](images/metal1D_scheme_large.png)

The statistical analysis performed is focuses on the LINK records present in deposited PDB structures. In particular, this section of a PDB file specifies the connectivity between zinc (or any other metal ion) and the amino acids of the protein, and each LINK record specifies one linkage. After having extracted the LINK records, data were filtered, only considering X-ray structures with resolution better than $2.5 \text{\AA}$. Aiming to obtain a statistically significant sample, for proteins which present more than one deposited structure, only the best-resolution one was used. This selection is based on the MOLECULE record of the deposited structures.
After the statistical filtering, a probability map is generated from the LINK records, extracting for each zinc ion the coordination experimentally observed. In this phase, LINK records involving only a single amino acid are excluded since they probably correspond to weak binding sites. In this phase, also LINK records containing water molecules are excluded. This is done for two reasons: in the first place, because it would be needed to restrict only to higher-resolution structure to have a reliable probability map including waters, which would drastically reduce the pool of structures used. Moreover, in the predicting phase, it would be needed to have protein structures with waters, which is not the general case due to experimental limitations and to difficulties to place water molecules a posteriori. Even if many tools are available, they are not able to handle the presence of metal ions. Due to these limitations, water molecules are considered only implicitly in this method.
<!---
Not sure if it is the best way to discuss why water molecules were not considered. My be a good idea to have a separate section in the methods, describing problems related to waters molecules (which would motivate why in both methods waters are considered only implicitly.)
-->

The probability map is used to predict the metal sites for a given protein structure. To make a prediction, each amino acid of the protein is scored. The score is assigned performing a geometrical serch from a reference point, defined as the coordinate of the most probable metal binding atom, within a search radius, considered as twice the typical distance between the metal ion and the binding atom of amino acids in proteins. This quantity is enlarged of an arbitrary factor, in order to be able to take into account deviations from the ideal cooridnation and give more flexibility to the search. In case of amino acids which present more than one atom which typically binds metals, such as Histidine, the mid-point is used and the search radius is enlarged accordingly. 
The score is assigned to each amino acid considering all the other reference points of other amino acids within the search radius, and summing the probabilities in the probability map for coordinations compatible with the one observed. In the ideal case, a score of $1$ would correspond to an amino acid surrounded by all possible coordinating amino acids observed in the probability map. In practice, scores result between $0$ and less than $1$.

Once all amino acids in the chain are scored, site predictions are made grouping the highest-scored amino acids in clusters, based on distance. In practice, highest-scored amino acids are the ones with a score within a given threshold of the highest-scored one. For each cluster, a site prediction is made as a weighted average between the coordinates of the reference point of each amino acid, using as weighting factor the amino acid score. In the case of clusters composed of only one amino acid, to be able anyway to perform a site prediction, a fictitious score equal to the single highest-scored amino acid in the cluster is assigned to the nearby amino acid  (within the search radius) with the highest score. It is important to note that this amino acid was not in the highest-scored ones used to make the clusters, but in this way also for clusters with only one amino acid, a prediciton is made. The potential artefacts introduced by assigning a fictitious score to another amino acid is resolved by a final re-scoring of all the predicted sites.
This final re-scoring also mitigates the errors which can be introduced by calculating the site simply as a weighted average. In particular, a final geometrical search is performed around each predicted site (within $0.6$ the search radius, which as explained is about twice the typical metal-amino acid distance) and a score is now assigned to the site. This score is assigned in the same way as the amino acid scores, based on the probability map, and has the advantage of being able to sort the predicted metal sites based on their probability, according to the method. In this last re-scoring, sites with a probability lower than a certain threshold with respect to the highest-scored one, are excluded, in order to focus only on predictions with the highest probability. 

### Metal 3D
>>>>>>> 38c43e44882c2e35996f298863d80521b91134b0


## Results

### Metal1D

<!---
The starting PDB structures used are the ones in the training set, discuss what's the best way to describe it: repeat what is also said later for Metal3D?  
-->

The statistical analysis performed is focuses on the LINK records present in deposited PDB structures. LINK data were filtered to obtain a statistically significant sample, restricting only to X-ray structures with resolution better than $2.5 \text{\AA}$, and considering only the best-resolution structure for proteins presenting more than one deposited structure. This filtered sample is used to extract a probability map for the coordinating amino acids, which is used to predict metal sites: for a given protein structure, each amino acid of the protein is scored based on the probability map, performing a geometrical serch around the amino acid. For Zn(II), the average LINK length resulted $(2.2 \pm 0.2)\text{\AA}$, and the default search radius was then set to $5.5\text{\AA}$. 

Once all amino acids in the chain are scored, site predictions are made grouping the highest-scored amino acids in clusters, based on distance. For each cluster, a site prediction is made as a weighted average between the coordinates of the reference point of each amino acid, using as weighting factor the amino acid score. A final re-scoring is performed, now assigning a score to the predicted locations for the sites, based on surrounding amino acids, in order to sort the predictions based on the probability, according to the method. This final re-scoring process also mitigates the possible artefacts originated by the prediction phase, excluding sites  which are placed in positions with a low-probable coordinating  environment. 

<!---
Not sure if this  two paragraphs are too long for the results section, I tried to highlight how the method works without focusing on the details (explained in the method section), but may need to make it shorter
-->
<script src="https://gist.github.com/lzhou1110/2a30a81cb8c175514ed627bc18016774.js"></script>

### Metal3D

The model was trained using a train/test/val split based on sequence identity. We used the MMSeqs2 clustered PDB at 30% sequence identity and used the highest resolution structure from each cluster that contained a zinc, did not contain DNA and had resolution <2.5A. In case no structure was found the cluster was discarded. 

The training examples were sampled from the chosen structure by choosing a balanced number of boxes from each protein that contained or did not contain a zinc. Each box was randomly rotated such that the model is insensitive to rotation. 

Metal3D predicts a per residue score that can then be averaged over all residues or used individually (e.g for protein design). 

Hyperparameter tuning using Ray tune


Evaluation of quality of predictions per box using discretized Jaccard Score (similarity of two sets). We noticed that at the edges often spurious density is predicted . 


![Discretized Jaccard index using different edge cutoffs and different probability cutoff](images/jaccard.png){#fig:jaccard}

<script src="https://gist.github.com/lzhou1110/2a30a81cb8c175514ed627bc18016774.js"></script>

#### Selectivity for other metals

![Recall for zinc testset and a 25 randomly drawn structures for other transition, alkali and earth-alkali metals for Metal3D](images/metal3d_metal_selectivity.jpg){#fig:selectivity-metal3d}

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

![Metal3D, AlphaFill and 3RZV zinc positions. Metal3D places the metal with high accuracy even if coordination is not perfectly predicted. Probability map colored at p=0.6](images/3rzv_alphafold_metal3d_alphafill.png){#fig:alphafold}

#### Hidden/transient metalsites

Metallodrugs are in important class of drugs that rely on binding inhibitors to a protein (or DNA). Metal3D can be used to screen the hidden metalloproteome by finding transient metal ion binding sites. 

The site where Rapta binds is detected with p=0.3 but in a high resolution structure without (1KX4) there is a salt bridge with a lysine that might occlude metal detection. One could weight by the rotamer/do MD simulation e.g similar as for cryptic pockets. 

![RAPTA (Ru-arene-phosphaadamantane) inhibitor bound to Nucleosome core particle. Metal3D identifies the two Ruthenium sites with low probabilities.](images/5xf6_rapta_vis.png){#fig:drugdesign}

#### HCA2: case study

HCA2 is the first enzyme where a catalytic zinc was discovered and is therefore one of the best studied metalloenzymes to date with a rich throve of mutational data available. 
On the wildtype enzyme crystal structure (2CBA) Metal3D perfectly recapitulates the location of the active site metal when using a high probability cutoff (p>0.4). The sites predicted with lower probability all look like reasonable transient binding sites at the surface of the protein. 

![Probability evolution in HCA2 for different probability cutoffs A) p=0.1 B) p=0.2 C) p=0.3 D) p=0.4 E) p= 0.75](images/2CBA_probabilities.png){#fig:2cba-probabilites}

We used in silico generated mutants matching mutants in the first and second shell of the active site zinc and probed the effect on the predicted metal probability. For mutants that decrease zinc binding also a drop in probability can be observed that correlates well with the experimentally measured $K_d$. We used a consistent set of $K_d$ values from the literature.

![Probability evolution in HCA2 for different probability cutoffs A) p=0.1 B) p=0.2 C) p=0.3 D) p=0.4 E) p= 0.75](images/kd_vs_p_nolog_newmethod.jpg){@fig:hca-kd}

### Comparison of Metal1D, Metal3D, MIB and BioMetAll

Many metal ion predictors exists that can be subdivided in two categories: binding site predictors and binding location predictors. The former label only the residues binding the ion, the latter also predict a location of the ion. 

In addition to Metal1D and Metal3D we also compared two recent predictors BioMetAll and MIB. MIB uses a fragment method to identify homologus binding sites to the motifs it finds in a given structure and will extract the location of the metal from the homologous structures in its database. The main performance regulator of MIB is the tscore cutoff which is a parameter for the template similiarity with higher values requiring higher similiarity. 
BioMetAll was calibrated on the PDB and places probes on a regular grid at all sites where they find the criteria to be fulfilled. For each collection of probes also a center of the probes is given which we used to assess performance as there is no individual ranking of the probes given by the program. The main parameter for BioMetAll is the cluster cutoff which indicates how many probes in reference to the largest cluster a specific cluster has. We used a cutoff of 0.5 requiring all chosen clusters to have at least 50% of the probes of the most popolous. 

For both tools the recommended settings match the accuracy of Metal3D p=0.75 with a lot more false positives. 
Metal1D offers high detection capabilites but also with a high number of false positives.
While MIB also offers high precision, BioMetAll (using the cluster center) is not very precise with a MAD for correctly identified sites of 2.8 +- XX. Metal1D which identifies more sites than BioMetAll is slightly more precise than BioMetAll. MIB detects less sites but does so with high precision because it can use homologues sites to correctly place the metal ligand. BioMetAll also often provides probes that correctly identify the metal but as there is no ranking of the probes any probe could be closest to the actual location. 

![Comparison of Metal1D, Metal3D, BioMetAll and MIB on the testset used to train Metal1D and Metal3D. Predicted sites are counted if within 5 A of true metal location. False positive probes are clustered and counted once per cluster. ](images/metal3d_biometall_comparison.jpg){@fig:comparison}

![Mean absolute deviation of predicted location of zincs using Metal1D, Metal3D, BioMetAll and MIB on the testset used to train Metal1D and Metal3D for all correctly identified sites. ](images/mad_violin.jpg){@fig:distances-testset-Metal3D}



## Discussion

3D CNN model accuracy. Recent work EquiDock uses no sidechains at all, gets ligand RMSD where only 25% are under a 2 threshold, https://arxiv.org/pdf/2202.05146.pdf Mean RMSD 8.3 A, Centroid 42.4


MIB discards all sites with less than 2 coordination partners so it will not be able to identify labile binding sites. 

Discuss 4L99 which is one of the FN for the model. Here are Lysine next to the zinc and therefore probability is low.

## Conclusion

## Supplement

### Metal site detection using Metal1D

Maybe ROC curve for Metal1D and Metal3D here


### Metal selectivity Metal1D

![Recall for zinc testset and a 25 randomly drawn structures for other transition, alkali and earth-alkali metals for Metal3D](images/metal3d_metal_selectivity.jpg){#fig:selectivity-metal3d}

![Distance distribution Metal1D](images/selectivity_metal1D_distances_violin.jpg){#fig:selectivity-distance-metal1d}

![Distance distribution Metal3D](images/distances_violin.jpg){#fig:selectivity-distance-metal3d}

![Probability distribution](images/probability_violin.jpg){#fig:selectivity-probability-metal3d}

### Comparison 

![MAD only 2 residue coordinated zincs](images/mad_violin.jpg){#fig:madonlyGoodZnmetal3d}

![Only 2 residue coordinated zincs](images/metal3d_biometall_comparison_goodZNonly.jpg){#fig:tpfpfnonlyGoodZn}

## References {.page_break_before}

<!-- Explicitly insert bibliography here -->
<div id="refs"></div>
