## Brief description
The Tumour Microenvironment (TME) is the collection of cells in and surrounding cancer cells in a tumour including a variety of immune cells. It often features considerable amounts of myeloid cells, especially neutrophils and monocyte-derived macrophages. In a tumour setting, macrophages cover a spectrum between a tumour suppressive (M1) or tumour promoting (M2) state. The biology of macrophages found in tumours (Tumour Associated Macrophages, TAMs) remains unclear, but understanding their impact on tumour progression is highly important.\
Many factors underlie the large variability in patients’ immune-infiltration patterns and their responses to therapy. Previous work by the group and others highlighted the extent of individual variability in haematopoietic cells in healthy individuals, and currently the group is investigating the relationship between this variability and chromatin 3D structure.

## Steps
This study focuses on detecting macrophage differentiation states (noted as M1, M2 and TAM) by applying a Boolean model on 
a literature and data-driven regulatory network. We organize the work into 3 main steps:
1. Building the macrophage polarization regulatory network: We start from Palma et al (2018) and we follow a literature based approach to extend the network. Microarray data were used to perform TF estimation and gene expression analysis to validate the extension.
2. Studying the dynamical behaviour of the regulatory system: We build a Boolean model to study the system behaviour and identify the main macrophage states (noted M1, M2 and NLC/TAM), by performing a supervised (literature based) and unsupervised (clustering the attractors according to their Jaccard index) method. 
3. Model validation: We perform several simulations mimicking knock-out experiments and different extracellular environments to validate the accurancy of the model. 

## Full text availability
For a full and more detailed description of the model, simulations, and other analysis, please refer to Marku et al (2020).\
Marku, M., Verstraete, N., Raynal, F., Madrid-Mencía, M., Domagala, M., Fournié, J.J., Ysebaert, L., Poupot, M. and Pancaldi, V., 2020. Insights on TAM Formation from a Boolean Model of Macrophage Polarization Based on In Vitro Studies. Cancers, 12(12), p.3664.

## Code description
The numerical simulations on the Boolean model of macrophage polarization are performed using the R library BoolNet. The R script run_macrophage.R in this repo contains all the steps of the simulations:
- Identification of all fized point attractors of the system;
- Attractor categorization.
- Estimation of Jaccard-Needham index and clusting.
