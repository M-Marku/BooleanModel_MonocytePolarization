## Brief description
The tumour microenvironment is the surrounding of a tumour, including blood vessels, fibroblasts, signaling molecules, the extracellular matrix and immune cells, especially neutrophils and monocyte-derived macrophages. In a tumour setting, macrophages encompass a spectrum between a tumour-suppressive (M1) or tumour-promoting (M2) state. The biology of macrophages found in tumours (Tumour Associated Macrophages) remains unclear, but understanding their impact on tumour progression is highly important. In this paper, we perform a comprehensive analysis of a macrophage polarization network, following two lines of enquiry: (i) we reconstruct the macrophage polarization network based on literature, extending it to include important stimuli in a tumour setting, and (ii) we build a dynamical model able to reproduce macrophage polarization in the presence of different stimuli, including the contact with cancer cells. Our simulations recapitulate the documented macrophage phenotypes and their dependencies on specific receptors and transcription factors, while also unravelling the formation of a special type of tumour associated macrophages in an in vitro model of chronic lymphocytic leukaemia. \
For a full and more detailed description of the model, simulations, and other analysis, please refer to Marku et al (2020).\
Marku, M., Verstraete, N., Raynal, F., Madrid-Mencía, M., Domagala, M., Fournié, J.J., Ysebaert, L., Poupot, M. and Pancaldi, V., 2020. Insights on TAM Formation from a Boolean Model of Macrophage Polarization Based on In Vitro Studies. Cancers, 12(12), p.3664.

## Steps
This study focuses on detecting macrophage differentiation states (noted as M1, M2 and TAM) by applying a Boolean model on 
a literature and data-driven regulatory network. We organize the work into 3 main steps:
1. Building the macrophage polarization regulatory network: We start from Palma et al (2018) and we follow a literature based approach to extend the network. Microarray data were used to perform TF estimation and gene expression analysis to validate the extension.
2. Studying the dynamical behaviour of the regulatory system: We build a Boolean model to study the system behaviour and identify the main macrophage states (noted M1, M2 and NLC/TAM), by performing a supervised (literature based) and unsupervised (clustering the attractors according to their Jaccard index) method. 
3. Model validation: We perform several simulations mimicking knock-out experiments and different extracellular environments to validate the accurancy of the model. 


## Code description
The numerical simulations on the Boolean model of macrophage polarization are performed using the R library BoolNet. The R script run_macrophage.R in this repo contains all the steps of the simulations:
- Identification of all fized point attractors of the system;
- Attractor categorization.
- Estimation of Jaccard-Needham index and clustering.
