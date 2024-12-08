# Zero_inflated_JSM
Code for implementation of Zero-inflated Joint species models
The script JSM_zibeta_TML.R contains three parts. 

(i) I simulate data of species presence and cover across an environmental gradient following a Hierarchical Model of Species Communities (HMSC). Please see https://onlinelibrary.wiley.com/doi/10.1111/ele.12757, for a full description of this model. In here we are not including latent factors that structure the species correlation error structure (as it is included in the Hmsc package).

(ii) I create a stan model to fit this data (i.e., presence of species and cover across an environmental gradient). The Hmsc package in R allow to HMSC models, but does not have zero-inflated beta distribution, which may be interesting for modeling  assembly in plant communities. See https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/1365-2745.14168 or https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/1365-2664.14529 for an implementation of this model. 
  *Some parts of HMSC in stan were based on scripts of Dr. Juan Manuel Morales (https://orcid.org/0000-0001-7269-7490) please see jmmorales github account for other examples.
  
(iii) I fit simulated data in STAN and evaluate if I could recover simualted parameters (i.e., model performance in the best-case scenario). 

To fit your own data you only need to copy the stan file that is generated in part ii, and load the data following the same structure as in the stan_data list. You can modify this code to include hierarhical effects, spatial autocorrelation etc.
