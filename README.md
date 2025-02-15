## WELCOME

This repository includes all scripts and custom libraries used in the work reported in [Rino-Silvestre, J. et al (2022)](https://ui.adsabs.harvard.edu/abs/2022arXiv221015400R/abstract) under the folder "Scripts", as well as all SKIRT simulations used under the folder "Sources". All of our results can be reproduced with the files present in those two folders.

Should the user have no desire to test the training script ("VDAE_Model_Train_Custom_Loss_v1.1.R"), the trained VDAE model, which results are described in the manuscript, is also provided in the folder "Training_Results/VDAE_Model".

## INSTRUCTIONS 

To reproduce the results reported in the manuscript:

	1) Download the EmulART folder tree.
	
	2) Adapt the path variables within the scripts to fit with your system.
	
	3) Train the VDAE model by running the script "VDAE_Model_Train_Custom_Loss_v1.1.R".	
	
	4) Create emulations by running the script "Encoder_Sampler_INLA_Decoder_v1.2.R". This script is set to use all 1e6 photon packet simulations, within the "Sources" folder tree, as inputs. This means it will perform 18 emulations (3 optical thickness values, 2 pov, 3 spatial sampling amounts), if you do not wish to do all these in the same run you must edit out the undesired values of the variable arrays "opt_depth_arr", "pov_str_arr", "pov_num_arr" and "naturals" (more information of the meaning of these can be found in the manuscript). If you wish to use another photon packet amount you must edit the value of the variable "counts_exp" (it can take values 4, 5, 6, 7 and 8). 
	
	5) Create statistics and residuals FITS files by running the script "Very_Statte.R".
	
	6) Create PDFs with plots regarding the residuals by running the script "Very_Plotter_v1.3.R".
	
The script "testing_norm.R" can be used at any point after step 2), and "Latent_Model_Features_Stats.R" can be used at any point after step 3).

## SCRIPTS 

Encoder_Sampler_INLA_Decoder_v1.2.R -> Emulate a SKIRT simulation using the pipeline described in the report. It assumes ‘VDAE_Model_Train_Custom_Loss_v1.1.R’ was run beforehand, and that VDAE model exists in a particular file directory. If you wish to use a different model you should then change the segment of this script concerned with loading the neural network architecture, as well as dependent variables. 

Latent_Model_Features_Stats.R -> Script used to calculate the Pearson correlation coefficients between the latent features and to create the latent features corner plot presented in the manuscript.

testing_norm.R -> Script used to validate the photon flux regularization procedure described in the manuscript.

VDAE_Model_Train_Custom_Loss_v1.1.R -> Train a VDAE model with the architecture and datasets described in the report. If you wish to train a model on different data you should change the initial part of the script concerned with loading data files, as well as dependent variables.

Very_Plotter_v1.3.R -> Creates several PDF files with plots regarding emulations residuals.

Very_Statte.R -> Creates SEDs and Residuals files, calculates and produces histograms for metrics to evaluate the quality of the emulations. It assumes ‘Encoder_Sampler_INLA_Decoder_v1.2.R’ was run beforehand and that a particular folder structure exists. 

## CUSTOM LIBRARIES 

List of libraries:

DAT_lib.R -> Library of functions used in the scripts above to perform different types of operations, from binning to re-scaling, on data.

SA_lib.R -> Library of functions used in the scripts above to identify different types of spaxels and create index maps based on those.

### END NOTES 

For an easier use of these scripts download the whole folder tree into your Desktop. If you prefer to place the folder tree somewhere else, or to arrange the folders in a different manner, you’ll have to adapt all “*path” variables defined in each script so as to reflect those changes.

More information regarding each script can be found in comments within them. 

### CONTACT 

If you have any questions feel free to contact us (please state “EmulART - Repo Questions” within the subject field) via: joao.silvestre@tecnico.ulisboa.pt
