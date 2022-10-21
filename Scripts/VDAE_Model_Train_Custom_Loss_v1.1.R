rm(list = ls())
options(digits = 2)
require(keras)
require(FITSio)
require(jjb)
source("~/Desktop/SKIRT-VDAE-INLA/Scripts/SA_lib.R")

################################################################################
# This command ensures the proper functioning of the vae structure in tf > 2.0 #
################################################################################
tensorflow::tf$compat$v1$disable_eager_execution()                             #
################################################################################

input_1e8_path <- paste0("~/Desktop/SKIRT-VDAE-INLA/Sources/sphShell/FITS/t9-0",
                         ".05/face/sphShell_t9-0.05_aniso_pp1e8_i00_secondaryd",
                         "irect.fits")
input_1e7_path <- paste0("~/Desktop/SKIRT-VDAE-INLA/Sources/sphShell/FITS/t9-0",
                         ".05/face/sphShell_t9-0.05_aniso_pp1e7_i00_secondaryd",
                         "irect.fits")
input_1e6_path <- paste0("~/Desktop/SKIRT-VDAE-INLA/Sources/sphShell/FITS/t9-0",
                         ".05/face/sphShell_t9-0.05_aniso_pp1e6_i00_secondaryd",
                         "irect.fits")
input_1e5_path <- paste0("~/Desktop/SKIRT-VDAE-INLA/Sources/sphShell/FITS/t9-0",
                         ".05/face/sphShell_t9-0.05_aniso_pp1e5_i00_secondaryd",
                         "irect.fits")
input_1e4_path <- paste0("~/Desktop/SKIRT-VDAE-INLA/Sources/sphShell/FITS/t9-0",
                         ".05/face/sphShell_t9-0.05_aniso_pp1e4_i00_secondaryd",
                         "irect.fits")
target_path <- "~/Desktop/SKIRT-VDAE-INLA/Training_Results/VDAE_Model/"
scale_path <- "~/Desktop/SKIRT-VDAE-INLA/Sources/Ref/Ref_Scale.dat"

##########################################################
################### Loading input files ##################
##########################################################
input_1e8 <- readFITS(input_1e8_path)                    #
input_1e7 <- readFITS(input_1e7_path)                    #
input_1e6 <- readFITS(input_1e6_path)                    #
input_1e5 <- readFITS(input_1e5_path)                    #
input_1e4 <- readFITS(input_1e4_path)                    #
                                                         #
dims <-  dim(input_1e8$imDat)                            #
first_useful_feature <- 40                               #
features <- first_useful_feature:dims[3]                 #
                                                         #
data_1e8 <- input_1e8$imDat[, , features]                #
data_1e7 <- input_1e7$imDat[, , features]                #
data_1e6 <- input_1e6$imDat[, , features]                #
data_1e5 <- input_1e5$imDat[, , features]                #
data_1e4 <- input_1e4$imDat[, , features]                #
                                                         #
rm(input_1e8, input_1e7, input_1e6, input_1e5, input_1e4)#
##########################################################

#########################################################################
############# Load um scale of the file this is thought for #############
#########################################################################
myScale <- read.csv(scale_path, header = FALSE, sep = " ")[1:dims[3], 1]#
myScale <- myScale[features]                                            #
#########################################################################

##################################################
# Some variables that will be handy further down #
##################################################
feature_size <- length(features)                 #
spatial_size <- dims[1] * dims[2]                #
fraction <- 5/6                                  #
##################################################

###############################################
############### Re-scaling data ###############
###############################################
TP_1e8 <- sum(data_1e8, na.rm = TRUE)         #
TP_1e7 <- sum(data_1e7, na.rm = TRUE)         #
TP_1e6 <- sum(data_1e6, na.rm = TRUE)         #
TP_1e5 <- sum(data_1e5, na.rm = TRUE)         #
TP_1e4 <- sum(data_1e4, na.rm = TRUE)         #
                                              #
data_1e7 <- data_1e7/TP_1e7 * TP_1e8 * 1/10   #
data_1e6 <- data_1e6/TP_1e6 * TP_1e8 * 1/100  #
data_1e5 <- data_1e5/TP_1e5 * TP_1e8 * 1/1000 #
data_1e4 <- data_1e4/TP_1e4 * TP_1e8 * 1/10000#
###############################################

#########################################################################
######################## Turning maps into lines ########################
#########################################################################
data_1e8 <- array_reshape(data_1e8, dim = c(spatial_size, feature_size))#
data_1e8[which(is.infinite(data_1e8))] <- NA                            #
data_1e8[which(is.na(data_1e8))] <- 0                                   #
                                                                        #
data_1e7 <- array_reshape(data_1e7, dim = c(spatial_size, feature_size))#
data_1e7[which(is.infinite(data_1e7))] <- NA                            #
data_1e7[which(is.na(data_1e7))] <- 0                                   #
                                                                        #
data_1e6 <- array_reshape(data_1e6, dim = c(spatial_size, feature_size))#
data_1e6[which(is.infinite(data_1e6))] <- NA                            #
data_1e6[which(is.na(data_1e6))] <- 0                                   #
                                                                        #
data_1e5 <- array_reshape(data_1e5, dim = c(spatial_size, feature_size))#
data_1e5[which(is.infinite(data_1e5))] <- NA                            #
data_1e5[which(is.na(data_1e5))] <- 0                                   #
                                                                        #
data_1e4 <- array_reshape(data_1e4, dim = c(spatial_size, feature_size))#
data_1e4[which(is.infinite(data_1e4))] <- NA                            #
data_1e4[which(is.na(data_1e4))] <- 0                                   #
#########################################################################

####################################################################
############## Setting up the input and reference sets #############
####################################################################
prelim_x <- rbind(data_1e8, data_1e7, data_1e6, data_1e5, data_1e4)#
prelim_y <- rbind(data_1e8, data_1e8, data_1e8, data_1e8, data_1e8)#
mult <- dim(prelim_x)[1] / spatial_size                            #
rm(data_1e8, data_1e7, data_1e6, data_1e5, data_1e4)               #
####################################################################

#########################################
###### Map with zero flux spaxels #######
#########################################
null_spx <- null_spaxel_map(prelim_x, 0)#
#########################################

###################################################################
##### To hold positions of spaxels with spectral values =/= 0 #####
###################################################################
non_null_array_indeces <- which(null_spx == FALSE, arr.ind = TRUE)#
non_null_size <- length(non_null_array_indeces)                   #
###################################################################

#####################################################
#### Map with estimate of real zero flux spaxels ####
#####################################################
true_Null <- true_null_spaxel_map(null_spx, 2, mult)#
#####################################################

################################################################
#### To hold positions of spaxels with spectral values == 0 ####
################################################################
null_array_indeces <- which(true_Null == TRUE, arr.ind = TRUE) #
################################################################

#####################################################################
####################### SAMPLING OF TRUE NULLS ######################
############### To be used to avoid overflowing INLA ################
######################## AND SPAXEL SELECTION #######################
##################### To train the DVAE model #######################
#####################################################################
null_amount <- floor(0.1 * length(null_array_indeces))              #
sample_null_indeces <- sample(which(true_Null == TRUE), null_amount)#
input_data_indeces <- c(non_null_array_indeces, sample_null_indeces)#
input_dataset_size <- length(input_data_indeces)                    #
rm(non_null_array_indeces, null_array_indeces, sample_null_indeces) #
                                                                    #
print('Spaxels have been selected.')                                #
#####################################################################

############################################################
############ Setting up train and test data sets ###########
############################################################
train_size <- round(input_dataset_size * fraction)         #
test_size <- input_dataset_size - train_size               #
                                                           #
xtrain <- array(NA, dim = c(train_size, feature_size))     #
ytrain <- array(NA, dim = c(train_size, feature_size))     #
print(paste0("Training set has ", train_size, " spaxels."))#
                                                           #
xtest <- array(NA, dim = c(test_size, feature_size))       #
ytest <- array(NA, dim = c(test_size, feature_size))       #
print(paste0("Test set has ", test_size, " spaxels."))     #
                                                           #
input_Data <- prelim_x[input_data_indeces, ]               #
ref_Data <- prelim_y[input_data_indeces, ]                 #
############################################################

################################################################################
###### Random sampling indexes for construction of training and test sets ######
####### and for later printing some (n) model SED reconstruction results #######
######## Building train set from the sampled indexes of the spaxel set  ########
######## Building test set from the remaining indexes of the spaxel set ########
################################################################################
set.seed(134775813)                                                            #
                                                                               #
n <- 20                                                                        #
                                                                               #
print_index_sample <- sample(1:test_size, n, replace = FALSE)                  #
                                                                               #
train_index_sample <- sample(1:input_dataset_size, train_size, replace = FALSE)#
test_index_sample <- setdiff(1:input_dataset_size, train_index_sample)         #
                                                                               #
xtrain <- input_Data[train_index_sample, ]                                     #
ytrain <- ref_Data[train_index_sample, ]                                       #
xtest <- input_Data[test_index_sample, ]                                       #
ytest <- ref_Data[test_index_sample, ]                                         #
                                                                               #
intYTest <- apply(ytest, 2, sum, na.rm = TRUE)                                 #
                                                                               #
# Residuals' plots won't handle 0's                                            #
intYTest_r <- intYTest                                                         #
intYTest_r[which(intYTest_r == 0)] <- 1                                        #
################################################################################

##################
# Def Grid Param #
##################
bsize <- 32      #
num_epo <- 4500  #
latent_size <- 8 #
max_bias <- 0.95 #
actf <- "selu"   #
optf <- "adam"   #
lf <- "custom"   #
p_red_lr <- 500  #
p_stop <- 3000   #
red_lr_f <- 0.25 #
##################

################################################################################
################## Building folder structure for output files ##################
################################################################################
LS_folder <- as.character(latent_size)                                         #
                                                                               #
if(!dir.exists(paste0(target_path, LS_folder))){                               #
  mkdir(paste0(target_path, LS_folder))                                        #
}                                                                              #
                                                                               #
MB_folder <- as.character(max_bias)                                            #
                                                                               #
if(!dir.exists(paste0(target_path, LS_folder, "/", MB_folder))){               #
  mkdir(paste0(target_path, LS_folder, "/", MB_folder))                        #
}                                                                              #
                                                                               #
p_folder <- as.character(p_red_lr)                                             #
                                                                               #
if(!dir.exists(paste0(target_path, LS_folder, "/", MB_folder, "/", p_folder))){#
  mkdir(paste0(target_path, LS_folder, "/", MB_folder, "/", p_folder))         #
}                                                                              #
                                                                               #
vdae_struct <- '32-16-8-16-32'                                                 #
                                                                               #
run_name <- paste0(vdae_struct, '_ep=', num_epo, '_bias<', max_bias, '_pat=',  #
                   p_red_lr, '_actf=', actf, '_lf=', lf, '_optf=', optf, '_bs=',
                   bsize, '_wZeros')                                           #
complete_path <- paste0(target_path, LS_folder, "/", MB_folder, "/", p_folder, #
                        "/", run_name)                                         #
path_and_name <- paste0(target_path, run_name)                                 #
################################################################################

###########################################################################
######################### Structuring the Network #########################
###########################################################################
enc_input = layer_input(shape = c(feature_size))                          #
                                                                          #
l1 = layer_dense(enc_input, units = 32, activation = actf,                #
                 bias_constrain = constraint_maxnorm(max_bias))           #
l2 = layer_dense(l1, units = 16, activation = actf,                       #
                 bias_constrain = constraint_maxnorm(max_bias))           #
z_mean = layer_dense(l2, latent_size)                                     #
z_log_var = layer_dense(l2, latent_size)                                  #
                                                                          #
encoder_u = keras_model(enc_input, z_mean)                                #
encoder_s = keras_model(enc_input, z_log_var)                             #
                                                                          #
sampling <- function(arg){                                                #
  z_mean <- arg[, 1:(latent_size)]                                        #
  z_log_var <- arg[, (latent_size + 1):(2 * latent_size)]                 #
  epsilon <- k_random_normal(shape = c(k_shape(z_mean)[[1]]), mean=0)     #
  z_mean + k_exp(z_log_var / 2) * epsilon                                 #
}                                                                         #
                                                                          #
z <- layer_concatenate(list(z_mean, z_log_var)) %>% layer_lambda(sampling)#
                                                                          #
dec_l2 = layer_dense(units = 16, activation = actf,                       #
                     bias_constrain = constraint_maxnorm(max_bias))       #
dec_l1 = layer_dense(units = 32, activation = actf,                       #
                     bias_constrain = constraint_maxnorm(max_bias))       #
dec_mean = layer_dense(units = feature_size, activation = "sigmoid")      #
l2_dec = dec_l2(z)                                                        #
l1_dec = dec_l1(l2_dec)                                                   #
mean_dec = dec_mean(l1_dec)                                               #
                                                                          #
dec_input = layer_input(shape = latent_size)                              #
l2_dec_2 = dec_l2(dec_input)                                              #
l1_dec_2 = dec_l1(l2_dec_2)                                               #
mean_dec_2 = dec_mean(l1_dec_2)                                           #
                                                                          #
decoder = keras_model(dec_input, mean_dec_2)                              #
                                                                          #
vdaen = keras_model(enc_input, mean_dec)                                  #
###########################################################################

################################################################################
######################## Customizing VDAE Loss Function ########################
################################################################################
vae_loss <- function(enc_input, mean_dec){                                     #
  xent_loss <- (feature_size / 1.0) *                                          #
    loss_mean_absolute_percentage_error(enc_input, mean_dec)                   #
  kl_loss <- -0.5 * k_mean(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var),#
                           axis = -1L)                                         #
  xent_loss + kl_loss                                                          #
}                                                                              #
################################################################################

#############################################
####### Saving structure info to file #######
#############################################
sink(paste0(complete_path,'_Structure.txt'))#
summary(vdaen)                              #
sink()                                      #
#############################################

################################################################################
################# Setting up parameters and training the model #################
################################################################################
cat(paste0("\nUsing:\n Latent Size -> ", latent_size,                          #
           "\n Activation Function -> ", actf, "\n Optimization Function -> ", #
           optf, "\n Loss Function -> ",lf,"\n Bias -> ", max_bias,            #
           "\n Patience -> ", p_red_lr,"\n"))                                  #
print("Starting model training:")                                              #
                                                                               #
vdaen %>% compile(optimizer = optf, loss = vae_loss)                           #
vdaen %>%                                                                      #
  fit(xtrain, ytrain, epochs = num_epo, batch_size = bsize,                    #
      validation_split = 1/5,                                                  #
      callbacks = c(callback_reduce_lr_on_plateau(monitor = "val_loss",        #
                                                  factor = red_lr_f,           #
                                                  patience = p_red_lr,         #
                                                  min_delta = 1e-3),           #
                    callback_csv_logger(filename = paste0(complete_path,       #
                                                          "_training_log.csv"),#
                                        separator = ";"),                      #
                    callback_terminate_on_naan(),                              #
                    callback_early_stopping(monitor = "val_loss",              #
                                            min_delta = 1e-3,                  #
                                            patience = p_stop)))               #
################################################################################

###############################################################################
######################## Saving the weights and models ########################
###############################################################################
save_model_weights_hdf5(encoder_u, filepath = paste0(path_and_name,           #
                                                   "_ENCODER_U_WEIGHTS.h5"))  #
save_model_weights_hdf5(encoder_s, filepath = paste0(path_and_name,           #
                                                     "_ENCODER_S_WEIGHTS.h5"))#
save_model_weights_hdf5(decoder, filepath = paste0(path_and_name,             #
                                                   "_DECODER_WEIGHTS.h5"))    #
###############################################################################

#########################################################################
####### Running test data from input set through the encoder model ######
#########################################################################
encoded_u_xtest <- encoder_u %>% predict(xtest)                         #
encoded_s_xtest <- encoder_s %>% predict(xtest)                         #
z_layer <- layer_concatenate(list(encoded_u_xtest, encoded_s_xtest)) %>%#
  layer_lambda(sampling)                                                #
                                                                        #
# This step serves merely to transform z from a tensor layer            #
# into an array for which regular R operation are available             #
formatter = keras_model(dec_input, dec_input*1)                         #
encoded_z_xtest = formatter %>% predict(z_layer, steps = 1)             #
                                                                        #
encoded_xtest_df <- data.frame(t(aperm(encoded_z_xtest, c(2, 1))))      #
                                                                        #
if(latent_size == 1){                                                   #
}else{                                                                  #
  switch(log(latent_size,2),                                            #
         colnames(encoded_xtest_df) <- c("Z1", "Z2"),                   #
         colnames(encoded_xtest_df) <- c("Z1", "Z2", "Z3", "Z4"),       #
         colnames(encoded_xtest_df) <- c("Z1", "Z2", "Z3", "Z4",        #
                                         "Z5", "Z6", "Z7", "Z8"))       #
}                                                                       #
                                                                        #
decoded_xtest <- decoder %>% predict(encoded_z_xtest)                   #
intDec <- apply(decoded_xtest, 2, sum, na.rm = TRUE)                    #
                                                                        #
# Residuals plots won't handle 0's                                      #
xtest_res <- abs(decoded_xtest - ytest) / ytest_r * 100                 #
intDec_res <- abs(intDec - intYTest) / intYTest_r * 100                 #
xtest_res[which(is.infinite(xtest_res), arr.ind = TRUE)] <- NA          #
intDec_res[which(is.infinite(intDec_res), arr.ind = TRUE)] <- NA        #
#########################################################################

################################################################################
############## Printing to file some of the reconstructed results ##############
########################### for performance review #############################
################################################################################
pdf(paste0(complete_path, '_Singles.pdf'), width = 16, height = 21)            #
par(mfrow = c(2, 1))                                                           #
par(cex.main = 1.5)                                                            #
                                                                               #
for (i in print_index_sample){                                                 #
                                                                               #
  input_coord <- test_index_sample[i]                                          #
                                                                               #
  array_coord <- input_data_indeces[input_coord]                               #
                                                                               #
  m <- array_coord %% spatial_size                                             #
  mi <- array_coord %/% spatial_size                                           #
  mf <- -1 * (mi + 1) + 9                                                      #
                                                                               #
  matrix_coord <- c((m - 1) %/% dims[2] + 1, (m - 1) %% dims[2] + 1)           #
                                                                               #
  latent_feature_string <- toString(encoded_z_xtest[i, ], sep = ",")           #
                                                                               #
  temp <- log10(c(ytest[i, ], decoded_xtest[i, ]))                             #
                                                                               #
  min_dec <- min(temp, na.rm = TRUE)                                           #
  max_dec <- max(temp, na.rm = TRUE)                                           #
                                                                               #
  if(is.infinite(min_dec)){ min_dec <- -50 }                                   #
  if(is.infinite(max_dec)){ max_dec <- 0 }                                     #
                                                                               #
  matplot(log10(myScale), cbind(log10(decoded_xtest[i, ]), log10(ytest[i, ])), #
          type = "b", pch = 1,                                                 #
          sub = paste0("E-D Prediction (black), True (red)"),                  #
          xlab = "log10(Wavelength (um))", ylab="log10(Power) (Normalized)",   #
          ylim = c(min_dec, max_dec))                                          #
  axis(3, at = log10(myScale), labels = round(myScale))                        #
  legend(2.7, median(temp), legend = c("E-D Pred", "True"), col = c(1, 3),     #
         lty = c(1, 3))                                                        #
  mtext(paste0("Spaxel (r,c): (", matrix_coord[1], ",", matrix_coord[2],       #
               "), from Sim 1e", mf, "."), side = 1, line = -1)                #
  mtext(paste0("Latent Feature's Values: {", latent_feature_string, "}"),      #
        side = 1, line = -2)                                                   #
  title("Predictions and True spaxel SEDs", line = 2.5)                        #
                                                                               #
  temp <- xtest_res[i, ]                                                       #
                                                                               #
  min_res <- min(temp, na.rm = TRUE)                                           #
  max_res <- max(temp, na.rm = TRUE)                                           #
                                                                               #
  if(is.infinite(min_res)){ min_res <- 0 }                                     #
  if(is.infinite(max_res)){ max_res <- 300 }                                   #
                                                                               #
  matplot(log10(myScale), cbind(log10(xtest_res[i, ])), type = "b", pch = 1,   #
          sub = paste0("Prediction residuals."),                               #
          xlab = "log10(Wavelength (um))", ylab="Residuals (%)",               #
          ylim = c(min_res, max_res))                                          #
  axis(3, at = log10(myScale), labels = round(myScale))                        #
  legend(2.7, median(temp), legend = "E-D Abs Res", col = c(1, 3),             #
         lty = c(1, 3))                                                        #
  title("Residuals of Predictions vs True spaxel SEDs", line = 2.5)            #
}                                                                              #
dev.off()                                                                      #
                                                                               #
pdf(paste0(complete_path,'_Integrated.pdf'), width = 16, height = 40)          #
                                                                               #
if(latent_size == 1){                                                          #
  par(mfrow = c(2, 1))                                                         #
}else{                                                                         #
  par(mfrow = c(3, 1))                                                         #
}                                                                              #
par(cex.main = 1.5)                                                            #
                                                                               #
temp <- log10(c(intDec, intYTest))                                             #
                                                                               #
min_int <- min(temp, na.rm = TRUE)                                             #
max_int <- max(temp, na.rm = TRUE)                                             #
                                                                               #
if(is.infinite(min_int)){ min_int <- -50 }                                     #
if(is.infinite(max_int)){ max_int <- 0 }                                       #
                                                                               #
matplot (log10(myScale), cbind(log10(intDec), log10(intYTest)), type = "b",    #
         pch = 1, sub = paste0("E-D Prediction (black), True (red)"),          #
         xlab = "log10(Wavelength (um))", ylab = "log10(Power) (Normalized)",  #
         ylim = c(min_int, max_int))                                           #
axis(3, at = log10(myScale), labels = round(myScale))                          #
legend(2.7, median(temp), legend = c("Prediction","True"), col = c(1, 3),      #
       lty = c(1, 3))                                                          #
title("Predictions and True spatial integration SEDs", line = 2.5)             #
                                                                               #
temp <- intDec_res                                                             #
                                                                               #
min_int_res <- min(temp, na.rm = TRUE)                                         #
max_int_res <- max(temp, na.rm = TRUE)                                         #
                                                                               #
if(is.infinite(min_int_res)){ min_int_res <- 0 }                               #
if(is.infinite(max_int_res)){ max_int_res <- 300 }                             #
                                                                               #
matplot(log10(myScale), cbind(intDec_res),                                     #
        type = "b", pch = 1, sub = paste0("Prediction residuals."),            #
        xlab = "log10(Wavelength (um))", ylab = "Residuals (%)",               #
        ylim = c(min_int_res, max_int_res))                                    #
axis(3, at = log10(myScale), labels = round(myScale))                          #
legend(2.7, median(temp), legend = "E-D Abs Res", col = 1, lty = 1)            #
title("Residuals of Prediction vs True integration of testset SEDs", line = 2.5)
                                                                               #
if(latent_size == 1){                                                          #
}else{                                                                         #
  switch(log(latent_size, 2),                                                  #
         pairs(~Z1 + Z2, data = encoded_xtest_df, upper.panel = NULL),         #
         pairs(~Z1 + Z2 + Z3 + Z4, data = encoded_xtest_df,                    #
               upper.panel = NULL),                                            #
         pairs(~Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8,                         #
               data = encoded_xtest_df, upper.panel = NULL))                   #
}                                                                              #
dev.off()                                                                      #
                                                                               #
pcc_m < matrix(NA, nrow = latent_size, ncol = latent_size)                     #
                                                                               #
for(i in 1:latent_size){                                                       #
  for(j in  i:latent_size){                                                    #
    if(i == j){                                                                #
      pcc_m[i, j] <- 1                                                         #
    }else pcc_m[i, j] <- pearsonCC(as.matrix(encoded_xtest_df[i]),             #
                                   as.matrix(encoded_xtest_df[j]))             #
  }                                                                            #
}                                                                              #
                                                                               #
write.table(pccm, file = paste0(complete_path, "pearsonCC_table.csv"))         #
################################################################################