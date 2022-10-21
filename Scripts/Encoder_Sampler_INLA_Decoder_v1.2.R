rm(list = ls())
options(digits = 16)
require(keras)
require(INLA)
require(sp)
require(fields)
require(spam)
require(ggplot2)
require(viridis)
require(yaml)
require(colorspace)
require(spatstat)
require(data.table)
require(gtools)
require(FITSio)
require(IDPmisc)
require(jjb)
set.seed(134775813)
source("/media/joaomfras/TOSHIBA EXT/EmulaRT/Scripts/SA_lib.R")
source("/media/joaomfras/TOSHIBA EXT/EmulaRT/Scripts/DAT_lib.R")

counts_exp_ref <- 8
counts_exp <- 7#Done: 6

# This parameter should reflect the ratio between the amount of photons in the
# input and in the reference file
# Example: input = 10^6, reference = 10^8 -> energy_scaling = 10^-2
energy_scaling <- 10^(counts_exp - counts_exp_ref)

opt_depth_arr <- c("0.05", "0.1", "1.0") 
pov_str_arr <- c("face", "edge")
pov_num_arr <- c("00", "90")
naturals <- c(2, 3, 5)
case_total <- length(opt_depth_arr) * length(pov_str_arr) * length(naturals)

time_array <- array(NA, dim = c(case_total, 2), 
                    dimnames = list(NULL, c("Case", "Time (min)")))

main_folder <- "/media/joaomfras/TOSHIBA EXT/EmulaRT/"
vdae_model_folder <- paste0(main_folder, "Training_Results/VDAE_Model/")
input_path <- paste0(main_folder, "Sources/FITS/")
scale_path <- paste0(main_folder, "Sources/Ref_Scale.dat")
results_folder <- paste0(main_folder, "Results/")
if(!dir.exists(results_folder)){
  mkdir(results_folder)
}

################################################################################
# This command ensures the proper functioning of the vae structure in tf > 2.0 #
################################################################################
tensorflow::tf$compat$v1$disable_eager_execution()                             #
################################################################################

################################################################################
################ Paths to the trained VDAE architecture files ##################
################################################################################
encoder_wu_path <- paste0(vdae_model_folder, "32-16-8-16-32_ep=4341_bias_0.95_",
                          "pat=500_actf=selu_lf=custom_optf=adam_bs=32_wZeros_",
                          "ENCODER_U_WEIGHTS.h5")                              #
encoder_ws_path <- paste0(vdae_model_folder, "32-16-8-16-32_ep=4341_bias_0.95_",
                          "pat=500_actf=selu_lf=custom_optf=adam_bs=32_wZeros_",
                          "ENCODER_S_WEIGHTS.h5")                              #
decoder_w_path <- paste0(vdae_model_folder, "32-16-8-16-32_ep=4341_bias_0.95_",#
                         "pat=500_actf=selu_lf=custom_optf=adam_bs=32_wZeros_",#
                         "DECODER_WEIGHTS.h5")                                 #
################################################################################

#############################################################################
########################## Structuring the Network ##########################
#############################################################################
feature_size <- 64                                                          #
latent_size <- 8                                                            #
max_bias <- 0.95                                                            #
actf <- "selu"                                                              #
ae_struct <- '32-16-8-16-32'                                                #
                                                                            #
enc_input = layer_input(shape = feature_size)                               #
                                                                            #
l1 = layer_dense(enc_input, units = 32, activation = actf,                  #
                 bias_constrain = constraint_maxnorm(max_bias))             #
l2 = layer_dense(l1, units = 16, activation = actf,                         #
                 bias_constrain = constraint_maxnorm(max_bias))             #
z_mean = layer_dense(l2, latent_size,                                       #
                     bias_constrain = constraint_maxnorm(max_bias))         #
z_log_var = layer_dense(l2, latent_size,                                    #
                        bias_constrain = constraint_maxnorm(max_bias))      #
                                                                            # 
encoder_mean = keras_model(enc_input, z_mean)                               #
encoder_var = keras_model(enc_input, z_log_var)                             #
                                                                            #
sampling <- function(arg){                                                  #
  z_mean <- arg[, 1:(latent_size)]                                          #
  z_log_var <- arg[, (latent_size + 1):(2 * latent_size)]                   #
  epsilon <- k_random_normal(shape = c(k_shape(z_mean)[[1]]), mean=0)       #
  z_mean + k_exp(z_log_var / 2) * epsilon                                   #
}                                                                           #
                                                                            #
dec_input = layer_input(shape = latent_size)                                #
dec_l2 = layer_dense(dec_input, units = 16, activation = actf,              #
                     bias_constrain = constraint_maxnorm(max_bias))         #
dec_l1 = layer_dense(dec_l2, units = 32, activation = actf,                 #
                     bias_constrain = constraint_maxnorm(max_bias))         #
dec_mean = layer_dense(dec_l1, units = feature_size, activation = "sigmoid")#
                                                                            #
decoder = keras_model(dec_input, dec_mean)                                  #
                                                                            #
load_model_weights_hdf5(encoder_mean, encoder_wu_path)                      #
load_model_weights_hdf5(encoder_var, encoder_ws_path)                       #
load_model_weights_hdf5(decoder, decoder_w_path)                            #
print("VDAE model is set. Weights have been loaded.")                       #
#############################################################################

case_n <- 1
for(opt_depth in opt_depth_arr){
  for(pov_str in pov_str_arr){

    pov_num <- pov_num_arr[which(pov_str_arr == pov_str, arr.ind = TRUE)]

    for(nat in naturals){
      SKIP <- nat^2
      non_null_samp_rat <- 1 / SKIP
      samp_str <- paste0(round(non_null_samp_rat * 100), "p100")
      
      null_skip <- 10

      ########################################################################
      ######## Defining names and folders for input and output files #########
      ########################################################################
      depth_path <- paste0("t9-", opt_depth, "/")                            #
      pov_path <- paste0(pov_str, "/")                                       #
      case_path <- paste0("1e", counts_exp, "-", samp_str, "/")              #
      input_name <- paste0("sphShell_t9-", opt_depth, "_aniso_pp1e",         #
                           counts_exp, "_i", pov_num, "_secondarydirect.fits")
      input_file <- paste0(input_path, depth_path, pov_path, input_name)     #
                                                                             #
      optD_folder <- paste0(results_folder, "reg-grid_", depth_path)         #
      if(!dir.exists(optD_folder)){                                          #
        mkdir(optD_folder)                                                   #
      }                                                                      #
                                                                             #
      pov_folder <- paste0(optD_folder, pov_path)                            #
      if(!dir.exists(pov_folder)){                                           #
        mkdir(pov_folder)                                                    #
      }                                                                      #
                                                                             #
      case_folder <- paste0(pov_folder, case_path)                           #
      if(!dir.exists(case_folder)){                                          #
        mkdir(case_folder)                                                   #
      }                                                                      #
                                                                             #
      emul_folder <- paste0(case_folder, "Emulation/")                       #
      if(!dir.exists(emul_folder)){                                          #
        mkdir(emul_folder)                                                   #
      }                                                                      #
                                                                             #
      input_name <- strsplit(input_name,'.fits')[[1]][1]                     #
      out_path <- paste0(emul_folder, input_name)                            #
      ########################################################################

      ##########################################################################
      ########################## Loading input files ###########################
      ##########################################################################
      input_fits <- readFITS(input_file)                                       #
      first_useful_feature <- 40                                               #
      dims <- dim(input_fits$imDat)                                            #
      spatial_size <- dims[1] * dims[2]                                        #
      features <- first_useful_feature:dims[3]                                 #
      feat_input_size <- length(features)                                      #
                                                                               #
      if(feature_size != feat_input_size){                                     #
        print(paste0("ERROR: Input number of features is not the expected (",  #
                     feature_size, ") but is instead ", feat_input_size, ". ", #
                     "Change the definition the vector 'features' and re-run."))
        stop()                                                                 #
      }                                                                        #
      rm(feat_input_size)                                                      #
      ##########################################################################

      #########################################################################
      ##### Loading the scale file (the file I thought this for is in um) #####
      #########################################################################
      myScale <- read.csv(scale_path, header = FALSE, sep = " ")[1:dims[3], 1]#
      myScale <- myScale[features]                                            #
      #########################################################################

      #########################################################################
      ## Building a Sampling Grid of spacial locations spaced by SKIP pixels ##
      ########### For better results SKIP should be a perfect square ##########
      #########################################################################
      grid_mask <- NULL                                                       #
      #Building sub-grid of pixels by selecting anchor indexes                #
      #Defining sequence of indexes to select based on how matrix are indexed #
      JUMP <- sqrt(SKIP)                                                      #
      ind_seq <- round(seq(1, dims[2], JUMP)) - 1                             #
      for(ind in ind_seq){                                                    #
        grid_mask <- c(grid_mask, round(seq(1, dims[1], JUMP) + dims[2] * ind))
      }                                                                       #
      rm(ind_seq)                                                             #
      grid_size <- length(grid_mask)                                          #
                                                                              #
      cell_pix <- array(NA, dim = c(grid_size, SKIP))                         #
      #Listing indexes of pixels belonging to the same grid-cell              #
      #Defining the center of grid-cells and updating sub-grid indexes        #
      for(g in 1:grid_size){                                                  #
        pix_cols <- grid_mask[g] + 0:(round(JUMP) - 1)                        #
        pix_rows_inc <- dims[2] * (0:(round(JUMP) - 1))                       #
        col_size <- length(pix_cols)                                          #
        ipix <- NULL                                                          #
                                                                              #
        for(j in 1:col_size){                                                 #
          ipix <- c(ipix, pix_cols + pix_rows_inc[j])                         #
        }                                                                     #
        cell_pix[g, ] <- ipix                                                 #
                                                                              #
        if(SKIP %% 2 == 1){                                                   #
          c_ind <- (1 + SKIP) / 2                                             #
        }else{                                                                #
          c_ind <- (SKIP - round(JUMP)) / 2                                   #
        }                                                                     #
        grid_mask[g] <- cell_pix[g, c_ind]                                    #
      }                                                                       #
      rm(c_ind, ipix, pix_rows_inc, pix_cols, col_size)                       #
      #########################################################################

      ################################################################
      ##### Turning map (without unwanted wavelengths) into line #####
      ################################################################
      input_hold <- input_fits$imDat[, , features] * energy_scaling  #
      input_hold[which(is.infinite(input_hold))] <- 0                #
      input_hold[which(is.null(input_hold))] <- 0                    #
      input_hold[which(is.na(input_hold))] <- 0                      #
      input_temp <- array_reshape(input_hold, dim = c(spatial_size,  #
                                                      feature_size)) #
      ################################################################

      ##################################################################
      # Creating a FITS file to hold a preliminary map of null spaxels #
      ##################################################################
      null_spx <- null_spaxel_map(input_temp, 0)                       #
      null_map <- array_reshape(null_spx, dim = c(dims[1], dims[2]))   #
      temp <- null_spx                                                 #
      temp[which(temp == TRUE, arr.ind = TRUE)] <- 1                   #
      temp <- (temp - 1) * (-1)                                        #
      temp <- array_reshape(temp, dim = c(dims[1], dims[2]))           #
                                                                       #
      writeFITSim(temp, file = paste0(out_path, '-pseudo-Nulls.fits')) #
      rm(input_temp, temp)                                             #
      print('Preliminary Zero Map has been constructed.')              #
      ##################################################################
      
      #################################################################
      ########## Map with estimate of real zero flux spaxels ##########
      #################################################################
      null_spx <- array_reshape(null_spx, dim = c(dims[2], dims[1]))  #
      true_Null <- true_null_spaxel_map(null_spx, 2, 1)               #
      true_Null <- array_reshape(true_Null, dim = c(dims[2], dims[1]))#
      #################################################################
      
      ################################################################
      ###### Creating a FITS file to hold a map of null spaxels ######
      ################################################################
      temp <- true_Null                                              #
      temp[which(temp == TRUE, arr.ind = TRUE)] <- 1                 #
      temp <- (temp - 1) * (-1)                                      #
                                                                     #
      writeFITSim(temp, file = paste0(out_path, '-true-Nulls.fits')) #
      rm(temp)                                                       #
      print('True Zero Map has been constructed.')                   #
      ################################################################

      #####################################
      # To hold positions of null spaxels #
      #####################################
      null_ind <- which(true_Null == TRUE)#
      rm(true_Null)                       #
      #####################################
      
      ###############################################################
      ################### SAMPLING OF TRUE NULLS ####################
      ############ To be used to avoid overflowing INLA #############
      ###############################################################
      grid_null_ind <- intersect(null_ind, grid_mask)               #
      samp_null_size <- round(length(grid_null_ind) / null_skip)    #
      samp_null_ind <- sample(grid_null_ind, samp_null_size)        #
                                                                    #
      out_null_ind <- setdiff(grid_null_ind, samp_null_ind)         #
      grid_mask <- grid_mask[-out_null_ind]                         #
      grid_size <- length(grid_mask)                                #
      rm(grid_null_ind, out_null_ind)                               #
      ###############################################################
      
      #################################################################
      ## Creating a FITS file to hold a map of sampled null spaxels ###
      #################################################################
      true_Null_Mask <- array(FALSE, dim = c(dims[2], dims[1]))       #
      true_Null_Mask[samp_null_ind] <- TRUE                           #
      temp <- true_Null_Mask                                          #
      temp[which(temp == TRUE)] <- 1                                  #
      temp[which(temp == FALSE)] <- 0                                 #
      temp <- (temp - 1) * (-1)                                       #
                                                                      #
      writeFITSim(temp, file = paste0(out_path, '-sample-Nulls.fits'))#
      rm(samp_null_ind, temp, true_Null_Mask)                         #
      print('True Zero Sample Map has been constructed.')             #
      #################################################################

      ###############################################################
      ############# Setting up data to input to encoder #############
      ###############################################################
      input_Net <- array_reshape(input_hold, dim = c(spatial_size,  #
                                                     feature_size)) #
      ###############################################################
      
      ###############################################################
      ######### Running test data through the encoder model #########
      ###############################################################
      encoded_u_input = encoder_mean %>% predict(input_Net)         #
      encoded_s_input = encoder_var %>% predict(input_Net)          #
      print("Encoding: Done.")                                      #
                                                                    #
      encoded_z_layer = layer_concatenate(list(encoded_u_input,     #
                                               encoded_s_input)) %>%#
        layer_lambda(sampling)                                      #
                                                                    #
      formatter = keras_model(dec_input, dec_input * 1)             #
      encoded_z = formatter %>% predict(encoded_z_layer, steps = 1) #
      print("Sampling: Done.")                                      #
      ###############################################################
      
      ########################################################################
      ############ Saving the spatial maps of the latent features ############
      ########################################################################
      encoded_output <- array(NaN, dim = c(spatial_size, latent_size))       #
                                                                             #
      encoded_map <- array_reshape(encoded_z, dim = c(dims[1], dims[2],      #
                                                      latent_size))          #
      writeFITSim(encoded_map, paste0(out_path, "_Latent_Input.fits"))       #
                                                                             #
      # for(g in 1:grid_size){                                                 #
      #   encoded_output[grid_mask[g]] <- median(encoded_z[cell_pix[g,]],      #
      #                                          na.rm = TRUE)                 #
      # }                                                                      #
      encoded_output[grid_mask, ] <- encoded_z[grid_mask, ]                  #
      encoded_samp <- array_reshape(encoded_output, dim = c(dims[1], dims[2],#
                                                            latent_size))    #
      writeFITSim(encoded_samp, paste0(out_path, "_Latent_Sample.fits"))     #
      rm(encoded_samp)                                                       #
      ########################################################################

      ################################################
      ######## Variables to be used with INLA ########
      ################################################
      nSize <- max(c(dims[1], dims[2]), na.rm = TRUE)#
      X <- rep(1:dims[1], each = dims[2])            #
      Y <- rep(1:dims[2], times = dims[1])           #
      zoom = 1                                       #
      p_range = c(2, 0.2)                            #
      p_sigma = c(2, 0.2)                            #
      inla.setOption(fmesher.timeout = 10)           #
      ################################################
      
      ####################################################################
      ######### Arrays to hold data for and from INLA procedures #########
      ####################################################################
      time_holder <- array(NaN, dim = latent_size)                       #
      encoded_mean <- array(NaN, dim = c(dims[1], dims[2], latent_size)) #
      encoded_res <- array(NaN, dim = c(dims[1], dims[2], latent_size))  #
      ####################################################################
      
      ##########################################################################
      ########################## MAIN INLA WORK CYCLE ##########################
      ##########################################################################
      print(paste0("Starting to reconstruct latent feature maps. Sampling ",   #
                   samp_str, " of useful data from case: ", opt_depth,         #
                   " thickness ", pov_str, "-on."))                            #
                                                                               #
      for(s in 1:latent_size){                                                 #
                                                                               #
        # DATA TRANSFORM                                                       #
        Z <- transf_by_method(encoded_output[, s], "None")                     #
                                                                               #
        fits_df <- data.frame(x = X, y = Y, z = Z)                             #
        new_exmag <- na.omit(fits_df[grid_mask, ])                             #
                                                                               #
        exCoord <- cbind(new_exmag$x, new_exmag$y)                             #
                                                                               #
        try_flag <- TRUE                                                       #
        seg_flag <- FALSE                                                      #
        man_flag <- FALSE                                                      #
        cutOff <- JUMP - 0.01                                                  #
        while(try_flag || seg_flag || man_flag){                               #
          try_flag <- FALSE                                                    #
          seg_flag <- FALSE                                                    #
          man_flag <- FALSE                                                    #
          # INLA Procedures                                                    #
          start_time <- Sys.time()                                             #
          meshError <- tryCatch(                                               #
            mesh <- try(inla.mesh.2d(exCoord, max.n = nSize, cutoff = cutOff)),#
            error = function(e)                                                #
              e)                                                               #
          end_time <- Sys.time()                                               #
          mesh_time <- end_time - start_time                                   #
          if(inherits(mesh, "try-error")) {                                    #
            print(paste0("INLA fmesher error, time exceeded. ",                #
                         "Increasing cutoff distance and trying again..."))    #
            try_flag <- TRUE                                                   #
            cutOff <- (cutOff + .01) * sqrt(2) - 0.01                          #
            next                                                               #
          }                                                                    #
          if(inherits(meshError, "error")) {                                   #
            print(paste0("INLA fmesher error, 'manifold' issue. ",             #
                         "Increasing cutoff distance and trying again..."))    #
            man_flag <- TRUE                                                   #
            cutOff <- (cutOff + .01) * sqrt(2) - 0.01                          #
            next                                                               #
          }                                                                    #
          print(paste0('INLA mesh.2d() took: ', mesh_time, ' seconds'))        #
                                                                               #
          start_time <- Sys.time()                                             #
          A <- inla.spde.make.A(mesh = mesh, loc = exCoord)                    #
          end_time <- Sys.time()                                               #
          make_time <- end_time - start_time                                   #
          print(paste0('INLA spde.make.A() took: ', make_time, ' seconds'))    #
                                                                               #
          start_time <- Sys.time()                                             #
          projection <- inla.mesh.projector(mesh, xlim = c(1, dims[1]),        #
                                            ylim = c(1, dims[2]),              #
                                            dim = zoom * c(dims[1] + 1,        #
                                                           dims[2] + 1))       #
          end_time <- Sys.time()                                               #
          proj_time <- end_time - start_time                                   #
          print(paste0('INLA mesh.projector() took: ', proj_time, ' seconds')) #
                                                                               #
          start_time <- Sys.time()                                             #
          spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 2,                  #
                                      prior.range = p_range,                   #
                                      prior.sigma = p_sigma)                   #
          end_time <- Sys.time()                                               #
          pcmatern_time <- end_time - start_time                               #
          print(paste0('INLA spde2.pcmatern() took: ', pcmatern_time,          #
                       ' seconds'))                                            #
                                                                               #
          start_time <- Sys.time()                                             #
          stk <- inla.stack(                                                   #
            data = list(par = new_exmag$z), A = list(A, 1),                    #
            effects = list(i = 1:spde$n.spde,                                  #
                           intercept = rep(1, length(new_exmag$x))),           #
            tag = 'est')                                                       #
          end_time <- Sys.time()                                               #
          stack_time <- end_time - start_time                                  #
          print(paste0('INLA stack() took: ', stack_time, ' seconds'))         #
                                                                               #
          start_time <- Sys.time()                                             #
          possibleError <- tryCatch(                                           #
            res <- inla(new_exmag$z ~ -1 + intercept + f(i, model = spde),     #
                        data = inla.stack.data(stk),                           #
                        control.predictor = list(A = inla.stack.A(stk),        #
                                                 compute = TRUE), verbose = F),#
            error = function(e)                                                #
              e)                                                               #
          end_time <- Sys.time()                                               #
          inla_time <- end_time - start_time                                   #
          if(inherits(possibleError, "error")) {                               #
            print(paste0("INLA reconstruction error, segmentation fault. ",    #
                         "Increasing cutoff distance and trying again..."))    #
            seg_flag <- TRUE                                                   #
            cutOff <- cutOff * sqrt(2)                                         #
            next                                                               #
          }                                                                    #
        }                                                                      #
        print(paste0('INLA inla() took: ', inla_time, ' seconds'))             #
                                                                               #
        start_time <- Sys.time()                                               #
                                                                               #
        # OUTPUT MEAN VALUE                                                    #
        output_mean <- inla.mesh.project(                                      #
          inla.mesh.projector(mesh, xlim = c(1, dims[1]), ylim = c(1, dims[2]),#
                              dim = zoom * c(dims[1], dims[2])),               #
          res$summary.random$i$mean) +                                         #
          t(matrix(as.numeric(res$summary.fixed$mean[1]),                      #
                   nrow = zoom * (dims[1]), ncol = zoom * (dims[2])))          #
                                                                               #
        # OUTPUT SD                                                            #
        output_sd <- inla.mesh.project(                                        #
          inla.mesh.projector(mesh, xlim = c(1, dims[1]), ylim = c(1, dims[2]),#
                              dim = c(dims[1], dims[2])),                      #
          res$summary.random$i$sd)                                             #
                                                                               #
        end_time <- Sys.time()                                                 #
        project_time <- end_time - start_time                                  #
        print(paste0('INLA projecting solution took: ', project_time,          #
                     ' seconds'))                                              #
                                                                               #
        total_inla_time <- project_time + inla_time + stack_time +             #
          pcmatern_time + proj_time + make_time + mesh_time                    #
        print(paste0('INLA in total took: ', total_inla_time, ' seconds'))     #
                                                                               #
        time_holder[s] <- total_inla_time                                      #
        # End of INLA Procedures                                               #
                                                                               #
        print(paste0("End of iteration ", s," out of ", latent_size, "."))     #
                                                                               #
        # encoded_mean[, , s] <- output_mean - minZ_3  + minZ                  #
        encoded_mean[ , , s] <- inv_transf_by_method(output_mean,              #
                                                     encoded_output[,s], "None")
        # encoded_sd[, , s] <- 10^output_sd                                    #
      }                                                                        #
                                                                               #
      print(paste0("Total Execution Time: ", sum(time_holder)/60, " minutes."))#
      ##########################################################################

      time_array[case_n, "Case"] <- paste0("optDepth=", opt_depth, "_pov=",
                                           pov_str, "_pp=1e", counts_exp, 
                                           "_samp=", samp_str)
      time_array[case_n, "Time (min)"] <- toString(sum(time_holder)/60)
      case_n <- case_n + 1

      ####################################################################
      ###### Extract SEDs from INLA reconstructions on latent space ######
      ############# and carry them back through decoder model ############
      ####################################################################
      encoded_mean_v <- array_reshape(encoded_mean,                      #
                                      dim = c(spatial_size, latent_size))#
      decoded_mean <- decoder %>% predict(encoded_mean_v)                #
      decoded_mean[null_ind, ] <- 0                                      #
      ####################################################################
      
      ######################################################################
      ##### Calculating map of residuals and reshaping arrays into maps ####
      ######################################################################
      # Percentual residuals plots won't handle 0's on encoded_map         #
      encoded_res <- abs((encoded_mean - encoded_map) / encoded_map) * 100 #
      encoded_res[which(is.infinite(encoded_res), arr.ind = TRUE)] <- NaN  #
      decoded_mean <- array_reshape(decoded_mean, dim = c(dims[1], dims[2],#
                                                          feature_size))   #
      ######################################################################
      
      #####################################################################
      ################# Saving All fits cubes of interest #################
      #####################################################################
      writeFITSim(encoded_mean, file = paste0(out_path,                   #
                                              '_Latent_INLA_MEAN.fits'))  #
      writeFITSim(encoded_res, file = paste0(out_path,                    #
                                             '_Latent_INLA_Dif.fits'))    #
      writeFITSim(decoded_mean, file = paste0(out_path,                   #
                                              '_Decoded_INLA_MEAN.fits')) #
      #####################################################################
    }
  }
}

write.table(time_array, file = paste0(results_folder, "emulation_times_exp", 
                                      counts_exp, ".csv"))