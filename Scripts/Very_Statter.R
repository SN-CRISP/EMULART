rm(list = ls())
options(digits = 2)
require(ggplot2)
require(gtools)
require(FITSio)
require(keras)
require(jjb)
source("/media/joaomfras/TOSHIBA EXT/EmulaRT/Scripts/SA_lib.R")


main_path <- "/media/joaomfras/TOSHIBA EXT/EmulaRT/"
input_path <- paste0(main_path, "Sources/FITS/")
results_path <- paste0(main_path, "Results/")
scale_path <- paste0(main_path, "Sources/Ref_Scale.dat")

opt_depth_arr <- c("0.05", "0.1", "1.0")
pov_str_arr <- c("edge", "face")
pov_num_arr <- c("90", "00")
counts <- c(4, 5, 6, 7)
counts_exp_ref <- 8
naturals <- c(2, 3, 5)#, 6, 7, 8

cnt_size <- length(counts)
nat_size <- length(naturals)
loop_size <- cnt_size * nat_size

for(opt_depth in opt_depth_arr){
  for(pov_str in pov_str_arr){
      
    pov_num <- pov_num_arr[which(pov_str_arr == pov_str, arr.ind = TRUE)]
    depth_path <- paste0("/t9-", opt_depth, "/", pov_str, "/")
    results_depth_fold <- paste0("reg-grid_t9-", opt_depth, "/", pov_str, "/")
    
    #######################################################################
    ## Change if you have a new ground truth / reference simulation cube ##
    #######################################################################
    ref_path <- paste0(input_path, depth_path, "sphShell_t9-", opt_depth, #
                       "_aniso_pp1e", counts_exp_ref, "_i", pov_num,      #
                       "_secondarydirect.fits")                           #
    #######################################################################
    
    ref_fits <- readFITS(ref_path)
    first_useful_feature <- 40
    dims <-  dim(ref_fits$imDat)
    features <- first_useful_feature:dims[3]
    spatial_size <- dims[1] * dims[2]
    feature_size <- length(features)
    
    ref <- ref_fits$imDat[, , features]
    ref_line <- array_reshape(ref, dim = c(spatial_size, feature_size))
    ref_tag <- "1e8 ref"
    intRef <- apply(ref, 3, sum, na.rm = TRUE)
    rm(ref_fits)
    
    print(paste0("Starting to process case: Optical Depth = ", opt_depth,
                 ", Point of View = ", pov_str, "-on..."))
    for(counts_exp in counts){
      
      src_tag <- paste0("1e", counts_exp, " src")
      emul_tag <- paste0("1e", counts_exp, " emul")
      spx_tags <- c(ref_tag, emul_tag, src_tag)
      
      for(nat in naturals){
        SKIP <- nat^2
        non_null_samp_rat <- 1 / SKIP
        samp_str <- paste0(round(non_null_samp_rat * 100), "p100")

        share_name <- paste0("sphShell_t9-", opt_depth, "_aniso_pp1e",
                             counts_exp, "_i", pov_num, "_secondarydirect")

        print(paste0("Loading data files for subcase: photon count = 10^", 
                     counts_exp, ", sampling = ", 
                     round(non_null_samp_rat * 100), "%..."))

        ########################################################################
        ############### Change according to out_path last folder ###############
        ########################################################################
        source_path <- paste0(input_path, depth_path, share_name, ".fits")     #
        case_path <- paste0(results_path, results_depth_fold, "1e", counts_exp,#
                            "-", samp_str, "/")                                #
        emul_path <- paste0(case_path, "Emulation/")                           #
        latent_map_path <- paste0(emul_path, share_name, "_Latent_Input.fits") #
        latent_inla_path <- paste0(emul_path, share_name,                      #
                                   "_Latent_INLA_MEAN.fits")                   #
        emul_path <- paste0(emul_path, share_name, "_Decoded_INLA_MEAN.fits")  #
                                                                               #
        target_sed_dir <- paste0(case_path, "SEDs/")                           #
        if(!dir.exists(target_sed_dir)){                                       #
          mkdir(target_sed_dir)                                                #
        }                                                                      #
                                                                               #
        target_res_dir <- paste0(case_path, "Residuals/")                      #
        if(!dir.exists(target_res_dir)){                                       #
          mkdir(target_res_dir)                                                #
        }                                                                      #
                                                                               #
        file_first_name <- paste0("t9-", opt_depth, "_i", pov_num, "_pp1e",    #
                                  counts_exp, "_", samp_str, "_")              #
        ########################################################################
        
        ############################################
        ########### Loading input files ############
        ############################################
        source_fits <- readFITS(source_path)       #
        lat_input_fits <- readFITS(latent_map_path)#
        lat_inla_fits <- readFITS(latent_inla_path)#
        emul_fits <- readFITS(emul_path)           #
        ############################################
        
        ########################################################################
        ######## Preparing input data to calculate performance metrics #########
        ########################################################################
        latent_size <- dim(lat_input_fits$imDat)[3]                            #
                                                                               #
        source <-source_fits$imDat[, , features]                               #
        emul <- emul_fits$imDat                                                #
        lat_input <- array_reshape(lat_input_fits$imDat, dim = c(spatial_size, #
                                                                 latent_size)) #
        lat_inla <- array_reshape(lat_inla_fits$imDat, dim = c(spatial_size,   #
                                                               latent_size))   #
                                                                               #
        intSource <- apply(source, 3, sum, na.rm = TRUE)                       #
        intEmul <- apply(emul, 3, sum, na.rm = TRUE)                           #
                                                                               #
        rm(emul_fits, lat_inla_fits, lat_input_fits)                           #
        ########################################################################

        ########################################################################
        ############ Load um scale of the file this is thought for #############
        ########################################################################
        myScale <- read.csv(scale_path, header = FALSE, sep = " ")[features, 1]#
        ########################################################################
        
        print("Calculating performance metrics...")
        ########################################################################
        ################## Calculation of performance metrics ##################
        ########################################################################
        res_emul <- abs(ref - emul) / ref * 100                                #
        res_source <- abs(ref - source) / ref * 100                            #
        intRes_emul <- abs(intRef - intEmul) / intRef * 100                    #
        intRes_source <- abs(intRef - intSource) / intRef * 100                #
        intRes_emul[which(is.infinite(intRes_emul), arr.ind = TRUE)] <- NA     #
        intRes_source[which(is.infinite(intRes_source), arr.ind = TRUE)] <- NA #
        res_emul[which(is.infinite(res_emul), arr.ind = TRUE)] <- NA           #
        res_source[which(is.infinite(res_source), arr.ind = TRUE)] <- NA       #
                                                                               #
        #Median - Median of Residuals per wavelength                           #
        #MAD - Mean Absolute Deviation of Residuals per wavelength             #
        #Int - Residuals of spatial sum per wavelength                         #
        res_emul_stats_pF<- array(NA, dim = c(feature_size, 4),                #
                                   dimnames = list(NULL, c("Wavelength (μm)",  #
                                                           "Median", "MAD",    #
                                                           "Int")))            #
        res_source_stats_pF<- array(NA, dim = c(feature_size, 4),              #
                                     dimnames = list(NULL, c("Wavelength (μm)",#
                                                             "Median", "MAD",  #
                                                             "Int")))          #
        res_emul_stats_pF[, "Wavelength (μm)"] <- myScale                      #
        res_source_stats_pF[, "Wavelength (μm)"] <- myScale                    #
        res_emul_stats_pF[, "Median"] <- apply(res_emul, 3, median, na.rm = T) #
        res_source_stats_pF[, "Median"] <- apply(res_source, 3, median,        #
                                                 na.rm = T)                    #
        res_emul_stats_pF[, "MAD"] <- apply(res_emul, 3, mad, na.rm = TRUE)    #
        res_source_stats_pF[, "MAD"] <- apply(res_source, 3, mad, na.rm = TRUE)#
        res_emul_stats_pF[, "Int"] <- intRes_emul                              #
        res_source_stats_pF[, "Int"] <- intRes_source                          #
                                                                               #
        #Median - Median of Residuals                                          #
        #MAD - Mean Absolute Deviation of Residuals                            #
        #TIR - Total Information Ratio                                         #
        #Null_SPX - Number of empty spaxels                                    #
        #Incomp_SPX - Number of incomplete spaxels                             #
        #Full_SPX - Number of complete spaxels                                 #
        perf_emul_stats <- array(NA, dim = 6,                                  #
                                 dimnames = list(c("Median", "MAD", "TIR",     #
                                                   "Null_SPX", "Incomp_SPX",   #
                                                   "Full_SPX")))               #
        perf_source_stats <- array(NA, dim = 6,                                #
                                   dimnames = list(c("Median", "MAD", "TIR",   #
                                                     "Null_SPX", "Incomp_SPX", #
                                                     "Full_SPX")))             #
        perf_emul_stats["Median"] <- median(res_emul, na.rm = TRUE)            #
        perf_source_stats["Median"] <- median(res_source, na.rm = TRUE)        #
        perf_emul_stats["MAD"] <- mad(res_emul, na.rm = TRUE)                  #
        perf_source_stats["MAD"] <- mad(res_source, na.rm = TRUE)              #
        perf_emul_stats["TIR"] <- info_size(emul, 0) / info_size(ref, 0)       #
        perf_source_stats["TIR"] <- info_size(source, 0) / info_size(ref, 0)   #
        perf_emul_stats["Null_SPX"] <- null_spx_size(emul, 1, 2, 0)            #
        perf_source_stats["Null_SPX"] <- null_spx_size(source, 1, 2, 0)        #
        perf_emul_stats["Incomp_SPX"] <- incomp_spx_size(emul, 1, 2, 0)        #
        perf_source_stats["Incomp_SPX"] <- incomp_spx_size(source, 1, 2, 0)    #
        perf_emul_stats["Full_SPX"] <- full_spx_size(emul, 1, 2, 0)            #
        perf_source_stats["Full_SPX"] <- full_spx_size(source, 1, 2, 0)        #
        ########################################################################
        
        print("Creating CSVs and FITS for performance statistics...")
        ########################################################################
        ########### Creating files to hold the normalized residuals ############
        ########################################################################
        write.table(res_emul_stats_pF,                                         #
                    file = paste0(target_res_dir, file_first_name,             #
                                  "Res-Emul-Stats-pWL.csv"))                   #
        write.table(perf_emul_stats,                                           #
                    file = paste0(target_res_dir, file_first_name,             #
                                  "Perf-Emul-Stats.csv"))                      #
        write.table(res_source_stats_pF,                                       #
                    file = paste0(target_res_dir, file_first_name,             #
                                  "Res-Source-Stats-pWL.csv"))                 #
        write.table(perf_source_stats,                                         #
                    file = paste0(target_res_dir, file_first_name,             #
                                  "Perf-Source-Stats.csv"))                    #
        ########################################################################
        
        ########################################################################
        ########### Creating a data cube of the normalized residuals ###########
        ########################################################################
        res_emul4Fits <- array_reshape(res_emul, dim = c(dims[1], dims[2],     #
                                                         feature_size))        #
        res_source4Fits <- array_reshape(res_source, dim = c(dims[1], dims[2], #
                                                             feature_size))    #
        writeFITSim(res_emul, file = paste0(target_res_dir, file_first_name,   #
                                            "Res-Emul.fits"))                  #
        writeFITSim(res_source, file = paste0(target_res_dir, file_first_name, #
                                              "Res-Source.fits"))              #
        ########################################################################
        
        ########################################################################
        ########## Reshaping Arrays to simplify plotting instructions ##########
        ########################################################################
        src_line <- array_reshape(source, dim = c(spatial_size, feature_size)) #
        emul_line <- array_reshape(emul, dim = c(spatial_size, feature_size))  #
        res_emul <- array_reshape(res_emul, dim = c(spatial_size, feature_size))
        res_source <- array_reshape(res_source, dim = c(spatial_size,          #
                                                        feature_size))         #
        ########################################################################
        
        print("Creating PDF for histogram of emulation residuals...")          #
        ########################################################################
        ######## Printing to file histograms of the performance metrics ########
        ########################################################################
        filename <- paste0(target_res_dir, file_first_name, 'Res_pF.pdf')      #
        pdf(filename, width = 22, height = 24)                                 #
        par(mfrow = c(4, 2))                                                   #
        par(cex.main = 1.5)                                                    #
                                                                               #
        for(f in 1:feature_size){                                              #
                                                                               #
          temp <- log10(res_emul[, f])                                         #
          temp[which(is.infinite(temp))] <- NaN                                #
          min_t <- min(temp, na.rm = TRUE)                                     #
          if(min_t >= 0){                                                      #
            sig_t <- -1                                                        #
          }else{                                                               #
            sig_t <- 1                                                         #
          }                                                                    #
          max_t <- max(temp, na.rm = TRUE)                                     #
                                                                               #
          hist(na.omit(temp),                                                  #
               breaks = seq(sig_t * 1.25 * min_t, 1.25 * max_t, 0.25),         #
               col = "red", freq = TRUE, xlim = c(-3, 10),                     #
               xaxp  = c(-3, 10, 26), cex.axis = 2, cex.lab = 2, cex.sub = 2,  #
               xlab = "log10(Residuals[%])", ylab = "Frequency",               #
               main = paste0("Emulation Residuals for Feature ", f))           #
          lines(rep(log10(res_emul_stats_pF[f, "Median"]), 2), c(0, 90000),    #
                type = "h", col = "blue")                                      #
        }                                                                      #
        dev.off()                                                              #
                                                                               #
        filename <- paste0(target_res_dir, file_first_name, 'Res_int.pdf')     #
        pdf(filename, width = 24, height = 16)                                 #
        par(mfrow = c(1, 1))                                                   #
        par(cex.main = 1.5)                                                    #
                                                                               #
        limRes <- 10^ceiling(log10(max(intRes_emul, na.rm = TRUE)))            #
        byRes <- limRes / 40                                                   #
                                                                               #
        hist(na.omit(intRes_emul), breaks = seq(0, limRes, byRes), col = "red",#
             freq = TRUE, xlim = c(0, 200), xaxp  = c(0, 200, 80),             #
             cex.axis = 2, cex.lab = 2, xlab = "Residuals[%]",                 #
             ylab = "Frequency", main = "Emulation Residuals")                 #
        lines(rep(perf_emul_stats["Median"], 2), c(0, 100), type = "h",        #
              col = "blue")                                                    #
                                                                               #
        dev.off()                                                              #
        ########################################################################
        
        ########################################################################
        ########## Printing to file some of the reconstructed results ##########
        ########################################################################
        print_index_sample <- c(44850, 44910, 44970)                           #
                                                                               #
        ax_x <- myScale                                                        #
        ax_y <- 10^seq(-50, 0, 1)                                              #
        ax_y_res <- seq(0, 100, 20)                                            #
                                                                               #
        print("Creating PDF for individual spaxels...")                        #
        filename <- paste0(target_sed_dir, file_first_name, '_Spx_SEDs.pdf')   #
        pdf(filename, width = 10, height = 6)                                  #
        par(mar=c(5.1, 5.6, 5.1, 0))                                           #
                                                                               #
        for (i in print_index_sample){                                         #
                                                                               #
          matrix_coord <- c((i - 1) %/% dims[2] + 1, (i - 1) %% dims[2] + 1)   #
          latent_feature_string <- toString(lat_input[i, ], sep = ",")         #
          latent_inla_string <- toString(lat_inla[i, ], sep = ",")             #
                                                                               #
          spx_bind <- log10(cbind(ref_line[i, ], emul_line[i, ], src_line[i, ]))
          spx_bind[which(is.infinite(spx_bind))] <- NaN                        #
          min_dec <- min(spx_bind, na.rm = TRUE)                               #
          max_dec <- max(spx_bind, na.rm = TRUE)                               #
                                                                               #
          matplot(log10(myScale), spx_bind, type = "b", pch = 1:3, cex.lab = 2,#
                  xlab = expression(paste("Wavelength (", mu, "m)")),          #
                  ylab = "Flux Density [W/m²]", ylim = c(min_dec, max_dec),    #
                  col = c("black", "goldenrod", "hotpink"), axes = FALSE)      #
          axis(1, at = log10(ax_x), labels = sprintf(ax_x, fmt = '%#.1f'),     #
               cex.axis = 1.5)                                                 #
          axis(2, at = log10(ax_y), labels = ax_y, cex.axis = 1.5)             #
          legend(2.2, max_dec, legend = spx_tags, cex = 1.5,                   #
                 col = c("black", "goldenrod", "hotpink"), lty = 1:3, pch = 1:3)
                                                                               #
          mtext(paste0("Spaxel (row, col): (", matrix_coord[1], ", ",          #
                       matrix_coord[2], ")"), side = 3, col = "black",         #
                line = -1, at = 2.2)                                           #
                                                                               #
          res_bind <- cbind(res_emul[i, ], res_source[i, ])                    #
          res_bind[which(is.infinite(res_bind))] <- NaN                        #
          min_res <- min(res_bind, na.rm = TRUE)                               #
          max_res <- max(res_bind, na.rm = TRUE)                               #
                                                                               #
          matplot(log10(myScale), res_bind, type = "b", pch = 2:3, cex.lab = 2,#
                  xlab = expression(paste("Wavelength [", mu, "m]")),          #
                  ylab="Normalized Residuals [%]",                             #
                  col = c("goldenrod", "hotpink"), axes = FALSE)               #
          axis(1, at = log10(ax_x), labels = sprintf(ax_x, fmt = '%#.1f'),     #
               cex.axis = 1.5)                                                 #
          axis(2, at = ax_y_res, labels = ax_y_res, cex.axis = 1.5)            #
          legend(2.2, max_res, legend = c(emul_tag, src_tag), cex = 1.5,       #
                 col = c("goldenrod", "hotpink"), lty = 2:3, pch = 2:3)        #
          title("Normalized Residuals of Spaxel Emulation", line = 2.5)        #
        }                                                                      #
        dev.off()                                                              #
                                                                               #
        print("Creating PDF for Integrated SEDs...")                           #
        filename <- paste0(target_sed_dir, file_first_name, '_Int_SEDs.pdf')   #
        pdf(filename, width = 10, height = 10)                                 #
        par(mar=c(5.1, 5.6, 5.1, 0))                                           #
                                                                               #
        int_bind <- log10(cbind(intRef, intEmul, intSource))                   #
        int_bind[which(is.infinite(int_bind))] <- NaN                          #
        min_int <- min(int_bind, na.rm = TRUE)                                 #
        max_int <- max(int_bind, na.rm = TRUE)                                 #
                                                                               #
        matplot (log10(myScale), int_bind, type = "b", pch = 1:3, cex.lab = 2, #
                 col = c("black", "goldenrod", "hotpink"),                     #
                 xlab = expression(paste("Wavelength [", mu, "m]")),           #
                 ylab = "Flux Density [W/m²]", axes = F)                       #
        axis(1, at = log10(ax_x), labels = sprintf(ax_x, fmt = '%#.1f'),       #
             cex.axis = 1.5)                                                   #
        axis(2, at = log10(ax_y), labels = ax_y, cex.axis = 1.5)               #
        legend(2.2, max_int, legend = spx_tags, cex = 1.5,                     #
               col = c("black", "goldenrod", "hotpink"), lty = 1:3, pch = 1:3) #
        title("Truth and Emulatation Spatial Integration SEDs", line = 2.5)    #
                                                                               #
        int_res_bind <- cbind(intRes_emul, intRes_source)                      #
        int_res_bind[which(is.infinite(int_res_bind))] <- NaN                  #
        min_int_res <- min(int_res_bind, na.rm = TRUE)                         #
        max_int_res <- max(int_res_bind, na.rm = TRUE)                         #
                                                                               #
        matplot(log10(myScale), int_res_bind, type = "b", pch = 2:3,           #
                col = c("goldenrod", "hotpink"), cex.lab = 2,                  #
                xlab = expression(paste("Wavelength [", mu, "m]")),            #
                ylab = "Normalized Residuals [%]", axes = FALSE)               #
        axis(1, at = log10(ax_x), labels = sprintf(ax_x, fmt = '%#.1f'),       #
             cex.axis = 1.5)                                                   #
        axis(2, at = ax_y_res, labels = ax_y_res, cex.axis = 1.5)              #
        legend(2.2, max_int_res, legend = c(emul_tag, src_tag), cex = 1.5,     #
               col = c("goldenrod", "hotpink"), lty = 2:3, pch = 2:3)          #
        title("Normalized Residuals of Spatial Integration of Emulation",      #
              line = 2.5)                                                      #
        dev.off()                                                              #
        ########################################################################
      }
    }
  }
}