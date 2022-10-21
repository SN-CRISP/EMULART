rm(list = ls())
options(digits = 2)
require(ggplot2)
require(gtools)
require(FITSio)
require(jjb)
###########################################################################
# Key folder structures, change according to your need. Bare in mind that #
# the present script assumes a specific folder tree structure, so you may #
# have to change paths further down the script if you do not abide by it. #
###########################################################################
main_path <- "/media/joaomfras/TOSHIBA EXT/EmulaRT/"                      #
reslt_path <- paste0(main_path, "Results/reg-grid_t9-")                   #
source_path <- paste0(main_path, "Sources/FITS/t9-")                      #
scale_path <- paste0(main_path, "/Sources/Ref_Scale.dat")                 #
###########################################################################

#############################################################
## Key parameters regarding emulations change according to ##
######################## your needs. ########################
#############################################################
opt_depth_arr <- c("0.05", "0.1", "1.0")                    #
pov_str_arr <- c("edge", "face")                            #
pov_num_arr <- c("90", "00")                                #
counts <- c(4, 5, 6, 7)                                     #
counts_ref <- 8                                             #
naturals <- c(2, 3, 5)                                      #
                                                            #
cnt_size <- length(counts)                                  #
nat_size <- length(naturals)                                #
loop_size <- cnt_size * nat_size                            #
first_feature <- 40                                         #
feature_size <- 64                                          #
features <- first_feature:(first_feature + feature_size - 1)#
#############################################################

###############################################################
############ Setting up labels for arrays that will ###########
############## hold data and performance metrics ##############
###############################################################
cnt_prcnt_labels <- NULL                                      #
cnt_prcnt_tags <- NULL                                        #
cnt_labels <- NULL                                            #
for(counts_exp in counts){                                    #
                                                              #
  new_cnt_label <- paste0("1e", counts_exp)                   #
  cnt_labels <- c(cnt_labels, new_cnt_label)                  #
                                                              #
  for(nat in naturals){                                       #
                                                              #
    SKIP <- nat^2                                             #
    non_null_samp_rat <- 1 / SKIP                             #
    samp_str <- paste0(round(non_null_samp_rat * 100), "p100")#
    samp_tag <- paste0(round(non_null_samp_rat * 100), "%")   #
                                                              #
    new_comb_label <- paste0(new_cnt_label, "_", samp_str)    #
    new_comb_tag <- paste0(new_cnt_label, "  ", samp_tag)     #
    cnt_prcnt_labels <- c(cnt_prcnt_labels, new_comb_label)   #
    cnt_prcnt_tags <- c(cnt_prcnt_tags, new_comb_tag)         #
  }                                                           #
}                                                             #
###############################################################

########################################################################
############ Load um scale of the file this is thought for #############
########################################################################
myScale <- read.csv(scale_path, header = FALSE, sep = " ")[features, 1]#
########################################################################

################################################################################
############ Loop over all parameters used values for the emulations ###########
### Calculate metrics for both the emulations and sources and generate files ###
################################################################################
#These two outside loops iterate over the general source case                  #
#Picking the optical depth and the point-of-view of the observer               #
################################                                               #
for(opt_depth in opt_depth_arr){                                               #
  for(pov_str in pov_str_arr){                                                 #
                                                                               #
    #########################################################################  #
    ########## Defining common values to be frequently used in the ##########  #
    ######################## two loops nested bellow ########################  #
    #########################################################################  #
    pov_num <- pov_num_arr[which(pov_str_arr == pov_str,                    #  #
                                 arr.ind = TRUE)]                           #  #
                                                                            #  #
    ref <- readFITS(paste0(source_path, opt_depth, "/", pov_str,            #  #
                           "/sphShell_t9-", opt_depth, "_aniso_pp",         #  #
                           "1e", counts_ref, "_i", pov_num,                 #  #
                           "_secondarydirect.fits"))$imDat[, , features]    #  #
    dims <- dim(ref)                                                        #  #
    ref_SED <- apply(ref, 3, sum, na.rm = T)                                #  #
                                                                            #  #
    spatial_size <- dims[1] * dims[2]                                       #  #
                                                                            #  #
    reslt_com <- paste0(reslt_path, opt_depth, "/", pov_str, "/")           #  #
    emuls <- array(NA, dim = c(dims[1], dims[2], feature_size, loop_size),  #  #
                   dimnames = list(NULL, NULL, NULL, cnt_prcnt_labels))     #  #
    emul_pWl_resd <- array(NA, dim = c(feature_size, loop_size),            #  #
                           dimnames = list(NULL, cnt_prcnt_labels))         #  #
    emul_int_resd <- array(NA, dim = c(feature_size, loop_size),            #  #
                           dimnames = list(NULL, cnt_prcnt_labels))         #  #
    #########################################################################  #
                                                                               #
    print(paste0("Starting to process case: Optical Depth = ", opt_depth,      #
                 ", Point of View = ", pov_str, "-on..."))                     #
    # These two loops iterate over the specific source case and the            #
    # emulation that results from using it as input.                           #
    # Picking the amount of photons simulated to produce the source            #
    # and the amount of source data sampled along a binned grid                #
    ##########################                                                 #
    for(counts_exp in counts){                                                 #
                                                                               #
      cnt_ind <- paste0("1e", counts_exp)                                      #
                                                                               #
      print("Loading files...")                                                #
      for(nat in naturals){                                                    #
                                                                               #
        SKIP <- nat^2                                                          #
        non_null_samp_rat <- 1 / SKIP                                          #
        samp_str <- paste0(round(non_null_samp_rat * 100), "p100")             #
        samp_tag <- paste0(round(non_null_samp_rat * 100), "%")                #
        comb_ind <- paste0(cnt_ind, "_", samp_str)                             #
                                                                               #
        print(paste0("Loading data files for subcase: photon count = 10^",     #
                     counts_exp, ", sampling = ", samp_tag, "..."))            #
                                                                               #
        fit_name <- paste0("sphShell_t9-", opt_depth, "_aniso_pp", cnt_ind,    #
                           "_i", pov_num, "_secondarydirect_")                 #
        reslt_files <- paste0(reslt_com, cnt_ind, "-", samp_str, "/")          #
        oth_name <- paste0("t9-", opt_depth, "_i", pov_num, "_pp", cnt_ind,    #
                           "_", samp_str, "_")                                 #
        resd_fold <- paste0(reslt_files, "Residuals/")                         #
        emul_fold <- paste0(reslt_files, "Emulation/")                         #
        emul_resd_path <- paste0(resd_fold, oth_name, "Res-Emul-Stats-pWL.csv")#
        emul_path <- paste0(emul_fold, fit_name, "Decoded_INLA_MEAN.fits")     #
                                                                               #
        temp_emul <- read.csv(emul_resd_path, header = TRUE, sep = " ")        #
        emul_pWl_resd[, comb_ind] <- temp_emul$Median                          #
        emul_int_resd[, comb_ind] <- temp_emul$Int                             #
        rm(temp_emul)                                                          #
                                                                               #
        emuls[, , , comb_ind] <- readFITS(emul_path)$imDat                     #
        rm(emul_path)                                                          #
      }                                                                        #
    }                                                                          #
                                                                               #
    ############################################################               #
    ###### Creating Integrated SEDs by summing all pixels ######               #
    #################### at each wavelength ####################               #
    ############################################################               #
    emuls_SED <- apply(emuls, 3:4, sum, na.rm = T)             #               #
    ############################################################               #
                                                                               #
    ###############################                                            #
    ###### Axes information #######                                            #
    ###############################                                            #
    ax_x <- myScale               #                                            #
    ax_y <- 10^seq(-50, -2, 1)    #                                            #
    ax_y_res <- 10^seq(-6, 15, .5)#                                            #
    ###############################                                            #
                                                                               #
    print("Creating PDFs for integrated SEDs and respective residuals...")     #
    ######################################################################     #
    # Creating PDFs with SEDs of emulations performed with same sampling #     #
    ######### pattern/percentage but different input simulations #########     #
    ######################################################################     #
    for(nat in naturals){                                                #     #
                                                                         #     #
      SKIP <- nat^2                                                      #     #
      non_null_samp_rat <- 1 / SKIP                                      #     #
      samp_str <- paste0(round(non_null_samp_rat * 100), "p100")         #     #
      samp_tag <- paste0(round(non_null_samp_rat * 100), "%")            #     #
      check_inds <- NULL                                                 #     #
      check_tags <- NULL                                                 #     #
                                                                         #     #
      for(counts_exp in counts){                                         #     #
        cnt_ind <- paste0("1e", counts_exp)                              #     #
        check_inds <- c(check_inds, paste0(cnt_ind, "_", samp_str))      #     #
        check_tags <- c(check_tags, paste0(cnt_ind, "  ", samp_tag))     #     #
      }                                                                  #     #
                                                                         #     #
      sed_bind <- NULL                                                   #     #
      sed_bind <- cbind(log10(ref_SED))                                  #     #
      int_resd_bind <- NULL                                              #     #
      for(c in check_inds){                                              #     #
        sed_bind <- cbind(sed_bind, log10(emuls_SED[, c]))               #     #
        int_resd_bind <- cbind(int_resd_bind, log10(emul_int_resd[, c])) #     #
      }                                                                  #     #
                                                                         #     #
      pdf(paste0(reslt_com, samp_str, "-int_sed_resd.pdf"), width = 10,  #     #
          height = 5)                                                    #     #
      par(mar=c(5.1, 5.6, 5.1, 0))                                       #     #
      min_int_resd <- min(int_resd_bind, na.rm = TRUE)                   #     #
      max_int_resd <- max(int_resd_bind, na.rm = TRUE)                   #     #
      matplot(log10(myScale), int_resd_bind, type = "b",                 #     #
              pch = 2:(1 + cnt_size), col = 2:(1 + cnt_size),            #     #
              cex.lab = 2, xlab = expression(paste("Wavelength (",       #     #
                                                   mu, "m)")),           #     #
              ylab = "Normalized Residuals [%]",                         #     #
              ylim = c(min_int_resd, max_int_resd), axes = FALSE)        #     #
      axis(1, at = log10(ax_x), labels = sprintf(ax_x, fmt = '%#.1f'),   #     #
           cex.axis = 1.5)                                               #     #
      axis(2, at = log10(ax_y_res), labels = sprintf(ax_y_res, fmt = '%#.1g'), # 
           cex.axis = 1.5)                                               #     #
      legend(2.5, max_int_resd, legend = check_tags,                     #     #
             pch = 2:(1 + cnt_size), col = 2:(1 + cnt_size), cex = 1.5)  #     #
      title("Normalized Residuals of Integrated SEDs", line = 2.5)       #     #
      dev.off()                                                          #     #
                                                                         #     #
      pdf(paste0(reslt_com, samp_str, "-int_sed.pdf"), width = 10,       #     #
          height = 10)                                                   #     #
      par(mar=c(5.1, 5.6, 5.1, 0))                                       #     #
      min_sed <- min(sed_bind, na.rm = TRUE)                             #     #
      max_sed <- max(sed_bind, na.rm = TRUE)                             #     #
      matplot(log10(myScale), sed_bind, type = "b",                      #     #
              pch = 1:(1 + cnt_size), col = 1:(1 + cnt_size),            #     #
              cex.lab = 2, xlab = expression(paste("Wavelength (",       #     #
                                                   mu, "m)")),           #     #
              ylab = "Flux Density [W/m²]",                              #     #
              ylim = c(min_sed, max_sed), axes = FALSE)                  #     #
      axis(1, at = log10(ax_x), labels = sprintf(ax_x, fmt = '%#.1f'),   #     #
           cex.axis = 1.5)                                               #     #
      axis(2, at = log10(ax_y), labels = ax_y, cex.axis = 1.5)           #     #
      legend(2.5, max_sed, legend = c("1e8 (Ref)", check_tags),          #     #
             pch = 1:(1 + cnt_size), col = 1:(1 + cnt_size), cex = 1.5)  #     #
      title("Integrated SEDs", line = 2.5)                               #     #
      dev.off()                                                          #     #
    }                                                                    #     #
    ######################################################################     #
                                                                               #
    print("Creating PDFs of SED for wavelength median of spatial residuals...")#
    ########################################################################## #
    ####### Creating PDFs with plots comparing SEDs and SEDs residuals ####### #
    ########################################################################## #
    pWl_resd_bind <- NULL                                                    # #
    for(l in cnt_prcnt_labels){                                              # #
      pWl_resd_bind <- cbind(pWl_resd_bind, log10(emul_pWl_resd[, l]))       # #
    }                                                                        # #
                                                                             # #
    pdf(paste0(reslt_com, "med_resd_pWl.pdf"), width = 10, height = 5)       # #
    par(mar=c(5.1, 5.6, 5.1, 0))                                             # #
    min_pWl_resd <- min(pWl_resd_bind, na.rm = TRUE)                         # #
    max_pWl_resd <- max(pWl_resd_bind, na.rm = TRUE)                         # #
    matplot(log10(myScale),  pWl_resd_bind, type = "b", pch = rep(2:4, 4),   # #
            col = c(rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3)),             # #
            cex.lab = 2, xlab = expression(paste("Wavelength (", mu, "m)")), # #
            ylab = "Normalized Residuals [%]",                               # #
            ylim = c(min_pWl_resd, 2.5), axes = FALSE)                       # #
    axis(1, at = log10(ax_x), labels = sprintf(ax_x, fmt = '%#.1f'),         # #
         cex.axis = 1.5)                                                     # #
    axis(2, at = log10(ax_y_res), labels = sprintf(ax_y_res, fmt = '%#.1f'), # #
         cex.axis = 1.5)                                                     # #
    legend(2.5, 9, legend = cnt_prcnt_tags, pch = rep(2:4, 4),               # #
           col = c(rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3)), cex = 1.5)   # #
    title("Median of Normalized Residuals per Wavelength", line = 2.5)       # #
    dev.off()                                                                # #
    ########################################################################## #
    rm(emuls, emul_int_resd, emuls_SED, emul_pWl_resd, ref_SED)                #
  }                                                                            #
}                                                                              #
################################################################################