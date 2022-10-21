############################################
# Data Analysis and Transformation Library #
############################################


###################################################################
## Creates a histogram which bins can be sampled the same amount ##
## Inputs:                                                       ##
##        data -> array of data to be binned;                    ##
##        samp_rat -> ratio of the data that is to be sampled;   ##
##        thrsh_type -> expected to be either "max", "avg" or    ##
##                      "null" will determine how to calculate   ##
##                      the count threshold for bin merger;      ##
##        null_n -> must be positive if 'thrsh_type' is 'null',  ##
##                  it will be set as the count threshold for    ##
##                  bin mergers.                                 ##
## Output:                                                       ##
##        histogram binned in such a way that all bins can be    ##
##        sampled the same amount.                               ##
###################################################################
create_nice_hist <- function(data, samp_rat, thrsh_type = "avg"){ #
  if(!is.vector(data) && !is.array(data)){                        #
    print(paste0("ERROR: 'data' is expected to either be a vecto",#
                 "r or an array. Returning NULL."))               #
    return(NULL)                                                  #
  }                                                               #
  if(length(data) == 0){                                          #
    print("ERROR: 'data' has length 0. Returning NULL.")          #
    return(NULL)                                                  #
  }                                                               #
  if(!is.numeric(data)){                                          #
    print(paste0("ERROR: 'data' is expected to be numeric. ",     #
                 "Returning NULL."))                              #
    return(NULL)                                                  #
  }                                                               #
  if(length(samp_rat) != 1){                                      #
    print(paste0("ERROR: 'samp_rat' should have length 1, and ",  #
                 "has length ", length(samp_rat), "instead. ",    #
                 "Returning NULL."))                              #
    return(NULL)                                                  #
  }                                                               #
  if(!is.numeric(samp_rat)){                                      #
    print(paste0("ERROR: 'samp_rat' is expected to be numeric. ", #
                 "Returning NULL."))                              #
    return(NULL)                                                  #
  }                                                               #
  if(samp_rat <= 0 || samp_rat > 1){                              #
    print(paste0("ERROR: 'samp_rat' is ", samp_rat, ", and not ", #
                 "between 0 and 1 as is expected. ",              #
                 "Returning NULL."))                              #
    return(NULL)                                                  #
  }                                                               #
  if(thrsh_type != "avg" && thrsh_type != "max"){                 #
    print(paste0("ERROR: 'thrsh_type' has an unexpected value: ", #
                 thrsh_type, ", it should be 'avg' or 'max'. ",   #
                 "Returning NULL."))                              #
    return(NULL)                                                  #
  }                                                               #
                                                                  #
  #Binning data                                                   #
  hist_info <- hist(data)                                         #
  bin_breaks <- hist_info$breaks                                  #
  bin_counts <- hist_info$counts                                  #
                                                                  #
  #Defining some limit quantities                                 #
  if(thrsh_type == 'max'){                                        #
    max_cts <- max(bin_counts)                                    #
    max_cts_brk <- max(which(bin_counts == max_cts))              #
    cts_threshold <- ceiling(samp_rat * max_cts)                  #
  }                                                               #
  if(thrsh_type == 'avg'){                                        #
    sum_cts <- sum(bin_counts)                                    #
    cts_threshold <- ceiling(samp_rat * sum_cts /                 #
                               length(bin_counts))                #
  }                                                               #
                                                                  #
  #Finding bins with less counts than samp_rat the bin max counts #
  brk2rm <- which(bin_counts < cts_threshold)                     #
                                                                  #
  #Creating temporary bins to be merged until meet cts_threshold  #
  temp_breaks <- bin_breaks[-length(bin_breaks)]                  #
  temp_counts <- bin_counts                                       #
  temp_brk2rm <- brk2rm                                           #
                                                                  #
  #Merging bins until they meet cts_threshold                     #
  while(length(temp_brk2rm) > 0){                                 #
                                                                  #
    bin <- temp_brk2rm[1]                                         #
    #Neighbors of bin to be merged                                #
    ln <- bin - 1                                                 #
    rn <- bin + 1                                                 #
    border <- FALSE                                               #
                                                                  #
    if(ln < 1){                                                   #
      temp_counts[bin] <- temp_counts[rn] + temp_counts[bin]      #
      temp_breaks <- temp_breaks[-rn]                             #
      temp_counts <- temp_counts[-rn]                             #
      border <- TRUE                                              #
    }                                                             #
    if(rn > length(temp_breaks)){                                 #
      temp_counts[ln] <- temp_counts[ln] + temp_counts[bin]       #
      temp_breaks <- temp_breaks[-bin]                            #
      temp_counts <- temp_counts[-bin]                            #
      border <- TRUE                                              #
    }                                                             #
    if(!border){                                                  #
      if(temp_counts[ln] <= temp_counts[rn]){                     #
        temp_counts[ln] <- temp_counts[ln] + temp_counts[bin]     #
        temp_breaks <- temp_breaks[-bin]                          #
        temp_counts <- temp_counts[-bin]                          #
      }else{                                                      #
        temp_counts[rn] <- temp_counts[rn] + temp_counts[bin]     #
        temp_breaks <- temp_breaks[-rn]                           #
        temp_counts <- temp_counts[-rn]                           #
      }                                                           #
    }                                                             #
    cts_threshold <- ceiling(samp_rat * sum_cts /                 #
                               length(temp_counts))               #
    temp_brk2rm <- which(temp_counts < cts_threshold)             #
  }                                                               #
                                                                  #
  mrg_bin_breaks <- append(temp_breaks,                           #
                           bin_breaks[length(bin_breaks)])        #
  rm(bin, border, ln, rn, temp_brk2rm, temp_breaks, temp_counts)  #
                                                                  #
  return(hist(data, breaks = mrg_bin_breaks, freq = TRUE))        #
}                                                                 #
###################################################################

#/////////////////////////////////////////////////////////////////#

###################################################################
#### Determines how many data points to sample from each bin. #####
############# That value being the same for every bin. ############
## Inputs:                                                       ##
##        counts -> vector containing each bins counts;          ##
##        samp_rat -> ratio of total counts to be sampled.       ##
## Output:                                                       ##
##        integer amount of counts to be sampled from each bin.  ##
###################################################################
equal_sampling_counts <- function(counts, samp_rat){              #
  if(!is.vector(counts) && !is.array(counts)){                    #
    print(paste0("ERROR: 'counts' is expected to either be ",     #
                 "a vector or an array. Returning NULL."))        #
    return(NULL)                                                  #
  }                                                               #
  if(!is.numeric(counts)){                                        #
    print(paste0("ERROR: 'counts' is expected to be numeric.",    #
                 " Returning NULL."))                             #
    return(NULL)                                                  #
  }                                                               #
  if(length(counts) == 0){                                        #
    print(paste0("ERROR: 'counts' has length 0. ",                #
                 "Returning NULL."))                              #
    return(NULL)                                                  #
  }                                                               #
  if(length(samp_rat) != 1){                                      #
    print(paste0("ERROR: 'samp_rat' should have length 1, and ",  #
                 "has length ", length(samp_rat), "instead. ",    #
                 "Returning NULL."))                              #
    return(NULL)                                                  #
  }                                                               #
  if(!is.numeric(samp_rat)){                                      #
    print(paste0("ERROR: 'samp_rat' is expected to be numeric. ", #
                 "Returning NULL."))                              #
    return(NULL)                                                  #
  }                                                               #
  if(samp_rat <= 0 || samp_rat > 1){                              #
    print(paste0("ERROR: 'samp_rat' is ", samp_rat, ", and not ", #
                 "between 0 and 1 as is expected. ",              #
                 "Returning NULL."))                              #
    return(NULL)                                                  #
  }                                                               #
  return(rep(round(sum(counts) * samp_rat / length(counts)),      #
             length(counts)))                                     #
}                                                                 #
###################################################################

#/////////////////////////////////////////////////////////////////#

###################################################################
#### Determines how many data points to sample from each bin. #####
#### That value being proportional to the counts of each bin. #####
## Inputs:                                                       ##
##        counts -> vector containing each bins counts;          ##
##        samp_rat -> ratio of total counts to be sampled.       ##
## Output:                                                       ##
##        integer vector of counts to be sampled from each bin.  ##
###################################################################
dist_sampling_counts <- function(counts, samp_rat){               #
  if(!is.vector(counts) && !is.array(counts)){                    #
    print(paste0("ERROR: 'counts' is expected to either be ",     #
                 "a vector or an array. Returning NULL."))        #
    return(NULL)                                                  #
  }                                                               #
  if(!is.numeric(counts)){                                        #
    print(paste0("ERROR: 'counts' is expected to be numeric.",    #
                 " Returning NULL."))                             #
    return(NULL)                                                  #
  }                                                               #
  if(length(counts) == 0){                                        #
    print(paste0("ERROR: 'counts' has length 0. ",                #
                 "Returning NULL."))                              #
    return(NULL)                                                  #
  }                                                               #
  if(length(samp_rat) != 1){                                      #
    print(paste0("ERROR: 'samp_rat' should have length 1, and ",  #
          "has length ", length(samp_rat), "instead. ",           #
          "Returning NULL."))                                     #
    return(NULL)                                                  #
  }                                                               #
  if(!is.numeric(samp_rat)){                                      #
    print(paste0("ERROR: 'samp_rat' is expected to be numeric. ", #
                 "Returning NULL."))                              #
    return(NULL)                                                  #
  }                                                               #
  if(samp_rat <= 0 || samp_rat > 1){                              #
    print(paste0("ERROR: 'samp_rat' is ", samp_rat, ", and not ", #
                 "between 0 and 1 as is expected. ",              #
                 "Returning NULL."))                              #
    return(NULL)                                                  #
  }                                                               #
  return(round(counts * samp_rat))                                #
}                                                                 #
###################################################################

#/////////////////////////////////////////////////////////////////#

########################################################################
####### Determines how many data points to sample from each bin. #######
###### That value being the same for every bin that has a minimum ######
#### number of counts. Every bin with a lower count amount is fully ####
##### sampled. The remaining samples (so as to achieve 'samp_rat') #####
######## are distributed equally among the remaining bins. #############
## Inputs:                                                            ##
##        counts -> vector containing each bins counts;               ##
##        samp_rat -> ratio of total counts to be sampled.            ##
## Output:                                                            ##
##        integer amount of counts to be sampled from each bin.       ##
########################################################################
base_equal_sampling_counts <- function(counts, samp_rat){              #
  if(!is.vector(counts) && !is.array(counts)){                         #
    print(paste0("ERROR: 'counts' is expected to either be ",          #
                 "a vector or an array. Returning NULL."))             #
    return(NULL)                                                       #
  }                                                                    #
  if(!is.numeric(counts)){                                             #
    print(paste0("ERROR: 'counts' is expected to be numeric.",         #
                 " Returning NULL."))                                  #
    return(NULL)                                                       #
  }                                                                    #
  if(length(counts) == 0){                                             #
    print(paste0("ERROR: 'counts' has length 0. Returning NULL."))     #
    return(NULL)                                                       #
  }                                                                    #
  if(length(samp_rat) != 1){                                           #
    print(paste0("ERROR: 'samp_rat' should have length 1, and ",       #
                 "has length ", length(samp_rat), "instead. ",         #
                 "Returning NULL."))                                   #
    return(NULL)                                                       #
  }                                                                    #
  if(!is.numeric(samp_rat)){                                           #
    print(paste0("ERROR: 'samp_rat' is expected to be numeric. ",      #
                 "Returning NULL."))                                   #
    return(NULL)                                                       #
  }                                                                    #
  if(samp_rat <= 0 || samp_rat > 1){                                   #
    print(paste0("ERROR: 'samp_rat' is ", samp_rat, ", and not ",      #
                 "between 0 and 1 as is expected. ",                   #
                 "Returning NULL."))                                   #
    return(NULL)                                                       #
  }                                                                    #
                                                                       #
  samp_counts <- rep(0, length(counts))                                #
                                                                       #
  temp_cts <- round(sum(counts) * samp_rat / length(counts))           #
  temp_dif <- temp_cts - counts                                        #
  temp_inds <- which(temp_dif > 0, arr.ind = T)                        #
                                                                       #
  excess <- sum(temp_dif[temp_inds])                                   #
  temp_cts <- temp_cts +                                               #
    round(excess / (length(counts) - length(temp_inds)))               #
                                                                       #
  bagged_inds <- temp_inds                                             #
                                                                       #
  while(excess != 0){                                                  #
    temp_dif <- temp_cts - counts                                      #
    temp_inds <- setdiff(which(temp_dif > 0, arr.ind = T), bagged_inds)#
    bagged_inds <- append(bagged_inds, temp_inds)                      #
    excess <- sum(temp_dif[temp_inds])                                 #
    temp_cts <- temp_cts +                                             #
      round(excess / (length(counts) - length(bagged_inds)))           #
  }                                                                    #
                                                                       #
  samp_counts[bagged_inds] <- counts[bagged_inds]                      #
  samp_counts[-bagged_inds] <- temp_cts                                #
                                                                       #
  return(samp_counts)                                                  #
}                                                                      #
########################################################################

#//////////////////////////////////////////////////////////////////////#

############################################################################
######### Determines how many data points to sample from each bin. #########
######### That value being proportional to the counts of each bin. #########
### Bins which proportional sample is smaller than a minimum number are ####
# are either fully sampled or sampled that minimum (whichever is smaller). #
### The sampling excess (so as to achieve 'samp_rat') is proportionally ####
##################### removed from the remaining bins. #####################
## Inputs:                                                                ##
##        counts -> vector containing each bins counts;                   ##
##        samp_rat -> ratio of total counts to be sampled;                ##
##        min_counts -> minimum number of samples per bin.                ##
## Output:                                                                ##
##        integer amount of counts to be sampled from each bin.           ##
############################################################################
base_dist_sampling_counts <- function(counts, samp_rat, min_counts){       #
  if(!is.vector(counts) && !is.array(counts)){                             #
    print(paste0("ERROR: 'counts' is expected to either be a vector or",   #
                 " an array. Returning NULL."))                            #
    return(NULL)                                                           #
  }                                                                        #
  if(!is.numeric(counts)){                                                 #
    print(paste0("ERROR: 'counts' is expected to be numeric.",             #
                 " Returning NULL."))                                      #
    return(NULL)                                                           #
  }                                                                        #
  if(length(counts) == 0){                                                 #
    print("ERROR: 'counts' has length 0. Returning NULL.")                 #
    return(NULL)                                                           #
  }                                                                        #
  if(length(samp_rat) != 1){                                               #
    print(paste0("ERROR: 'samp_rat' should have length 1, and has length ",#
                 length(samp_rat), "instead. Returning NULL."))            #
    return(NULL)                                                           #
  }                                                                        #
  if(!is.numeric(samp_rat)){                                               #
    print(paste0("ERROR: 'samp_rat' is expected to be numeric. ",          #
                 "Returning NULL."))                                       #
    return(NULL)                                                           #
  }                                                                        #
  if(samp_rat <= 0 || samp_rat > 1){                                       #
    print(paste0("ERROR: 'samp_rat' is ", samp_rat, ", and not ",          #
                 "between 0 and 1 as is expected. Returning NULL."))       #
    return(NULL)                                                           #
  }                                                                        #
  if(!is.numeric(min_counts)){                                             #
    print(paste0("ERROR: 'min_counts' is expected to be numeric.",         #
                 " Returning NULL."))                                      #
    return(NULL)                                                           #
  }                                                                        #
  if(min_counts <= 0){                                                     #
    print(paste0("ERROR: 'min_counts' is ", min_counts, ", and is ",       #
                 "expected to be greater than 0. Returning NULL."))        #
    return(NULL)                                                           #
  }                                                                        #
                                                                           #
  #How many counts should be sampled in total                              #
  total_amount <- round(sum(counts) * samp_rat)                            #
                                                                           #
  #Amount to be proportionally distributed among bins                      #
  dist_amount <- 0                                                         #
                                                                           #
  while(dist_amount <= 0){                                                 #
                                                                           #
    #Output holding number of counts to be sampled from each Bin           #
    samp_counts <- counts - counts                                         #
                                                                           #
    #Which Bins don't have minimum counts                                  #
    flag_ind <- rep(FALSE, length(counts))                                 #
    ini_inds <- which(counts <= min_counts, arr.ind = T)                   #
    flag_ind[ini_inds] <- TRUE                                             #
                                                                           #
    #Those Bins are totally sampled                                        #
    samp_counts[ini_inds] <- counts[ini_inds]                              #
                                                                           #
    #The remaining Bins are at least sampled by the minimum amount         #
    samp_counts[-ini_inds] <- min_counts                                   #
                                                                           #
    #How many counts can still be sampled from each of the remaining bins  #
    rem_per_bin <- counts[-ini_inds] - min_counts                          #
                                                                           #
    #How many total counts can still be sampled                            #
    rem_counts <- sum(rem_per_bin)                                         #
                                                                           #
    #How many total counts are already attributed                          #
    ini_amount <- sum(samp_counts)                                         #
                                                                           #
    #How many counts are left to be distributed                            #
    dist_amount <- total_amount - ini_amount                               #
                                                                           #
    #Checking whether min_counts have to be adjusted                       #
    if(dist_amount <=0){                                                   #
      min_counts <- min_counts / 2                                         #
    }                                                                      #
  }                                                                        #
                                                                           #
  #Update sampling ratio                                                   #
  updt_samp_rat <- dist_amount / rem_counts                                #
                                                                           #
  #How many counts are sampled from each of the remaining bins             #
  rem_samp <- ceiling(rem_per_bin * updt_samp_rat)                         #
                                                                           #
  #Update output                                                           #
  samp_counts[-ini_inds] <- samp_counts[-ini_inds] + rem_samp              #
                                                                           #
  return(samp_counts)                                                      #
}                                                                          #
############################################################################


#//////////////////////////////////////////////////////////////////////////#


#######################################################################
############# Offsets value range of an array to ]0; +Inf[ ############
#######################################################################
myShift <- function(data){                                            #
  if(!is.vector(data) && !is.array(data)){                            #
    print(paste0("ERROR: 'data' is expected to either be a vector",   #
                 " or an array. Returning NULL."))                    #
    return(NULL)                                                      #
  }                                                                   #
  if(!is.numeric(data)){                                              #
    print("ERROR: 'data' is expected to be numeric. Returning NULL.") #
    return(NULL)                                                      #
  }                                                                   #
  if(length(data) == 0){                                              #
    print("ERROR: 'data' has length 0. Returning NULL.")              #
    return(NULL)                                                      #
  }                                                                   #
  dat <- data[which(!is.infinite(data), arr.ind = TRUE)]              #
  mins <- myMins(dat)                                                 #
  old_min <- mins[1]                                                  #
  new_min <- mins[2]                                                  #
  return(data - old_min + new_min)                                    #
}                                                                     #
#######################################################################

#/////////////////////////////////////////////////////////////////////#

#######################################################################
######## Offsets value range of "new_data" to that of "source" ########
#######################################################################
myShiftBack <- function(new_data, source){                            #
  if(!is.vector(source) && !is.array(source)){                        #
    print(paste0("ERROR: 'source' is expected to either be a vector ",#
                 "or an array. Returning NULL."))                     #
    return(NULL)                                                      #
  }                                                                   #
  if(!is.numeric(source)){                                            #
    print(paste0("ERROR: 'source' is expected to be numeric. ",       #
                 "Returning NULL."))                                  #
    return(NULL)                                                      #
  }                                                                   #
  if(length(source) == 0){                                            #
    print("ERROR: 'source' has length 0. Returning NULL.")            #
    return(NULL)                                                      #
  }                                                                   #
  if(!is.vector(new_data) && !is.array(new_data)){                    #
    print(paste0("ERROR: 'new_data' is expected to either be a vector",
                 " or an array. Returning NULL."))                    #
    return(NULL)                                                      #
  }                                                                   #
  if(!is.numeric(new_data)){                                          #
    print(paste0("ERROR: 'new_data' is expected to be numeric.",      #
                 "Returning NULL."))                                  #
    return(NULL)                                                      #
  }                                                                   #
  if(length(new_data) == 0){                                          #
    print("ERROR: 'new_data' has length 0. Returning NULL.")          #
    return(NULL)                                                      #
  }                                                                   #
  dat <- source[which(!is.infinite(source), arr.ind = TRUE)]          #
  mins <- myMins(dat)                                                 #
  old_min <- mins[1]                                                  #
  new_min <- mins[2]                                                  #
  return(new_data + old_min - new_min)                                #
}                                                                     #
#######################################################################

#/////////////////////////////////////////////////////////////////////#

#######################################################################
#### Determines the minimum of "data" and calculates a new strictly ###
##### positive minimum. Returns a vector containing both minima. ######
#######################################################################
myMins <- function(data){                                             #
  if(!is.vector(data) && !is.array(data)){                            #
    print(paste0("ERROR: 'data' is expected to either be a vector or",#
                 " an array. Returning NULL."))                       #
    return(NULL)                                                      #
  }                                                                   #
  if(!is.numeric(data)){                                              #
    print("ERROR: 'data' is expected to be numeric. Returning NULL.") #
    return(NULL)                                                      #
  }                                                                   #
  if(length(data) == 0){                                              #
    print("ERROR: 'data' has length 0. Returning NULL.")              #
    return(NULL)                                                      #
  }                                                                   #
  old_min <- min(data, na.rm = TRUE)                                  #
  R <- max(data, na.rm = TRUE) - old_min                              #
  data <- data - old_min                                              #
  new_min <- min(data[which(data > 0, arr.ind = TRUE)], na.rm = TRUE) #
  new_min <- new_min / R^2                                            #
  return(c(old_min, new_min))                                         #
}                                                                     #
#######################################################################

#/////////////////////////////////////////////////////////////////////#

#######################################################################
# Given numerical "source_data" and character "method", this function #
#  applies function "myShift" to "source_data" followed by a another  #
#                  function selected by "method".                     #
#######################################################################
transf_by_method <- function(source_data, method){                    #
  if(!is.vector(source_data) && !is.array(source_data)){              #
    print(paste0("ERROR: 'source_data' is expected to either be a ve",#
                 "ctor or an array. Returning NULL."))                #
    return(NULL)                                                      #
  }                                                                   #
  if(!is.numeric(source_data)){                                       #
    print(paste0("ERROR: 'source_data' is expected to be numeric. ",  #
          "Returning NULL."))                                         #
    return(NULL)                                                      #
  }                                                                   #
  if(length(source_data) == 0){                                       #
    print("ERROR: 'source_data' has length 0. Returning NULL.")       #
    return(NULL)                                                      #
  }                                                                   #
  if(!is.character(method)){                                          #
    print(paste0("ERROR: 'method' is expected to be character. ",     #
                 "Returning NULL."))                                  #
    return(NULL)                                                      #
  }                                                                   #
  if(method != 'None' && method != 'Norm' &&                          #
     method != 'Log10' && method != 'Ln') {                           #
    print(paste0("ERROR: 'method' is expected to be 'None', 'Norm', ",# 
                 "'Log10' or 'Ln'. Returning NULL."))                 #
    return(NULL)                                                      #
  }                                                                   #
                                                                      #
  shifted_data <- myShift(source_data)                                #
                                                                      #
  switch(method,                                                      #
         None = shifted_data,                                         #
         Norm = shifted_data / mean(shifted_data, na.rm = T),         #
         Log10 = myShift(log10(shifted_data)),                        #
         Ln = myShift(log(shifted_data)))                             #
}                                                                     #
#######################################################################

#/////////////////////////////////////////////////////////////////////#

#######################################################################
#### Given numerical "new_data" generated from "source_data", this ####
#### function reverses the transformation selected by "method" and ####
############### applied by function "transf_by_method". ###############
#######################################################################
inv_transf_by_method <- function(new_data, source_data, method){      #
  if(!is.vector(source_data) && !is.array(source_data)){              #
    print(paste0("ERROR: 'source_data' is expected to either be a ve",#
                 "ctor or an array. Returning NULL."))                #
    return(NULL)                                                      #
  }                                                                   #
  if(!is.numeric(source_data)){                                       #
    print(paste0("ERROR: 'source_data' is expected to be numeric.",   #
                 "Returning NULL."))                                  #
    return(NULL)                                                      #
  }                                                                   #
  if(length(source_data) == 0){                                       #
    print("ERROR: 'source_data' has length 0. Returning NULL.")       #
    return(NULL)                                                      #
  }                                                                   #
  if(!is.vector(new_data) && !is.array(new_data)){                    #
    print(paste0("ERROR: 'new_data' is expected to either be a ve",   #
                 "ctor or an array. Returning NULL."))                #
    return(NULL)                                                      #
  }                                                                   #
  if(!is.numeric(new_data)){                                          #
    print(paste0("ERROR: 'new_data' is expected to be numeric.",      #
          "Returning NULL."))                                         #
    return(NULL)                                                      #
  }                                                                   #
  if(length(new_data) == 0){                                          #
    print("ERROR: 'new_data' has length 0. Returning NULL.")          #
    return(NULL)                                                      #
  }                                                                   #
  if(!is.character(method)){                                          #
    print(paste0("ERROR: 'method' is expected to be character.",      #
                 "Returning NULL."))                                  #
    return(NULL)                                                      #
  }                                                                   #
  if(method != 'None' && method != 'Norm' &&                          #
     method != 'Log10' && method != 'Ln') {                           #
    print(paste0("ERROR: 'method' is expected to be 'None', 'Norm', ",#
                 "'Log10' or 'Ln'. Returning NULL."))                 #
    return(NULL)                                                      #
  }                                                                   #
                                                                      #
  shifted_data <- myShift(source_data)                                #
                                                                      #
  switch(method,                                                      #
    None = myShiftBack(new_data, source_data),                        #
    Norm = myShiftBack(new_data * mean(shifted_data, na.rm = T),      #
                       source_data),                                  #
    Log10 = myShiftBack(10^myShiftBack(new_data, log10(shifted_data)),#
                        source_data),                                 #
    Ln = myShiftBack(exp(1)^myShiftBack(new_data, log(shifted_data)), #
                     source_data))                                    #
}                                                                     #
#######################################################################