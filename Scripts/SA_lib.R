###########################
# Spaxel Analysis Library #
###########################

require(keras)
require(rapportools)

##############################################################
##############################################################
# Calculates the Pearson Correlation between two data arrays #
##############################################################
##############################################################
pearsonCC <- function(x, y){                                 #
                                                             #
  return(cov(x, y) / (sd(x) * sd(y)))                        #
}                                                            #
##############################################################
##############################################################

##//////////////////////////////////////////////////////////##

###############################################################
###############################################################
########## Determines if k conforms to expected types #########
# Returns: 1 if NA or NaN; 2 if NULL; 3 if numeric or boolean # 
###############################################################
###############################################################
null_type <- function(k){                                     #
                                                              #
  type_flag <- NULL                                           #
                                                              #
  if(is.na(k))                                                #
    type_flag <- 1                                            #
  if(is.null(k))                                              #
    type_flag <- 2                                            #
  if(is.finite(k) || is.boolean(k) || is.infinite(k))         #
    type_flag <- 3                                            #
                                                              #
  if(is.null(type_flag)){                                     #
    print(paste0("ERROR: Unexpected k value. This function ", #
                 "expects k to be real, NA, NaN, or NULL."))  #
  }                                                           #
  return(type_flag)                                           #
}                                                             #
###############################################################
###############################################################

##///////////////////////////////////////////////////////////##

###############################################################################
###############################################################################
######################### Checks if an array x is N-d #########################
############################## Returns a boolean ##############################
###############################################################################
###############################################################################
is_Nd_array <- function(x, N){                                                #
                                                                              #
  ds <- length(dim(x))                                                        #
                                                                              #
  if(ds != N){                                                                #
    print(paste0("ERROR: The provided input x has ", ds, " dimensions.",      #
                 "This function expects ", N, " dimensions. Returning NULL."))#
    return(FALSE)                                                             #
  }                                                                           #
  else return(TRUE)                                                           #
}                                                                             #
###############################################################################
###############################################################################

##///////////////////////////////////////////////////////////////////////////##

################################################################################
################################################################################
######### Given a 3-d array x and position of its spatial dimensions ###########
###### Inputs: x - 3-d array; a - 1st spatial dimension; b - 2nd spatial d #####
########## Returns the [1] spatial and [2] spectral size of the array ##########
########## or NULL if the array is not 3-d, or a and b are not valid. ##########
################################################################################
################################################################################
array_3d_sizes <- function(x, a, b){                                           #
                                                                               #
  if(!is_Nd_array(x, 3)){                                                      #
    print("The array is not 3d, returning NULL.")                              #
    return(c(NULL, NULL))                                                      #
  }                                                                            #
                                                                               #
  if(a < 1 || a > 3 || b < 1 || b > 3){                                        #
    print(paste0("ERROR: The provided spatial dimensions a = ", a, " and b = ",#
                 b, " do not conform with the expectations. This function exp",#
                 "ects a and b to be either 1, 2 or 3. Returning NULL."))      #
    return(c(NULL, NULL))                                                      #
  }                                                                            #
                                                                               #
  if(a == b){                                                                  #
    print(paste0("ERROR: The provided spatial dimensions a = ", a, " and b = ",#
                 b, " do not conform with the expectations. This function exp",#
                 "ects a != b. Returning NULL."))                              #
    return(c(NULL, NULL))                                                      #
  }                                                                            #
                                                                               #
  c <- setdiff(1:3, c(a, b))                                                   #
  dims <- dim(x)                                                               #
  spatial <- dims[a] * dims[b]                                                 #
  spectral <- dims[c]                                                          #
                                                                               #
  return(c(spatial, spectral))                                                 #
}                                                                              #
################################################################################
################################################################################

##////////////////////////////////////////////////////////////////////////////##

######################################################################
######################################################################
###### Counts how many positions in an array have information, #######
####### where information is considered to be any value but k ########
#### Inputs: x - N-d array; k - value (ex. 0, 1, NA, NULL, FALSE) ####
######################### Returns an integer #########################
######################################################################
######################################################################
info_size <- function(x, k){                                         #
                                                                     #
  type_flag <- null_type(k)                                          #
  switch(type_flag,                                                  #
         count <- length(x[which(!is.na(x), arr.ind = TRUE)]),       #
         count <- length(x[which(!is.null(x), arr.ind = TRUE)]),     #
         count <- length(x[which(x != k, arr.ind = TRUE)]))          #
                                                                     #
  return(count)                                                      #
}                                                                    #
######################################################################
######################################################################

##//////////////////////////////////////////////////////////////////##

######################################################################
######################################################################
############ Counts how many spaxels in 3d array are full ############
## Inputs: x - 3-d array; a - 1st spatial dim; b - 2nd spatial dim; ##
####### k - value for lack of info (ex. 0, 1, NA, NULL, FALSE) #######
######################### Returns an integer #########################
######################################################################
######################################################################
full_spx_size <- function(x, a, b, k){                               #
                                                                     #
  temp <- array_3d_sizes(x, a, b)                                    #
  spatial_size <- temp[1]                                            #
  feature_size <- temp[2]                                            #
  count <- spatial_size                                              #
  rm(temp)                                                           #
                                                                     #
  y <- array_reshape(x, dim = c(spatial_size, feature_size))         #
                                                                     #
  type_flag <- null_type(k)                                          #
  switch(type_flag,                                                  #
         for(i in 1:spatial_size){                                   #
           j <- 1                                                    #
           break_flag <- TRUE                                        #
                                                                     #
           while(j <= feature_size && break_flag){                   #
             if(is.na(y[i, j])){                                     #
               break_flag <- FALSE                                   #
               count <- count - 1                                    #
             }                                                       #
             else j <- j + 1                                         #
           }                                                         #
         },                                                          #
                                                                     #
         for(i in 1:spatial_size){                                   #
           j <- 1                                                    #
           break_flag <- TRUE                                        #
                                                                     #
           while(j <= feature_size && break_flag){                   #
             if(is.null(y[i, j])){                                   #
               break_flag <- FALSE                                   #
               count <- count - 1                                    #
             }                                                       #
             else j <- j + 1                                         #
           }                                                         #
         },                                                          #
                                                                     #
         for(i in 1:spatial_size){                                   #
           j <- 1                                                    #
           break_flag <- TRUE                                        #
                                                                     #
           while(j <= feature_size && break_flag){                   #
             if(y[i, j] == k){                                       #
               break_flag <- FALSE                                   #
               count <- count - 1                                    #
             }                                                       #
             else j <- j + 1                                         #
           }                                                         #
         })                                                          #
                                                                     #
  return(count)                                                      #
}                                                                    #
######################################################################
######################################################################

##//////////////////////////////////////////////////////////////////##

######################################################################
######################################################################
########### Counts how many spaxels in an array are empty, ###########
## Inputs: x - 3-d array; a - 1st spatial dim; b - 2nd spatial dim; ##
####### k - value for lack of info (ex. 0, 1, NA, NULL, FALSE) #######
######################### Returns an integer #########################
######################################################################
######################################################################
null_spx_size <- function(x, a, b, k){                               #
                                                                     #
  temp <- array_3d_sizes(x, a, b)                                    #
  spatial_size <- temp[1]                                            #
  feature_size <- temp[2]                                            #
  count <- spatial_size                                              #
  rm(temp)                                                           #
                                                                     #
  y <- array_reshape(x, dim = c(spatial_size, feature_size))         #
                                                                     #
  type_flag <- null_type(k)                                          #
  switch(type_flag,                                                  #
         for(i in 1:spatial_size){                                   #
           j <- 1                                                    #
           break_flag <- TRUE                                        #
                                                                     #
           while(j <= feature_size && break_flag){                   #
             if(!is.na(y[i, j])){                                    #
               break_flag <- FALSE                                   #
               count <- count - 1                                    #
             }                                                       #
             else j <- j + 1                                         #
           }                                                         #
         },                                                          #
                                                                     #
         for(i in 1:spatial_size){                                   #
           j <- 1                                                    #
           break_flag <- TRUE                                        #
                                                                     #
           while(j <= feature_size && break_flag){                   #
             if(!is.null(y[i, j])){                                  #
               break_flag <- FALSE                                   #
               count <- count - 1                                    #
             }                                                       #
             else j <- j + 1                                         #
           }                                                         #
         },                                                          #
                                                                     #
         for(i in 1:spatial_size){                                   #
           j <- 1                                                    #
           break_flag <- TRUE                                        #
                                                                     #
           while(j <= feature_size && break_flag){                   #
             if(y[i, j] != k){                                       #
               break_flag <- FALSE                                   #
               count <- count - 1                                    #
             }                                                       #
             else j <- j + 1                                         #
           }                                                         #
         })                                                          #
                                                                     #
  return(count)                                                      #
}                                                                    #
######################################################################
######################################################################

##//////////////////////////////////////////////////////////////////##

####################################################################
####################################################################
######## Counts how many spaxels are neither empty nor full ########
# Inputs: x - 3-d array; a - 1st spatial dim; b - 2nd spatial dim; #
###### k - value for lack of info (ex. 0, 1, NA, NULL, FALSE) ######
######################## Returns an integer ########################
####################################################################
####################################################################
incomp_spx_size <- function(x, a, b, k){                           #
                                                                   #
  count <- array_3d_sizes(x, a, b)[1]                              #
  null_count <- null_spx_size(x, a, b, k)                          #
  full_count <- full_spx_size(x, a, b, k)                          #
  count <- count - null_count - full_count                         #
                                                                   #
  return(count)                                                    #
}                                                                  #
####################################################################
####################################################################

##////////////////////////////////////////////////////////////////##

#######################################################################
#######################################################################
# Creates a boolean map with location of null spaxels (empty spaxels) #
### Inputs: x - 2d array (first dimension with spatial information, ###
## the second with spectral information); k - value of lack of info ###
##################### (ex. 0, 1, NA, NULL, FALSE) #####################
### Returns a 2-d array with TRUE in the coordinates of null spaxels ##
#######################################################################
#######################################################################
null_spaxel_map <- function(x, k){                                    #
                                                                      #
  if(!is_Nd_array(x, 2)) return(NULL)                                 #
                                                                      #
  dims <- dim(x)                                                      #
  spatial_size <- dims[1]                                             #
  feature_size <- dims[2]                                             #
                                                                      #
  if(is.null(spatial_size)){                                          #
    print("ERROR: spatial dimensions are null. Returning NULL.")      #
    return(NULL)                                                      #
  }                                                                   #
  if(is.null(feature_size)){                                          #
    print("ERROR: spectral dimension is null. Returning NULL.")       #
    return(NULL)                                                      #
  }                                                                   #
                                                                      #
  ############################################                        #
  ### Will hold a map of null flux spaxels ###                        #
  ############################################                        #
  null_spx <- array(TRUE, dim = spatial_size)#                        #
  ############################################                        #
                                                                      #
  ###############################################################     #
  #################Prelim ID NULL SPAXELS #######################     #
  #### defined as spaxels that are k for all wavelengths (l) ####     #
  ###############################################################     #
  type_flag <- null_type(k)                                     #     #
  switch(type_flag,                                             #     #
         for(i in 1:spatial_size){                              #     #
           j <- 1                                               #     #
           null_flag <- TRUE                                    #     #
                                                                #     #
           while(j <= feature_size && null_flag){               #     #
             if(!is.na(x[i, j])){                               #     #
               null_flag <- FALSE                               #     #
               null_spx[i] <- FALSE                             #     #
             }                                                  #     #
             j <- j + 1                                         #     #
           }                                                    #     #
         },                                                     #     #
                                                                #     #
         for(i in 1:spatial_size){                              #     #
           j <- 1                                               #     #
           break_flag <- TRUE                                   #     #
                                                                #     #
           while(j <= feature_size && null_flag){               #     #
             if(!is.null(x[i, j])){                             #     #
               null_flag <- FALSE                               #     #
               null_spx[i] <- FALSE                             #     #
             }                                                  #     #
             j <- j + 1                                         #     #
           }                                                    #     #
         },                                                     #     #
                                                                #     #
         for(i in 1:spatial_size){                              #     #
           j <- 1                                               #     #
           null_flag <- TRUE                                    #     #
                                                                #     #
           while(j <= feature_size && null_flag){               #     #
             if(x[i, j] != k){                                  #     #
               null_flag <- FALSE                               #     #
               null_spx[i] <- FALSE                             #     #
             }                                                  #     #
             j <- j + 1                                         #     #
           }                                                    #     #
         })                                                     #     #
  ###############################################################     #
                                                                      #
  return(null_spx)                                                    #
}                                                                     #
#######################################################################
#######################################################################

##///////////////////////////////////////////////////////////////////##

###############################################################################
###############################################################################
#### Creates a boolean map with an estimate location of TRUE null spaxels. ####
#### Those are spaxels that are empty because there is nothing coming from ####
#### that particular spatial region, instead of simply being a result of a ####
######### low event counts. This function takes a boolean map of null #########
######### spaxels and filters out fake nulls by analysing iteratively #########
######################### each spaxels neighbourhood. #########################
########### Inputs: x - 2d boolean array of signaling null spaxels; ###########
######## sBSize - half-length of square box centered around pixel x, y ########
############# Returns a 2-d array with TRUE for true null spaxels #############
###############################################################################
###############################################################################
true_null_spaxel_map <- function(z, sBSize, m){                               #
                                                                              #
  if(!is_Nd_array(z, 2)) return(NULL)                                         #
  dimz <- dim(z)                                                              #
                                                                              #
  spatial_zize <- dimz[1] * dimz[2] / m                                       #
                                                                              #
  #########################################################################   #
  # To eliminate probable 'fake' nulls (nulls that have less than 2       #   #
  # null l1-neighbours)                                                   #   #
  # 1- Make a map of preliminary nulls                                    #   #
  #                                                                       #   #
  # For each pixel:                                                       #   #
  # 2- Check if there are non-nulls, in a k-sized box around a pixel      #   #
  #    (x,y) in every direction ([x:x+k],[x-k:x],[y:y+k],[y-k:y])         #   #
  #                                                                       #   # 
  # 3- If there are non-nulls in every direction, that pixel is           #   #
  #    converted to non-null                                              #   #
  #                                                                       #   #
  # After all pixels have been checked                                    #   #
  # 4- Update map of nulls and repeat 2:3 until there are no changes      #   #
  #    to the map                                                         #   #
  #########################################################################   #
  for(mult in 1:m){                                                       #   #
                                                                          #   #
    curr_neighb_map <- matrix(0, nrow = dimz[2], ncol = dimz[1])          #   #
    temp_neighb_map <- matrix(0, nrow = dimz[2], ncol = dimz[1])          #   #
                                                                          #   #
    dif_sum <- 1                                                          #   #
                                                                          #   #
    while(dif_sum > 0 ){                                                  #   #
      for(x in 1:dimz[2]){                                                #   #
        for(y in 1:dimz[1]){                                              #   #
                                                                          #   #
          # left, x - sBSize                                              #   #
          # right, x + sBSize                                             #   #
          # north, y - sBSize                                             #   #
          # south, y + sBSize                                             #   #
          dir_flag <- c(0, 0, 0, 0)                                       #   #
                                                                          #   #
          q <- (mult - 1) * spatial_zize + (y - 1) * dimz[2] + x          #   #
                                                                          #   #
          # Determining if there are non-null neighbours                  #   #
          # in all directions of a null                                   #   #
          if(z[q] == TRUE){                                               #   #
                                                                          #   #
            xVal <- (x - sBSize):(x + sBSize)                             #   # 
            yVal <- (y - sBSize):(y + sBSize)                             #   #
                                                                          #   #
            for(k in xVal){                                               #   #
                                                                          #   #
              # if we're not beyond the map horizontal borders            #   #
              if(k <= dimz[2] && k > 0){                                  #   #
                for(j in yVal){                                           #   #
                                                                          #   #
                  # if we're not beyond the map vertical borders          #   #
                  if(j <= dimz[1] && j > 0){                              #   #
                                                                          #   #
                    p <- (mult - 1) * spatial_zize + (j - 1) * dimz[2] + k#   #
                                                                          #   #
                    if(z[p] == FALSE){                                    #   #
                      if(k < x){ dir_flag[1] <- 1 }                       #   #
                      if(k > x){ dir_flag[2] <- 1 }                       #   #
                      if(j < y){ dir_flag[3] <- 1 }                       #   #
                      if(j > y){ dir_flag[4] <- 1 }                       #   #
                    }                                                     #   #
                  }                                                       #   #
                }                                                         #   #
              }                                                           #   #
            }                                                             #   #
          }                                                               #   #
          curr_neighb_map[x, y] <- length(dir_flag[which(dir_flag != 0)]) #   #
        }                                                                 #   #
      }                                                                   #   #
                                                                          #   #
      dif_sum <- abs(sum(curr_neighb_map - temp_neighb_map, na.rm = TRUE))#   #
      temp_neighb_map <- curr_neighb_map                                  #   #
                                                                          #   #
      temp_Null <- z                                                      #   #
      temp_Null[which(curr_neighb_map == 4)] <- FALSE                     #   #
      z <- temp_Null                                                      #   #
    }                                                                     #   #
  }                                                                       #   #
  print('True Zero Map has been constructed.')                            #   #
  #########################################################################   #
                                                                              #
  return(z)                                                                   #
}                                                                             #
###############################################################################
###############################################################################