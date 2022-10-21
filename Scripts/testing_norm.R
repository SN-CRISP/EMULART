rm(list = ls())
options(digits = 2)
require(FITSio)

source_path <- "/media/joaomfras/TOSHIBA EXT/EmulaRT/Sources/FITS/t9-0.05/face/"

input_1e8_path <- paste0(source_path, "sphShell_t9-0.05_aniso_pp1e8_i00", 
                         "_secondarydirect.fits")
input_1e7_path <- paste0(source_path, "sphShell_t9-0.05_aniso_pp1e7_i00", 
                         "_secondarydirect.fits")
input_1e6_path <- paste0(source_path, "sphShell_t9-0.05_aniso_pp1e6_i00", 
                         "_secondarydirect.fits")
input_1e5_path <- paste0(source_path, "sphShell_t9-0.05_aniso_pp1e5_i00", 
                         "_secondarydirect.fits")
input_1e4_path <- paste0(source_path, "sphShell_t9-0.05_aniso_pp1e4_i00", 
                         "_secondarydirect.fits")

scale_path <- "/media/joaomfras/TOSHIBA EXT/EmulaRT/Sources/Ref_Scale.dat"

output_path <- "/media/joaomfras/TOSHIBA EXT/EmulaRT/Results/"

#Loading files
input_1e8 <- readFITS(input_1e8_path)
input_1e7 <- readFITS(input_1e7_path)
input_1e6 <- readFITS(input_1e6_path)
input_1e5 <- readFITS(input_1e5_path)
input_1e4 <- readFITS(input_1e4_path)

dims =  dim(input_1e8$imDat) 
features <- 40:dims[3]
feature_size <- length(features)
myScale <- read.csv(scale_path, header = FALSE, sep = " ")[1:dims[3], 1] 
myScale <- myScale[features]

#Selecting relevant spectral features
data8 <- input_1e8$imDat[, , features]
data7 <- input_1e7$imDat[, , features]
data6 <- input_1e6$imDat[, , features]
data5 <- input_1e5$imDat[, , features]
data4 <- input_1e4$imDat[, , features]
rm(input_1e8, input_1e7, input_1e6, input_1e5, input_1e4)

#Creating integrated SEDs from the sum of all pixels at each wavelength
sed8 <- apply(data8, 3, sum, na.rm = TRUE)
sed7 <- apply(data7, 3, sum, na.rm = TRUE)
sed6 <- apply(data6, 3, sum, na.rm = TRUE)
sed5 <- apply(data5, 3, sum, na.rm = TRUE)
sed4 <- apply(data4, 3, sum, na.rm = TRUE)

#Calculating integrated SEDs residuals
res7 <- abs((sed8 - sed7) / sed8) * 100
res6 <- abs((sed8 - sed6) / sed8) * 100
res5 <- abs((sed8 - sed5) / sed8) * 100
res4 <- abs((sed8 - sed4) / sed8) * 100
res7[which(is.infinite(res7), arr.ind = T)] <- NA
res6[which(is.infinite(res6), arr.ind = T)] <- NA
res5[which(is.infinite(res5), arr.ind = T)] <- NA
res4[which(is.infinite(res4), arr.ind = T)] <- NA

#Checking the total photon flux content of each cube
TP8 <- sum(data8, na.rm = TRUE)
TP7 <- sum(data7, na.rm = TRUE)
TP6 <- sum(data6, na.rm = TRUE)
TP5 <- sum(data5, na.rm = TRUE)
TP4 <- sum(data4, na.rm = TRUE)
print(TP8)
print(TP7)
print(TP6)
print(TP5)
print(TP4)

#Re-scaling data since total flux is the same for each cube
reg_data7 <- data7 * 10^(7-8)
reg_data6 <- data6 * 10^(6-8)
reg_data5 <- data5 * 10^(5-8)
reg_data4 <- data4 * 10^(4-8)

#Creating new integrated SEDs
reg_sed7 <- apply(reg_data7, 3, sum, na.rm = TRUE)
reg_sed6 <- apply(reg_data6, 3, sum, na.rm = TRUE)
reg_sed5 <- apply(reg_data5, 3, sum, na.rm = TRUE)
reg_sed4 <- apply(reg_data4, 3, sum, na.rm = TRUE)

#Calculating integrated SEDs residuals
reg_res7 <- abs((sed8 - reg_sed7) / sed8) * 100
reg_res6 <- abs((sed8 - reg_sed6) / sed8) * 100
reg_res5 <- abs((sed8 - reg_sed5) / sed8) * 100
reg_res4 <- abs((sed8 - reg_sed4) / sed8) * 100
reg_res7[which(is.infinite(reg_res7), arr.ind = T)] <- NA
reg_res6[which(is.infinite(reg_res6), arr.ind = T)] <- NA
reg_res5[which(is.infinite(reg_res5), arr.ind = T)] <- NA
reg_res4[which(is.infinite(reg_res4), arr.ind = T)] <- NA

#Printing SEDs (before and after total energy re-scaling)
ax_x <- myScale
ax_y_pre <- 10^seq(-17, -4, 1)
ax_y_post <- 10^seq(-13, -4, 1)
ax_y_pre_res <- seq(0, 50, 10)
ax_y_post_res <- seq(88, 102, 2)

pdf(paste0(output_path,'test_SEDs_pre.pdf'), width = 10, height = 10)
par(mar=c(5.1, 5.6, 5.1, 1))

matplot (log10(myScale), cbind(log10(sed8), log10(sed4), log10(sed5), 
                               log10(sed6), log10(sed7)), type = "b", 
         pch = 1:5, col = 1:5, 
         xlab = expression(paste("Wavelength [", mu, "m]")),
         ylab = "Flux Density [W/m²]", cex.lab = 2, axes = F)
axis(1, at = log10(ax_x), labels = sprintf(ax_x, fmt = '%#.1f'), cex.axis = 1.5)
axis(2, at = log10(ax_y_post), labels = ax_y_post, cex.axis = 1.5)
legend(2.6, -5, legend = c("1e8", "1e4", "1e5", "1e6", "1e7"),
       col = 1:5, pch = 1:5, cex = 1.5)
title("Integrated SEDs", line = 2.5)
dev.off()

pdf(paste0(output_path,'test_SEDs_pre_RES.pdf'), width = 10, height = 5)
par(mar=c(5.1, 5.6, 5.1, 1))

matplot (log10(myScale), cbind(res4, res5, res6, res7), type = "b", 
         pch = 2:5, col = 2:5, 
         xlab = expression(paste("Wavelength [", mu, "m]")),
         ylab = "Normalized Residuals [%]", cex.lab = 2, axes = F,
         ylim = c(0, 50))
axis(1, at = log10(ax_x), labels = sprintf(ax_x, fmt = '%#.1f'), cex.axis = 1.5)
axis(2, at = ax_y_pre_res, labels = ax_y_pre_res, cex.axis = 1.5)
legend(2.7, 50, legend = c("1e4", "1e5", "1e6", "1e7"),
       col = 2:5, pch = 2:5, cex = 1.5)
title("Integrated SEDs Residuals", line = 2.5)
dev.off()

pdf(paste0(output_path,'test_SEDs_post.pdf'), width = 10, height = 10)
par(mar=c(5.1, 5.6, 5.1, 1))

matplot (log10(myScale), cbind(log10(sed8), log10(reg_sed4), log10(reg_sed5), 
                               log10(reg_sed6), log10(reg_sed7)), type = "b", 
         pch = 1:5, col = 1:5, 
         xlab = expression(paste("Wavelength [", mu, "m]")),
         ylab = "Flux Density [W/m²]", cex.lab = 2, axes = F)
axis(1, at = log10(ax_x), labels = sprintf(ax_x, fmt = '%#.1f'), cex.axis = 1.5)
axis(2, at = log10(ax_y_pre), labels = ax_y_pre, cex.axis = 1.5)
legend(2.6, -6, legend = c("1e8", "1e4", "1e5", "1e6", "1e7"),
       col = 1:5, pch = 1:5, cex = 1.5)
title("Integrated SEDs", line = 2.5)

dev.off()

pdf(paste0(output_path,'test_SEDs_post_RES.pdf'), width = 10, height = 5)
par(mar=c(5.1, 5.6, 5.1, 1))

matplot (log10(myScale), cbind(reg_res4, reg_res5, reg_res6, reg_res7), 
         type = "b", pch = 2:5, col = 2:5, 
         xlab = expression(paste("Wavelength [", mu, "m]")),
         ylab = "Normalized Residuals [%]", cex.lab = 2, axes = F, 
         ylim = c(88, 102))
axis(1, at = log10(ax_x), labels = sprintf(ax_x, fmt = '%#.1f'), cex.axis = 1.5)
axis(2, at = ax_y_post_res, labels = ax_y_post_res, cex.axis = 1.5)
legend(2.7, 120, legend = c("1e4", "1e5", "1e6", "1e7"),
       col = 2:5, pch = 2:5, cex = 1.5)
title("Integrated Reg SEDs Residuals", line = 2.5)
dev.off()