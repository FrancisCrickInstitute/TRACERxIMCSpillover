################## SPILLOVER ESTIMATION FOR TRACERx WITH ADAPTIVE BINNING ########################

### Original script authored by Vito Zanotelli et al. 
### Original script available at https://github.com/BodenmillerGroup/cyTOFcompensation/blob/master/scripts/imc_generatespillmat_long.Rmd
### Adaptations to the script made by Emma Colliver (2021) 
### Adapted script made available at https://github.com/FrancisCrickInstitute/TRACERxIMCSpilllover 
### Any adaptations to individual functions or newly written functions are explicitly highlighted below
### Original functions not invoked in the adapted scripts are removed for brevity

### Aim
##This script shows how to estimate spillover from single metal spots on an agarose coated slide.
##Each spot should be imaged with a single acquisition. The name of the acquisition should be the metal that is used:
##E.g. PanormaA_1_Yb176_23.txt

### Parameters
#' @param fol_ss folder containing .txt acquisitions of IMC single stains
#' @param ssmetals_from_fn logical, Are the single stains correctly named xxx_x_metal_x.txt? (E.g. Dy161 1-1000_8_Dy161_8.txt
#' @param ssmass Vector of masses of the single stains used. Required if ssmetals_from_file_fn is False
#' @param remove_incorrect_bc Remove barcodes not matching the filename single stain annotation (requires either ssmetals_from_fn=T or fn2ssmetal )
#' @param minevents Minimal number of events (after debarcoding) that need to be present in a single stain in order that a spillover estimation is performed
#' @param ... Optional parameters will be passed to CATALYST::computeSpillmat 
#' @return a list containing the spillover matrix (sm), the data (data) and the debarcoded Catalyst object (re)

### Loading libraries
library(CATALYST)
library(data.table)
library(ggplot2)
library(flowCore)
library(dplyr)
library(dtplyr)
library(stringr)
library(ggpmisc)

### Sourcing helper functions
source('/PATH/TO/spillover_imc_helpers-release.R') #replace with path to spillover IMC helper functions

## Configuring the filepaths
# list of folders that contain each a complete single stain acquisition (e.g. in case that one wants to run and compare multiple single stains from different days)
fols_ss = c('/PATH/TO/SINGLE/STAIN/ACQUISITION/DATA/')
# output folder
fol_out = '/PATH/TO/OUTPUT/DATA/'
# name prefix for all output
prefix ='spillover_matrix_'

### Loading single stains
list_img_ss <-lapply(fols_ss, load_ss_fol)
names(list_img_ss) <- fols_ss

### Adapting the column names to be recognized metal names by CATALYST 
##CATALYST needs to have the metal names in the format (METAL)(MASS)Di
list_img_ss = lapply(list_img_ss, function(x) lapply(x, fixnames))
dats_raw = lapply(list_img_ss, imglist2dat)

### Extracting the single stain masses from the acquisition name 
for (dat in dats_raw){
  dat[, metal:= strsplit(.BY[[1]], '_')[[1]][3],by=file]
  dat[, mass:= as.numeric(str_extract_all(.BY[[1]], "[0-9]+")[[1]]),by=metal]
}

### Removing less bright pixels - NEW FUNCTION
dats_raw_bright = lapply(dats_raw, retain_brightest_signal)

############## ADAPTIVE APPROACH TO BINNING - ADAPTATION TO THE ORIGINAL SCRIPT

### Defining a vector with the number of pixels to bin over for each marker 
marker_pixel_bins <- data.frame(metal = unique(dat$metal), pixel_bin = NA)

### Loop over files and channels so that all channels and files have at >=250 
### median counts for that channel 
### Implemented to ensure enough median pixel counts for accurate spillover estimation while
### not 'over-binning' which could introduce noise
### Adheres to recommendation from original script 
### Original script: "If the median per-pixel intensities are to[o] low, it could be worth to sum up 
### some consecutive pixels to get a better accuracy for the estimation. 
### This is valid because for segmentation based quantitative image analysis usually anyways pixels 
### are aggregated. If the binning is choosen to[o] big, there is however a potential accumulation 
### of background noise."

for (i in unique(dat$metal)) {
  ### Calculate per-file medians
  dats_agg_sum = rbindlist(lapply(dats_raw_bright, calc_file_medians),idcol = T)
  dats_agg_sum_signal = dats_agg_sum[which(substr(gsub("Di","",dats_agg_sum$variable),3,5) == substr(dats_agg_sum$metal,3,5)),]
  pixelbin = 1 # initial state
  while (dats_agg_sum_signal[which(dats_agg_sum_signal$metal == i),]$med < 250) {
    pixelbin = pixelbin + 1
    dats_agg <- lapply(dats_raw_bright, function(x) aggregate_pixels(x, n=pixelbin))
    dats_agg_sum = rbindlist(lapply(dats_agg, calc_file_medians), idcol = T)
    dats_agg_sum_signal = dats_agg_sum[which(substr(gsub("Di","",dats_agg_sum$variable),3,5) == substr(dats_agg_sum$metal,3,5)),]
  }
    marker_pixel_bins[which(marker_pixel_bins$metal == i),]$pixel_bin = pixelbin 
}


### Applies number of pixels the aggregation should happen on a per-marker basis ('adaptive' approach) - NEW
dats_agg <- aggregate_pixels_by_marker(dats_raw_bright, m = as.character(marker_pixel_bins$metal), n = as.numeric(marker_pixel_bins$pixel_bin))
dats_agg_list <- list(fols_ss = dats_agg)

### CATALYST based compensation 
### Estimating the spillover
### To estimate the spillover, the (aggregated) pixel values are first debarcoded using CATALYST, treating them like single cells. This step acts as a quality filter to remove background/noisy/weak pixels as well as pixels with artefacts (e.g. specles with strong signal in many channels).
### If the true metal was correctly encoded in the filename, the 'remove_incorrect_bc' option will check the debarcoding and remove events assigned to the wrong barcode.
### Then this identified, strong single stain pixels will be used for the spillover estimation.

custom_isotope_list <- c(CATALYST::isotope_list, list(ArAr=80)) #adding Argon to the list

res = lapply(dats_agg_list, function(x) re_from_dat(x,
                                               ss_ms=x[!is.na(mass), unique(mass)],
                                               minevents = 20, #changed from default 40
                                               correct_bc = T))
sms = lapply(res, function(x) computeSpillmat(x))

### Saving the spillover matrices 
for (i in seq_along(sms)){
  outname = file.path(fol_out, paste0(prefix, basename(fols_ss[i]),'_sm_adaptive.csv')) #e.g. filename
  write.csv(sms[[i]],file = outname)
}

### Visualization of the spillover matrix 
for (i in seq_along(sms)){
  print(names(dats_agg_list)[i])
  ss_ms = dats_agg_list[[i]][!is.na(mass), unique(mass)]
  p = CATALYST::plotSpillmat(ss_ms,sms[[i]],isotope_list=custom_isotope_list)
  print(p)
}