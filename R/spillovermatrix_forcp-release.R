################## SPILLOVER ESTIMATION FOR TRACERx WITH ADAPTIVE BINNING ########################

### Original script authored by Vito Zanotelli et al. 
### Original script available at https://github.com/BodenmillerGroup/cyTOFcompensation/blob/master/scripts/imc_adaptsm.Rmd
### Adaptations to the script made by Emma Colliver (2021) 
### Adapted script made available at https://github.com/FrancisCrickInstitute/TRACERxIMCSpilllover 
### Functions below are unchanged from the original script apart from customisation of the isotope list
### Original functions not invoked in the adapted scripts are removed for brevity

### Aim
##Adapts the spillover matrix for use in CellProfiler

### Loading libraries
library(CATALYST)

### Setting up filepaths
fn_sm = '/PATH/TO/SPILLOVER/*sm_adaptive.csv'
fn_imc_metals = '/PATH/TO/metals.csv'
fol_out = dirname(fn_sm) # output folder
prefix_out =  'adaptive_'

#Write the spillover matrix for CellProfiler
sm = read.csv(fn_sm, row.names = 1) 
analysis_channels = read.csv(fn_imc_metals,header = F)
analysis_channels = paste(as.character(analysis_channels$V1), 'Di', sep = '')
custom_isotope_list <- c(CATALYST::isotope_list, list(ArAr=80)) #added Argon
sm_table = CATALYST::adaptSpillmat(input_sm = as.matrix(sm), out_chs = analysis_channels, isotope_list = custom_isotope_list)

#Writes out a 32bit tiff that can be used in CellProfiller together with the "CorrectSpilloverApply" module
imagename = "sm_pixel_adaptive.tiff"
tiff::writeTIFF(sm_table, file.path(fol_out, imagename), bits.per.sample = 32, reduce = T)