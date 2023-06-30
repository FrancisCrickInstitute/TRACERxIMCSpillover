################## SPILLOVER IMC HELPER FUNCTIONS ########################

### Original script authored by Vito Zanotelli et al. 
### Original script available at https://github.com/BodenmillerGroup/cyTOFcompensation/blob/master/scripts/spillover_imc_helpers.R
### Adaptations to the script made by Emma Colliver (2021) 
### Adapted script made available at https://github.com/FrancisCrickInstitute/TRACERxIMCSpilllover 
### Any adaptations to individual functions or newly written functions are explicitly highlighted below
### Original functions not invoked in the adapted scripts are removed for brevity

### Loading libraries
library(CATALYST)
library(data.table)
library(flowCore)
library(dplyr)
library(dtplyr)
library(stringi)

### Helper to load single stain .txt files 
load_ss_fol <- function(fol_ss){
  # helper function to load all .txt files from a folder
  fns_txt <- list.files(fol_ss,pattern = '*.[0-9]+.txt$') 
  imgs.ss <- lapply(fns_txt, function(x){
    fread(file.path(fol_ss, x))})
  names(imgs.ss) <- fns_txt
  return(imgs.ss)
}

### Extracting correct metal names from .txt file name 
get_metals_from_txtname <- function(nam){
  nam =  gsub('.*\\(', '', nam)
  nam = gsub('\\)', '',nam)
  return(nam)
}

### Extracting the correct metal names from the .txt files 
fixnames <- function(imgdat){
  imgdat = copy(imgdat)
  dat =imgdat
  colnames(dat) = sapply(colnames(dat), function(x) get_metals_from_txtname(x))
  return(dat)
}

### Creating a data file from the list of files
imglist2dat <- function(datlist){
  imgdat <- rbindlist(datlist, fill=T, idcol = 'file')
  return(imgdat)
}

### Calculating file medians 
calc_file_medians <- function(dat){
  # calculates medians per file
  tdat = dat %>%  
    dplyr::select(-c(Start_push, End_push, Pushes_duration,   X , Y  ,  Z)) %>%
    melt.data.table(id.vars = c('metal', 'mass','file')) %>%
    do(data.table(.)[, list(med=median(value)), by=.(variable, metal, mass, file)]) 
  return(tdat)
  
}

### Restricting analysis to the brightest 50% of pixels to minimise subsequent need for 
### pixel binning and associated introduction of noise - NEW FUNCTION
retain_brightest_signal <- function(dat) {
  retained_pixels <- dat
  retained_pixels <- retained_pixels[-c(1:nrow(retained_pixels)),] 
  for (i in unique(dat$metal)) {
    metal_under_study <- i
    ### specific exceptions captured for inconsistent isotope labelling in experimental data 
    if ((metal_under_study != "Nd144") & (metal_under_study != "Er160")) {
      lower_quantile <- quantile(dat[which(dat$metal == metal_under_study),][[paste0(metal_under_study,"Di")]],0.5)
      retained_pixels_metal <- dat[which((dat$metal == metal_under_study) & (dat[[paste0(metal_under_study,"Di")]] >= lower_quantile)),]
    } else if (metal_under_study == "Nd144")  {
      lower_quantile <- quantile(dat[which(dat$metal == metal_under_study),][["Sm144Di"]],0.5)
      retained_pixels_metal <- dat[which((dat$metal == metal_under_study) & (dat[["Sm144Di"]] >= lower_quantile)),]
    } else if (metal_under_study == "Er160") {
      lower_quantile <- quantile(dat[which(dat$metal == metal_under_study),][["Gd160Di"]],0.5)
      retained_pixels_metal <- dat[which((dat$metal == metal_under_study) & (dat[["Gd160Di"]] >= lower_quantile)),]
    }
    retained_pixels <- rbind(retained_pixels,retained_pixels_metal)
  }
  retained_pixels
}

############ Helpers to aggregate consecutive pixels 
get_consecutive_bin <- function(nel, nbin){
  # gets consecutive pixels
  idx = rep(1:ceiling(nel/nbin), each=nbin)
  return(idx[1:nel])
}

aggregate_pixels <- function(dat, n){
  # sums over n consecutive pixels
  tdat = dat[, rowsum(.SD, get_consecutive_bin(.N, n)) ,by=.(file, mass, metal)]
  return(tdat)
}

aggregate_pixels_by_marker <- function(dat, m, n){ 
  # sums over n consecutive pixels
  dat_reduced <- dat[[1]][which(dat[[1]][,"metal"]==m[1]),]
  tdat = dat_reduced[, rowsum(.SD, get_consecutive_bin(.N, n[1])) ,by=.(file, mass, metal)]
  for (i in 2:length(m)) {
    dat_reduced <- dat[[1]][which(dat[[1]][,"metal"]==m[i]),]
    dat_agg = dat_reduced[, rowsum(.SD, get_consecutive_bin(.N, n[i])) ,by=.(file, mass, metal)]
    tdat = rbind.data.frame(tdat,dat_agg)
  }
  return(tdat)
}

############ Performing compensation 
filter_rare_bc <- function(re, minevents){
  # allows filtering out of rare events
  stats = table(re@bc_ids)
  nonfreq = names(stats)[stats <minevents]
  re@bc_ids[re@bc_ids %in% nonfreq] = '0'
  return(re)
}

ensure_correct_bc <- function(re, mass){
  # enforces correct barcodes according to mass/metal identified from the file name
  re@bc_ids[re@bc_ids != as.character(mass)] = '0'
  return(re)
}

re_from_dat <- function(dat, ss_ms, minevents=10, correct_bc=NULL){
  # Debarcode the data, enforce some minimal quality
  # dat: data
  # bc_ms: list of masses used for the single stains
  # minevents: minimum number of events a metal needs to have to be considered
  # correct_bc: list of ground truth barcodes, e.g. from the filenames
  ff = dat %>%
    dplyr::select(-c(file, mass, metal)) %>%
    as.matrix.data.frame() %>%
    flowFrame()
  
  re <- CATALYST::assignPrelim(x=ff,y= ss_ms)
  re <- estCutoffs(re)
  re <- applyCutoffs(re)
  
  # filter for conditions with less then minevents
  if (!is.null(correct_bc) && correct_bc) {
    re = ensure_correct_bc(re, dat[, mass])
  }
  
  re = filter_rare_bc(re, minevents)
  return(re)
}

sm_from_dat <- function(dat, ss_ms, minevents=10, remove_incorrect_bc=T, ...){
  # Helper function to debarcode and calculate spillover
  # ss_ms: masses used for single stains
  # minevents: minimal number of events to consider spillover estimation
  # remove_incorrect_bc; enforce correct barcode by using
  re <- re_from_dat(dat, ss_ms, minevents,remove_incorrect_bc)
  sm <- computeSpillmat(re, ...)
  return(sm)
}