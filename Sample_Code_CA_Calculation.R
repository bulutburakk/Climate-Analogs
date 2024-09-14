
# This code Calculates Euclidean, Mahalanobis and Wasserstein Distances for a selected coordinates within the study area.
# also plot the results using Sigma Dissimilarity for comparison of the methods
rm(list=ls(all=TRUE))
library(raster)
library(ncdf4)
library(transport)
library(tidyverse)
library(adehabitatLT)   

infold  = # path for the Sample_CA.Rdata
outfold = infold # path for the outfold

load(paste0(infold,"Sample_CA.Rdata")) 
### europe_input_data [6797, 120, 16] # 6797 grid points within extent(c(-10, 50, 30, 70)), 120 years (1981 - 2100), 16 climate Variables 
### data is obtained from gfdl.esm4 model ssp370 scenario - ISIMIP3b

### landseamask_europe   - land sea mask for the study domain
### loc_pix_land_mask_eu - location of each pixel/grid in the input data array (we masked the study area to obtain only grid points over the lands)
### sample_cities        - sample city informations for the CA calculation
### sea_area_shp & study_area_shp - shapefiles for plots

################################################################################

season_name = c("DJF","MAM","JJA","SON")
var_name = c("tas","pr","tasmax","tasmin")

ind_selection_list_run = array(1,dim=c(4,4)) # array of climate variables to be used in the CA analysis. 1 for include, 0 for exclude from the analysis.

rownames(ind_selection_list_run) = var_name
colnames(ind_selection_list_run) = season_name

crop_extent_eur = extent(c(-10, 50, 30, 70)) # extent of the study area - required for preparing result rasters.

# head(europe_input_data[1,,])

################# Functions of the each method to calculate CA #################

CA_Wasser  <- function(pix,p,input_data,parallel){
  
  # pix0 = which(loc_pix_global_land_mask == pix)
  pix0 = which(loc_pix_land_mask_eu == pix)
  ind_length = length(which(ind_selection_list_run==1))
  num_land_pix = dim(input_data)[1]
  
  wasserstein_lapply <- function(var1,var2,p){
    var1pp = pp(var1)
    var2pp = pp(var2)
    wasserstein(var1pp,var2pp,p)
  }
  
  for (ind0 in 1:ind_length){
    
    ## this section standardized all climate variables according to selected reference grid in all periods data 
    ind = which(ind_selection_list_run==1)[ind0]
    
    dummy_HS = input_data[,1:30,ind]
    pix0_HS = dummy_HS[pix0,]
    if(sd(pix0_HS) == 0){pix0_HS[1] = pix0_HS[1] + 0.1 }       ### check sd values - if 0 - then add 0.01 to the first year data
    pix0_F1 = input_data[pix0,31:60,ind]
    if(sd(pix0_F1) == 0){pix0_F1[1] = pix0_F1[1] + 0.1 }       ### check sd values - if 0 - then add 0.01 to the first year data
    pix0_F2 = input_data[pix0,61:90,ind]
    if(sd(pix0_F2) == 0){pix0_F2[1] = pix0_F2[1] + 0.1 }       ### check sd values - if 0 - then add 0.01 to the first year data
    pix0_F3 = input_data[pix0,91:120,ind]
    if(sd(pix0_F3) == 0){pix0_F3[1] = pix0_F3[1] + 0.1 }       ### check sd values - if 0 - then add 0.01 to the first year data
    
    # Standardize all analog pools during the HS period using the reference location's mean and standard deviation for the corresponding period.
    norm_HS = (dummy_HS - mean(pix0_HS)) / sd(pix0_HS)  
    norm_F1 = (dummy_HS - mean(pix0_F1)) / sd(pix0_F1) 
    norm_F2 = (dummy_HS - mean(pix0_F2)) / sd(pix0_F2) 
    norm_F3 = (dummy_HS - mean(pix0_F3)) / sd(pix0_F3)  
    
    # Standardize all reference location climate variables for all periods (by using corresponding period's mean and standard deviation)
    norm_pix0_HS = (pix0_HS - mean(pix0_HS)) / sd(pix0_HS)  
    norm_pix0_F1 = (pix0_F1 - mean(pix0_F1)) / sd(pix0_F1)  
    norm_pix0_F2 = (pix0_F2 - mean(pix0_F2)) / sd(pix0_F2)  
    norm_pix0_F3 = (pix0_F3 - mean(pix0_F3)) / sd(pix0_F3)  
    
    norm_all      =  cbind(norm_HS,norm_F1,norm_F2,norm_F3) 
    norm_pix0_all =  cbind(norm_pix0_HS,norm_pix0_F1,norm_pix0_F2,norm_pix0_F3) 
    
    if(exists("all_input_norm_data")){  
      all_input_norm_data = cbind(all_input_norm_data,norm_all)
      all_input_ref_data  = cbind(all_input_ref_data,norm_pix0_all)}else{
        all_input_norm_data = norm_all
        all_input_ref_data  = norm_pix0_all
      }
    
  } # end ind loop 
  
  if(ind_length > 1){
    
    dim(all_input_norm_data) = c(num_land_pix,30,4,ind_length) # 30 each period's data length & 4 number of periods (HS,F1,F2,F3)
    df <- tibble(Pix_No = 1:num_land_pix)
    df$HS<-lapply(seq_along(df$Pix_No), function(x)  (all_input_norm_data[x,,1,]))
    df$F1<-lapply(seq_along(df$Pix_No), function(x)  (all_input_norm_data[x,,2,]))
    df$F2<-lapply(seq_along(df$Pix_No), function(x)  (all_input_norm_data[x,,3,]))
    df$F3<-lapply(seq_along(df$Pix_No), function(x)  (all_input_norm_data[x,,4,]))
    
    dim(all_input_ref_data) = c(30,4,ind_length)
    
    if(parallel){
      
      distance_HS = mclapply(df$HS,wasserstein_lapply,var2=all_input_ref_data[,1,],p,mc.cores = 16)
      distance_F1 = mclapply(df$F1,wasserstein_lapply,var2=all_input_ref_data[,2,],p,mc.cores = 16)
      distance_F2 = mclapply(df$F2,wasserstein_lapply,var2=all_input_ref_data[,3,],p,mc.cores = 16)
      distance_F3 = mclapply(df$F3,wasserstein_lapply,var2=all_input_ref_data[,4,],p,mc.cores = 16)
      
    }else{
      
      distance_HS = sapply(df$HS,wasserstein_lapply,var2=all_input_ref_data[,1,],p)
      distance_F1 = sapply(df$F1,wasserstein_lapply,var2=all_input_ref_data[,2,],p)
      distance_F2 = sapply(df$F2,wasserstein_lapply,var2=all_input_ref_data[,3,],p)
      distance_F3 = sapply(df$F3,wasserstein_lapply,var2=all_input_ref_data[,4,],p)
      
    }
    
  } else {
    
    dim(all_input_norm_data) = c(global_land_pix_len,30,4)
    df <- tibble(Pix_No = 1:global_land_pix_len)
    df$HS<-lapply(seq_along(df$Pix_No), function(x)  (all_input_norm_data[x,,1]))
    df$F1<-lapply(seq_along(df$Pix_No), function(x)  (all_input_norm_data[x,,2]))
    df$F2<-lapply(seq_along(df$Pix_No), function(x)  (all_input_norm_data[x,,3]))
    df$F3<-lapply(seq_along(df$Pix_No), function(x)  (all_input_norm_data[x,,4]))
    
    dim(all_input_ref_data) = c(30,4)
    
    if(parallel){
      
      distance_HS = mclapply(df$HS,wasserstein1d,b=all_input_ref_data[,1,],p,mc.cores = 16)
      distance_F1 = mclapply(df$F1,wasserstein1d,b=all_input_ref_data[,2,],p,mc.cores = 16)
      distance_F2 = mclapply(df$F2,wasserstein1d,b=all_input_ref_data[,3,],p,mc.cores = 16)
      distance_F3 = mclapply(df$F3,wasserstein1d,b=all_input_ref_data[,4,],p,mc.cores = 16)
      
    }else{
      
      distance_HS = sapply(df$HS,wasserstein1d,b=all_input_ref_data[,1,],p)
      distance_F1 = sapply(df$F1,wasserstein1d,b=all_input_ref_data[,2,],p)
      distance_F2 = sapply(df$F2,wasserstein1d,b=all_input_ref_data[,3,],p)
      distance_F3 = sapply(df$F3,wasserstein1d,b=all_input_ref_data[,4,],p)
      
    }
    
  }
  
  # this section creates raster for each period 
  if(1){
    
    nrows = nrow(landseamask_europe)
    ncols = ncol(landseamask_europe)
    
    new_raster_ar = array(NA, dim=c(nrows,ncols))
    new_raster_ar[loc_pix_land_mask_eu] = unlist(distance_F1)
    new_raster_mat  <- matrix(new_raster_ar,nrows,ncols,byrow=TRUE)
    distance_ras_F1 <- raster(new_raster_mat,
                              xmn=crop_extent_eur[1], xmx=crop_extent_eur[2], ymn=crop_extent_eur[3], ymx=crop_extent_eur[4],
                              crs="+proj=longlat +datum=WGS84 +no_defs")
    names(distance_ras_F1) = paste0("Wasser_F1")
    
    new_raster_ar = array(NA, dim=c(nrows,ncols))
    new_raster_ar[loc_pix_land_mask_eu] = unlist(distance_F2)
    new_raster_mat  <- matrix(new_raster_ar,nrows,ncols,byrow=TRUE)
    distance_ras_F2 <- raster(new_raster_mat,
                              xmn=crop_extent_eur[1], xmx=crop_extent_eur[2], ymn=crop_extent_eur[3], ymx=crop_extent_eur[4],
                              crs="+proj=longlat +datum=WGS84 +no_defs") 
    names(distance_ras_F2) = paste0("Wasser_F2")
    
    new_raster_ar = array(NA, dim=c(nrows,ncols))
    new_raster_ar[loc_pix_land_mask_eu] = unlist(distance_F3)
    new_raster_mat  <- matrix(new_raster_ar,nrows,ncols,byrow=TRUE)
    distance_ras_F3 <- raster(new_raster_mat,
                              xmn=crop_extent_eur[1], xmx=crop_extent_eur[2], ymn=crop_extent_eur[3], ymx=crop_extent_eur[4],
                              crs="+proj=longlat +datum=WGS84 +no_defs")
    names(distance_ras_F3) = paste0("Wasser_F3")
    
    new_raster_ar = array(NA, dim=c(nrows,ncols))
    new_raster_ar[loc_pix_land_mask_eu] = unlist(distance_HS)
    new_raster_mat  <- matrix(new_raster_ar,nrows,ncols,byrow=TRUE)
    distance_ras_HS <- raster(new_raster_mat,
                              xmn=crop_extent_eur[1], xmx=crop_extent_eur[2], ymn=crop_extent_eur[3], ymx=crop_extent_eur[4],
                              crs="+proj=longlat +datum=WGS84 +no_defs")
    names(distance_ras_HS) = paste0("Wasser_HS")
    
    if(exists("ras_CA_Wasser_Clim")){  
      ras_CA_Wasser_Clim = stack(ras_CA_Wasser_Clim,distance_ras_HS,distance_ras_F1,distance_ras_F2,distance_ras_F3)}else{
        ras_CA_Wasser_Clim = stack(distance_ras_HS,distance_ras_F1,distance_ras_F2,distance_ras_F3)
      }
    
  }
  
  return(ras_CA_Wasser_Clim)
  
}

CA_Euclid  <- function(pix,input_data){
  
  # pix0 = which(loc_pix_global_land_mask == pix)
  pix0 = which(loc_pix_land_mask_eu == pix)
  ind_length = length(which(ind_selection_list_run==1))
  num_land_pix = dim(input_data)[1]
  
  for (ind0 in 1:ind_length){
    
    ind = which(ind_selection_list_run==1)[ind0]
    dummy_HS = input_data[,1:30,ind]
    
    each_pix_mean_HS = rowMeans(dummy_HS)
    sd_each_pix_mean_HS = sd(each_pix_mean_HS)
    
    pix0_HS_mean = mean(dummy_HS[pix0,])
    pix0_F1_mean = mean(input_data[pix0,31:60,ind])
    pix0_F2_mean = mean(input_data[pix0,61:90,ind])
    pix0_F3_mean = mean(input_data[pix0,91:120,ind])
    
    pix0_HS_sd = sd(dummy_HS[pix0,])
    pix0_F1_sd = sd(input_data[pix0,31:60,ind])
    pix0_F2_sd = sd(input_data[pix0,61:90,ind])
    pix0_F3_sd = sd(input_data[pix0,91:120,ind])
    
    HS = ((each_pix_mean_HS-pix0_HS_mean)^2 / (pix0_HS_sd)^2) 
    F1 = ((each_pix_mean_HS-pix0_F1_mean)^2 / (pix0_F1_sd)^2) 
    F2 = ((each_pix_mean_HS-pix0_F2_mean)^2 / (pix0_F2_sd)^2) 
    F3 = ((each_pix_mean_HS-pix0_F3_mean)^2 / (pix0_F3_sd)^2) 
    
    sqrt_dist_all =  cbind(HS,F1,F2,F3) #* ind_weighting_list[ind] 
    rm(HS,F1,F2,F3)
    
    if(exists("all_sqrt_dist_all")){  
      all_sqrt_dist_all = cbind(all_sqrt_dist_all, sqrt_dist_all)}else{
        all_sqrt_dist_all = sqrt_dist_all
      }
    
  } # end ind
  
  if(ind_length > 1){
    
    dim(all_sqrt_dist_all) = c(num_land_pix,4,ind_length) # 30 each period's data length & 4 number of periods (HS,F1,F2,F3)
    
    distance_HS = sqrt(rowSums(all_sqrt_dist_all[,1,]))
    distance_F1 = sqrt(rowSums(all_sqrt_dist_all[,2,]))
    distance_F2 = sqrt(rowSums(all_sqrt_dist_all[,3,]))
    distance_F3 = sqrt(rowSums(all_sqrt_dist_all[,4,]))
    
    rm(all_sqrt_dist_all)
    
  } else {
    
    dim(all_sqrt_dist_all) = c(num_land_pix,4) # 30 each period's data length & 4 number of periods (HS,F1,F2,F3)
    
    distance_HS = sqrt(sum(all_sqrt_dist_all[,1]))
    distance_F1 = sqrt(sum(all_sqrt_dist_all[,2]))
    distance_F2 = sqrt(sum(all_sqrt_dist_all[,3]))
    distance_F3 = sqrt(sum(all_sqrt_dist_all[,4]))
    
    rm(all_sqrt_dist_all)
    
  }
  
  # this section creates raster for each period 
  if(1){
    
    nrows = nrow(landseamask_europe)
    ncols = ncol(landseamask_europe)
    
    new_raster_ar = array(NA, dim=c(nrows,ncols))
    new_raster_ar[loc_pix_land_mask_eu] = unlist(distance_F1)
    new_raster_mat <- matrix(new_raster_ar,nrows,ncols,byrow=TRUE)
    distance_ras_F1 <- raster(new_raster_mat, xmn=crop_extent_eur[1], xmx=crop_extent_eur[2], ymn=crop_extent_eur[3], ymx=crop_extent_eur[4], crs="+proj=longlat +datum=WGS84 +no_defs") # , transpose=FALSE
    names(distance_ras_F1) = paste0("Euclid_F1")
    
    new_raster_ar = array(NA, dim=c(nrows,ncols))
    new_raster_ar[loc_pix_land_mask_eu] = unlist(distance_F2)
    new_raster_mat <- matrix(new_raster_ar,nrows,ncols,byrow=TRUE)
    distance_ras_F2 <- raster(new_raster_mat, xmn=crop_extent_eur[1], xmx=crop_extent_eur[2], ymn=crop_extent_eur[3], ymx=crop_extent_eur[4], crs="+proj=longlat +datum=WGS84 +no_defs") # , transpose=FALSE
    names(distance_ras_F2) = paste0("Euclid_F2")
    
    new_raster_ar = array(NA, dim=c(nrows,ncols))
    new_raster_ar[loc_pix_land_mask_eu] = unlist(distance_F3)
    new_raster_mat <- matrix(new_raster_ar,nrows,ncols,byrow=TRUE)
    distance_ras_F3 <- raster(new_raster_mat, xmn=crop_extent_eur[1], xmx=crop_extent_eur[2], ymn=crop_extent_eur[3], ymx=crop_extent_eur[4], crs="+proj=longlat +datum=WGS84 +no_defs") # , transpose=FALSE
    names(distance_ras_F3) = paste0("Euclid_F3")
    
    new_raster_ar = array(NA, dim=c(nrows,ncols))
    new_raster_ar[loc_pix_land_mask_eu] = unlist(distance_HS)
    new_raster_mat <- matrix(new_raster_ar,nrows,ncols,byrow=TRUE)
    distance_ras_HS <- raster(new_raster_mat, xmn=crop_extent_eur[1], xmx=crop_extent_eur[2], ymn=crop_extent_eur[3], ymx=crop_extent_eur[4], crs="+proj=longlat +datum=WGS84 +no_defs") # , transpose=FALSE
    names(distance_ras_HS) = paste0("Euclid_HS")
    
    if(exists("ras_CA_Euclid_Clim")){  
      ras_CA_Euclid_Clim = stack(ras_CA_Euclid_Clim,distance_ras_HS,distance_ras_F1,distance_ras_F2,distance_ras_F3)}else{
        ras_CA_Euclid_Clim = stack(distance_ras_HS,distance_ras_F1,distance_ras_F2,distance_ras_F3)
      }
    rm(distance_ras_HS,distance_ras_F1,distance_ras_F2,distance_ras_F3)
    gc(verbose = FALSE)
    
  }
  
  return(ras_CA_Euclid_Clim)
  
}

CA_Mahalan <- function(pix,input_data){
  
  # pix0 = which(loc_pix_global_land_mask == pix)
  pix0 = which(loc_pix_land_mask_eu == pix)
  ind_length = length(which(ind_selection_list_run==1))
  num_land_pix = dim(input_data)[1]
  
  MD_func = function(data,ICV_data){
    
    trunc.SDs <- 0.1 # truncation 
    
    mean_ref = colMeans(ICV_data)
    sd_ref = apply(ICV_data,2,sd)
    
    var_no = dim(data)[2]
    
    res = data * NA
    std_ref = ICV_data * NA
    
    for(i in 1:var_no){
      res[,i] = (data[,i] - mean_ref[i]) / sd_ref[i]
      std_ref[,i] = (ICV_data[,i] - mean_ref[i]) / sd_ref[i]
    }
    
    PCA <- prcomp(std_ref)  
    PCs <- max(which(unlist(summary(PCA)[1])>trunc.SDs)) 
    
    X <- predict(PCA,std_ref)   # project the analog pool onto the PCs
    Yj<- predict(PCA,res)       # project the projected future conditions onto the PCs
    
    res0 = sqrt(mahalanobis(Yj[,1:PCs],center=FALSE,cov(X[,1:PCs],use='na.or.complete')))
    
    return(res0)
  }
  
  ind = which(ind_selection_list_run==1) # ind_selection_list_run is for runing the code since loaded file also has ind_selection_list !
  
  pix0_HS_all_inds = input_data[pix0,1:30,ind]
  pix0_EF_all_inds = input_data[pix0,31:60,ind]
  pix0_MF_all_inds = input_data[pix0,61:90,ind]
  pix0_FF_all_inds = input_data[pix0,91:120,ind]
  
  dummy_HS0 = input_data[,1:30,ind]
  dummy_HS  = aperm(dummy_HS0,c(1,3,2))
  
  input_data_HS = rowMeans(dummy_HS,dim=2)
  
  distance_HS = MD_func(input_data_HS,pix0_HS_all_inds)
  distance_F1 = MD_func(input_data_HS,pix0_EF_all_inds)
  distance_F2 = MD_func(input_data_HS,pix0_MF_all_inds)
  distance_F3 = MD_func(input_data_HS,pix0_FF_all_inds)
  
  # this section creates raster for each period 
  if(1){
    
    nrows = nrow(landseamask_europe)
    ncols = ncol(landseamask_europe)
    
    new_raster_ar = array(NA, dim=c(nrows,ncols))
    new_raster_ar[loc_pix_land_mask_eu] = unlist(distance_F1)
    new_raster_mat <- matrix(new_raster_ar,nrows,ncols,byrow=TRUE)
    distance_ras_F1 <- raster(new_raster_mat, xmn=crop_extent_eur[1], xmx=crop_extent_eur[2], ymn=crop_extent_eur[3], ymx=crop_extent_eur[4], crs="+proj=longlat +datum=WGS84 +no_defs") # , transpose=FALSE
    names(distance_ras_F1) = paste0("Mahalan_F1")
    
    new_raster_ar = array(NA, dim=c(nrows,ncols))
    new_raster_ar[loc_pix_land_mask_eu] = unlist(distance_F2)
    new_raster_mat <- matrix(new_raster_ar,nrows,ncols,byrow=TRUE)
    distance_ras_F2 <- raster(new_raster_mat, xmn=crop_extent_eur[1], xmx=crop_extent_eur[2], ymn=crop_extent_eur[3], ymx=crop_extent_eur[4], crs="+proj=longlat +datum=WGS84 +no_defs") # , transpose=FALSE
    names(distance_ras_F2) = paste0("Mahalan_F2")
    
    new_raster_ar = array(NA, dim=c(nrows,ncols))
    new_raster_ar[loc_pix_land_mask_eu] = unlist(distance_F3)
    new_raster_mat <- matrix(new_raster_ar,nrows,ncols,byrow=TRUE)
    distance_ras_F3 <- raster(new_raster_mat, xmn=crop_extent_eur[1], xmx=crop_extent_eur[2], ymn=crop_extent_eur[3], ymx=crop_extent_eur[4], crs="+proj=longlat +datum=WGS84 +no_defs") # , transpose=FALSE
    names(distance_ras_F3) = paste0("Mahalan_F3")
    
    new_raster_ar = array(NA, dim=c(nrows,ncols))
    new_raster_ar[loc_pix_land_mask_eu] = unlist(distance_HS)
    new_raster_mat <- matrix(new_raster_ar,nrows,ncols,byrow=TRUE)
    distance_ras_HS <- raster(new_raster_mat, xmn=crop_extent_eur[1], xmx=crop_extent_eur[2], ymn=crop_extent_eur[3], ymx=crop_extent_eur[4], crs="+proj=longlat +datum=WGS84 +no_defs") # , transpose=FALSE
    names(distance_ras_HS) = paste0("Mahalan_HS")
    
    if(exists("ras_CA_Mahalan_Clim")){  
      ras_CA_Mahalan_Clim = stack(ras_CA_Mahalan_Clim,distance_ras_HS,distance_ras_F1,distance_ras_F2,distance_ras_F3)}else{
        ras_CA_Mahalan_Clim = stack(distance_ras_HS,distance_ras_F1,distance_ras_F2,distance_ras_F3)
      }
    rm(distance_ras_HS,distance_ras_F1,distance_ras_F2,distance_ras_F3)
    gc(verbose = FALSE)
    
  }
  
  return(ras_CA_Mahalan_Clim)
  
}


dummy_coordinates = c(25.279652,54.687157) # The latitude of Vilnius, Lithuania is 54.687157, and the longitude is 25.279652.
dummy_cell        = cellFromXY(landseamask_europe,dummy_coordinates)
# pix = dummy_cell

pix = city_cells[sel_city]  # Grid/Pixel Number of the focal location  can be selected from sample cities or user can be find the pixel number using the coordinates within study extent

eur_ras_Euclid  = CA_Euclid(pix,europe_input_data)
eur_ras_Mahalan = CA_Mahalan(pix,europe_input_data)
eur_ras_Wasser  = CA_Wasser(pix,2,europe_input_data,parallel=FALSE)
plot(eur_ras_Wasser,zlim=range(values(eur_ras_Wasser),na.rm=T))  # plots the pure WD Results for the selected location

ras_names    = c("eur_ras_Euclid","eur_ras_Mahalan","eur_ras_Wasser")
# the loop calculates Sigma Dissimilarity values for all 3 methods
for(method in 1:3){
  
  if(method == 1 & exists("all_result_rasters")){rm(all_result_rasters)}   # check and deletes previous analysis results
  
  selected_method_raster = eval(parse(text = ras_names[method])) # rasters of the selected method
  all_distance_values    = c(values(selected_method_raster))     # obtains all distance values from the rasters (all 4 periods Hs to FF)
  
  all_SD_values = all_distance_values * NA                       # creates an empty array for Sigma Dissimilarity (SD) values
  loc_na = which(is.na(all_distance_values))                     # locates the NA values in the distance raster array (grids over sea)
  data = as.vector(all_distance_values[-loc_na])                 # selects grid points without NA data for the SD calculation.
  
  NN.chi_values <- pchi(data, df = 16)                           # pchi gives the distribution function (df value is 16 if all avaliable climate variables in use)
  NN.sigma_values  <- qchi(NN.chi_values,1)                      # qchi gives the quantile function
  all_SD_values[-loc_na] = NN.sigma_values
  
  n_cell   = ncell(selected_method_raster)
  end_per  = seq(n_cell,n_cell*4,n_cell)
  strt_per = c(1,(end_per + 1 )[-4])
  
  for(periods in 1:4){
    
    dummy_ras = selected_method_raster[[periods]] * NA                          # creates an empty raster for Sigma Dissimilarity (SD) values
    values(dummy_ras) = all_SD_values[strt_per[periods]:end_per[periods]]       # fill raster with corresponding period's SD values
    
    # below 3 lines converts NaN. Infinite and more than 6 sigma values into 7
    dummy_ras[dummy_ras>6] = 7
    dummy_ras[which(values(dummy_ras)=="NaN")] = 7
    dummy_ras[which(is.infinite(values(dummy_ras)))] = 7
    
    
    if(exists("all_result_rasters")){all_result_rasters = stack(all_result_rasters,dummy_ras)}else{all_result_rasters = dummy_ras}} # stacks all 12 rasters (4 Periods x 3 Methods)
} # prep of plot rasters using Sigma Dissimilarity 

################################################################################
### CA Maps Plot

method_names = c("Euclidean","Mahalanobis","Wasserstein")
pers         = c("HS (1981-2010)", "EF (2011-2040)", "MF (2041-2070)", "FF (2071-2100)")
tick_names   = c( expression(paste("0",sigma)),NA,expression(paste("2",sigma)),NA,expression(paste("4",sigma)),NA,expression(paste("6",sigma)) )

col_palraw  <- c('#313695','#a50026','#d73027','#f46d43','#fdae61','#fee090','#d9d9d9')
breakpoints <- c(seq(0,7,1))
ColScheme   <- colorRampPalette(col_palraw)(length(breakpoints)-1) 

if(exists("sel_city")){fig_name = city_names[sel_city]}else{fig_name = pix}     # obtains figure name is pixel selected from the sample cities

png(filename = paste0(outfold,"CA_",fig_name,"_All_Methods.png",sep=""),width = 29, height = 21, units = "cm", bg = "white",  res = 600)
par(mar=c(1.25,0,0,0),pty = "m",oma=c(1.75,1.75,1.75,5.50),mfrow=c(3,4)) 

for(method in 1:3){
  
  for(periods in 1:4){
    
    i = periods + (method-1)*4
    
    plot_raster = all_result_rasters[[i]]
    
    plot(dummy_ras, col=ColScheme, breaks=breakpoints,axes=FALSE,legend=FALSE,box=FALSE)
    plot(study_area_shp,add=TRUE,lwd=.50)
    plot(sea_area_shp,add=TRUE,col="#9ecae1",lwd=.50)
    
    if( method  == 3 ) {axis(1,at=seq(-10,50,10),line=-0.050,cex.axis=1.00,padj=-0.80)}
    if( periods == 4 ) {axis(4,at=seq(30,70,10),line=-0.20,cex.axis=1.00,padj=-1)}
    if( periods == 1 ) {mtext(method_names[method],side=2,line=0.5,font=2,cex=0.85)}  
    if( method  == 1 ) {mtext(pers[periods],side=3,line=0,font=2,cex=0.85)}
    
    SD_categories = cut(values(plot_raster),seq(0,7,1),include.lowest = TRUE)
    
    legend(x= -10 ,y= 70 ,legend = table(SD_categories), pch = 22, pt.bg = col_palraw, pt.cex = 1.5, horiz = FALSE, bty = "n",ncol=2, 
           x.intersp=0.5,y.intersp= 0.75,text.width = 4.5)
    
  } # end periods HS to FF
  
}# end method 

par(mfrow=c(1,1),new=FALSE, oma=c(0,29.5,0,0))
plot(plot_raster,legend.only=TRUE, col=ColScheme, breaks=breakpoints,legend.width=2, legend.shrink=.4, horizontal = FALSE,
     axis.args=list(cex.axis=0.85,line = -0.75,font=2,tick=FALSE,at=0:6,labels = tick_names),
     legend.args=list(text="Dissimilarity", side=2, font=2, line=.25, cex=0.85))

dev.off()

################################################################################