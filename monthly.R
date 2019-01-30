####functions to extract monthly covariate values
library(raster)

###function that constructs 

###function to extract monthly covariate values from monthly files
##inputs:
# - covariate name and resolution OR covariate start and end path
# - vector of years and vector of months
# - list of time lags corresponding to each covariate,
#   e.g a first element c(1,2) means time lags of 1 and 2 months for 
#   first covariate
covariate_options <- c("Rain",
                       "LST_day",
                       "LST_night",
                       "TCB",
                       "EVI")
resolutions <- c("1km", "5km", "10km")
file_paths <- list(c("Z:/mastergrids/Other_Global_Covariates/Rainfall/CHIRPS/", "/chirps-v2-0.", ".sum.", ".NN.tif"),
                   c("Z:/mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day/", "/LST_Day_v6.", ".mean.", ".mean.tif"),
                   c("Z:/mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Night/", "/LST_Night_v6.", ".mean.", ".mean.tif"),
                   c("Z:/mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/TCB_v6/", "/TCB_v6.", ".mean.", ".mean.tif"),
                   c("Z:/mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/EVI_v6/", "/EVI_v6.", ".mean.", ".mean.tif")
                   )
names(file_paths) <- covariate_options
monthly_covariates <- function(locations,
                               years,
                               months,
                               covariates = NULL, 
                               resolution = NULL, 
                               startpaths = NULL,
                               endpaths=NULL,
                               covariate_names = NULL,
                               timelags=NULL,
                               crop_extent=NULL,
                               n_blur=NULL){
  #checks
  match.arg(covariates, covariate_options, several.ok=TRUE)
  match.arg(resolution, resolutions)
  #if covariates exist...
  if(!is.null(covariates)){
    #should be same length as timelags
    if(length(covariates)!=length(timelags)) stop("covariates and timelags are different lengths")
    #resolution should exist
    if(is.null(resolution)) stop("resolution not supplied")
    #warn if crop extent is not supplied
    if(is.null(crop_extent)) warning("no extent supplied, this is likely to increase time taken")
  }else{
    #startpaths and endpaths should be supplied
    if(is.null(startpaths)) stop("startpath not supplied")
    if(is.null(endpaths)) stop("endpaths not supplied")
    #ideally names would be supplied
    if(is.null(covariate_names)) warning("no covariates names supplied")
  }
  
  
  N_t <- length(years)
  N_s <- dim(locations)[1]
  N_covs <- length(covariates)
  #if no timelags supplied create list assuming none
  if(is.null(timelags)){
    timelags <- as.list(rep(0, N_covs))
  }
  
  
  if(min(unlist(timelags) < 0)) warnings("negative time lags?")
    
    
  #if covariates are supplied but not names, add names
  if(is.null(covariate_names) & !is.null(covariates)){
    covariate_names <- covariates
  }
  
  out_vals <- list()
  
  #for each covariate
  for(i in 1:N_covs){
    year.month_list <- lapply(timelags[[i]], function(lag) year_months_make(years, months - lag))
    print(year.month_list)
    year.month <- unlist(year.month_list)
    year.month_unique <- unique(year.month)
    N_ymu <- length(year.month_unique)
    if(!is.null(covariate)){
      covariate_paths <- make_monthly_paths(covariates[i], resolution)
    }else{
      covariate_paths <- list(startpaths[[i]], endpaths[[i]])
    }
    
    covariate_value_list_unique <- list()
    for(j in 1:N_ymu){
      print(year.month_unique[j])
      #read raster
      print(paste0(covariate_paths[[1]], year.month_unique[j], covariate_paths[[2]]))
      cov_raster <- raster(paste0(covariate_paths[[1]], year.month_unique[j], covariate_paths[[2]]))
      if(!is.null(crop_extent)) cov_raster <- crop(cov_raster, crop_extent)
      if(!is.null(n_blur)) cov_raster <- focal(cov_raster, w = matrix(1,n_blur,n_blur), mean, na.rm=T)
      cov_vals <- extract(cov_raster, SpatialPoints(locations))
      covariate_value_list_unique[[year.month_unique[j]]] <- cov_vals
    }
    
    for(k in 1:length(timelags[[i]])){
      out_vals[[paste0(covariate_names[i], timelags[[i]][k])]] <- covariate_value_list_unique[match(year.month_list[[k]],
                                                                                          year.month_unique)]
      names(out_vals[[length(out_vals)]]) <- year_months_make(years, months)
    }
    
  }
  
  return(out_vals)
}



###function to construct start and end paths from 

make_monthly_paths <- function(covariate = NULL,
                               resolution = NULL){
  # covariate_options <- c("Rain",
  #                      "LST_day",
  #                      "LST_night",
  #                      "TCB")
  # resolutions <- c("1km", "5km", "10km")
  match.arg(covariate, covariate_options)
  match.arg(resolution, resolutions)
  file_path <- file_paths[[which(names(file_paths) == covariate)]]
  start_path <- paste0(file_path[1], resolution,  "/Monthly", file_path[2])
  end_path <- paste0(file_path[3], resolution, file_path[4])
  return(list(start_path, end_path))
}


###function that takes a year and month and constructs year.month for file path
###allows negative months
library(stringr)
year_month_make <- function(year, month){
  if(month > 12) stop("month should not be greater than 12")
  while(month < 1){
    month <- month + 12
    year <- year - 1
  }
  return(paste(year, str_pad(month, 2, pad="0"), sep = "."))
}

year_months_make <- function(years, months){
  if(length(years)!=length(months)) stop(paste0("different number of years (", length(years), ") and months (", length(months),")"))
  return(sapply(1:length(years), function(i) year_month_make(years[i], months[i])))
}













########test
mad_raster <- raster("Z:/Madagascar_incidence/CovariateData/5km/Monthly_variables/LST_day/LST_Day.2010.03.mean.5km.mean.tif")
# test1 <- monthly_covariates(hf_list[index_use[1:10], ], 
#                             2015:2016,
#                             1:2,
#                             c("Rain", "LST_day"), 
#                             "5km", 
#                             n_blur = 5,
#                             crop_extent=extent(mad_raster))

# test2 <- monthly_covariates(hf_list[index_use[1:10], ],
#                             c(2015, 2015),
#                             1:2,
#                             c("Rain", "LST_day", "EVI", "TCB"),
#                             "5km",
#                             n_blur = 5,
#                             crop_extent=extent(mad_raster),
#                             timelags=list(c(0, 1), 0, 0, 0))

