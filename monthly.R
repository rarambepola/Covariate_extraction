#get current working directory
current_working_directory <- getwd()

#set working directory
setwd(paste0(Sys.getenv("HOME"), "/Madagascar_v2"))

#get z drive path
load("zdrive_path.RData")

#reset to previous working directory
setwd(current_working_directory)

####functions to extract monthly covariate values
library(raster)
library(foreach)
library(doParallel)


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
                       "EVI",
                       "TSI_pf",
                       "TSI_pv")
resolutions <- c("1km", "5km", "10km")
modis_resolution_summaries <- c(".Data.tif", ".mean.tif", NA)
tsi_resolution_summaries <- c(NA, ".Data.tif", NA)
file_paths <- list(
  #rain
  list(paste0(zdrive_path, "mastergrids/Other_Global_Covariates/Rainfall/CHIRPS/"), "/chirps-v2-0.", ".sum.", 
       c(".NN.tif", ".NN.tif", ".Data.tif")),
  #LST_day
  list(paste0(zdrive_path, "mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day/"), "/LST_Day_v6.", ".mean.", 
       modis_resolution_summaries),
  #LST_night
  list(paste0(zdrive_path, "mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Night/"), "/LST_Night_v6.", ".mean.",
       modis_resolution_summaries),
  #TCB
  list(paste0(zdrive_path, "mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/TCB_v6/"), "/TCB_v6.", ".mean.", 
       modis_resolution_summaries),
  #EVI
  list(paste0(zdrive_path, "mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/EVI_v6/"), "/EVI_v6.", ".mean.",
       modis_resolution_summaries),
  #TSI_pf
  list(paste0(zdrive_path, "mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/"), "/TSI-Martens2-Pf.", ".Data.",
       tsi_resolution_summaries),
  #TSI_pv
  list(paste0(zdrive_path, "mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pv_Dynamic/"), "/TSI-Martens2-Pv.", ".Data.",
       tsi_resolution_summaries)
  )

for(i in 1:length(file_paths)){
  names(file_paths[[i]][[4]]) <- resolutions
}


names(file_paths) <- covariate_options

####monthly extraction function####
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
                               n_blur=NULL,
                               reverse_order=FALSE,
                               run_parallel=FALSE,
                               n_cores=NULL){
  #check if needs to be run in parallel
  if(run_parallel){
    library(doParallel)
    if(is.null(n_cores)) stop("number of cores not supplied")
  }
  
  #checks
  match.arg(covariates, covariate_options, several.ok=TRUE)
  match.arg(resolution, resolutions)
  #if covariates exist...
  if(!is.null(covariates)){
    #should be same length as timelags
    if(!is.null(timelags)){
      if(length(covariates)!=length(timelags)) stop("covariates and timelags are different lengths")
    }
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
  
  #if run in parallel, create cluster
  if(run_parallel){
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
  }
  `%do_seq_par%` <- ifelse(run_parallel, `%dopar%`, `%do%`)
  
  #for each covariate
  for(i in 1:N_covs){
    print(covariates[i])
    #print(i)
    year.month_list <- lapply(timelags[[i]], function(lag) year_months_make(years, months - lag))
    #print(year.month_list)
    year.month <- unlist(year.month_list)
    year.month_unique <- unique(year.month)
    N_ymu <- length(year.month_unique)
    if(!is.null(covariates)){
      covariate_paths <- make_monthly_paths(covariates[i], resolution)
      #covariate_paths <- make_monthly_paths("Rain", resolution)
      #print(covariate_paths)
    }else{
      covariate_paths <- list(startpaths[[i]], endpaths[[i]])
    }
    
    covariate_value_list_unique <- list()
    
    covariate_value_list_unique <- foreach(j=1:N_ymu,
                                           .packages=c("raster")) %do_seq_par% {
      cov_raster <- raster(paste0(covariate_paths[[1]], year.month_unique[j], covariate_paths[[2]]))
      if(!is.null(crop_extent)) cov_raster <- crop(cov_raster, crop_extent)
      if(!is.null(n_blur)) cov_raster <- focal(cov_raster, w = matrix(1,n_blur,n_blur), mean, na.rm=T)
      cov_vals <- extract(cov_raster, SpatialPoints(locations))
      cov_vals
    }

    names(covariate_value_list_unique) <- year.month_unique
    # for(j in 1:N_ymu){
    #   print(year.month_unique[j])
    #   #read raster
    #   print("covariate_paths")
    #   print(covariate_paths)
    #   print(paste0(covariate_paths[[1]], year.month_unique[j], covariate_paths[[2]]))
    #   cov_raster <- raster(paste0(covariate_paths[[1]], year.month_unique[j], covariate_paths[[2]]))
    #   if(!is.null(crop_extent)) cov_raster <- crop(cov_raster, crop_extent)
    #   #if(!is.null(n_blur)) cov_raster <- focal(cov_raster, w = matrix(1/(n_blur*n_blur),n_blur,n_blur))
    #   if(!is.null(n_blur)) cov_raster <- focal(cov_raster, w = matrix(1,n_blur,n_blur), mean, na.rm=T)
    #   cov_vals <- extract(cov_raster, SpatialPoints(locations))
    #   covariate_value_list_unique[[year.month_unique[j]]] <- cov_vals
    # }
    
    for(k in 1:length(timelags[[i]])){
      out_vals[[paste0(covariate_names[i], timelags[[i]][k])]] <- covariate_value_list_unique[match(year.month_list[[k]],
                                                                                          year.month_unique)]
      names(out_vals[[length(out_vals)]]) <- year_months_make(years, months)
    }
  }
  
  if(run_parallel) stopCluster(cl)
  
  #if reverse order, return a list indexed first by time and then by covariate
  if(reverse_order){
    out_vals_reverse <- list()
    for(i in 1:length(out_vals[[1]])){
      time_names <- names(out_vals[[1]])
      for(j in 1:length(out_vals)){
        cov_names <- names(out_vals)
        out_vals_reverse[[time_names[i]]][[cov_names[j]]] <- out_vals[[j]][[i]]
      }
    }
    return(out_vals_reverse)
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
  start_path <- paste0(file_path[[1]], resolution,  "/Monthly", file_path[[2]])
  end_path <- paste0(file_path[[3]], resolution, file_path[[4]][which(names(file_path[[4]]) == resolution)])
  

  #print(file_path[[4]][which(names(file_path[[4]]) == resolution)])
  #check if resolution is available
  if(is.na(file_path[[4]][which(names(file_path[[4]]) == resolution)])){
    stop(paste0("resolution ", resolution, " not available for covariate ", covariate))
  }
  
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
  if(length(years)!=length(months)) stop(paste0("different number of years (", length(years), 
                                                ") and months (", length(months),")"))
  return(sapply(1:length(years), function(i) year_month_make(years[i], months[i])))
}








####individual monthly extraction function####
##function that splits vector into list of equal length parts
vector_split <- function(full_vector, N_parts=NULL, N_length=NULL){

  N_vec <- length(full_vector)
  #checks
  if(is.null(N_parts) & is.null(N_length)) stop("need either number of parts of length of each part")
  if(!is.null(N_parts) & !is.null(N_length)) stop("both number and length of part given, only need one")
  #check that whichever one is given divides the size of the vector
  if(!is.null(N_parts)){
    if(N_vec %% N_parts) stop("number does not divide length of vector")
    N_length <- N_vec / N_parts
  }
  if(!is.null(N_length)){
    if(N_vec %% N_length) stop("length does not divide length of vector")
    N_parts <- N_vec / N_length
  }
  
  out_list <- list()
  for(i in 1:N_parts){
    out_list[[i]] <- full_vector[1:N_length + ((i-1) * N_length)]
  }
  return(out_list)
}


##function that extracts at each point at a given yearmonth
##instead of every location at every month
monthly_covariates_individual <- function(locations,
                                          years,
                                          months,
                                          covariates = NULL, 
                                          resolution = NULL, 
                                          startpaths = NULL,
                                          endpaths=NULL,
                                          covariate_names = NULL,
                                          timelags=NULL,
                                          crop_extent=NULL,
                                          n_blur=NULL,
                                          reverse_order=FALSE,
                                          run_parallel=FALSE,
                                          n_cores=NULL){
  #check if needs to be run in parallel
  if(run_parallel){
    library(doParallel)
    if(is.null(n_cores)) stop("number of cores not supplied")
  }
  
  #checks
  match.arg(covariates, covariate_options, several.ok=TRUE)
  match.arg(resolution, resolutions)
  #if covariates exist...
  if(!is.null(covariates)){
    #should be same length as timelags
    if(!is.null(timelags)){
      if(length(covariates)!=length(timelags)) stop("covariates and timelags are different lengths")
    }
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
  N_months <- length(months)
  
  #months should be same length as years and locations
  if(N_t != N_s) stop("different number of locations and years")
  if(N_t != N_months) stop("different number of years and months")
  
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
  
  #if run in parallel, create cluster
  if(run_parallel){
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
  }
  `%do_seq_par%` <- ifelse(run_parallel, `%dopar%`, `%do%`)
  
  out_vals <- list()
  
  #for each covariate
  for(i in 1:N_covs){
    N_timelags <- length(timelags[[i]])
    print(covariates[i])
    year.month_list <- lapply(timelags[[i]], function(lag) year_months_make(years, months - lag))
    #print(year.month_list)
    year.month <- unlist(year.month_list)
    year.month_loc <- cbind(year.month, rep(1:N_s, N_timelags))
    year.month_unique <- unique(year.month)
    N_ymu <- length(year.month_unique)
    
    #for each unique time, record which location (which index) needs to be extracted
    year.month_index <- lapply(year.month_unique, function(ym) which(ym == year.month))
    
    #make covariate paths
    if(!is.null(covariates)){
      covariate_paths <- make_monthly_paths(covariates[i], resolution)
      #covariate_paths <- make_monthly_paths("Rain", resolution)
      #print(covariate_paths)
    }else{
      covariate_paths <- list(startpaths[[i]], endpaths[[i]])
    }
    
    #for each unique time
    covariate_value_list_unique <- foreach(j=1:N_ymu,
                                           .packages=c("raster")) %do_seq_par% {
        cov_raster <- raster(paste0(covariate_paths[[1]], year.month_unique[j], covariate_paths[[2]]))
        # print(year.month_loc[year.month_index[[j]], 2])
        # print(dim(locations))
        ym_locations <- locations[as.numeric(year.month_loc[year.month_index[[j]], 2]), , drop=FALSE]
        if(!is.null(crop_extent)) cov_raster <- crop(cov_raster, crop_extent)
        if(!is.null(n_blur)) cov_raster <- focal(cov_raster, w = matrix(1,n_blur,n_blur), mean, na.rm=T)
        cov_vals <- extract(cov_raster, SpatialPoints(ym_locations))
        cov_vals
     }
    if(run_parallel) stopCluster(cl)
    
    # covariate_value_list_unique <- list()
    # for(j in 1:N_ymu){
    #   cov_raster <- raster(paste0(covariate_paths[[1]], year.month_unique[j], covariate_paths[[2]]))
    #   # print(year.month_loc[year.month_index[[j]], 2])
    #   # print(dim(locations))
    #   ym_locations <- locations[as.numeric(year.month_loc[year.month_index[[j]], 2]), , drop=FALSE]
    #   if(!is.null(crop_extent)) cov_raster <- crop(cov_raster, crop_extent)
    #   if(!is.null(n_blur)) cov_raster <- focal(cov_raster, w = matrix(1,n_blur,n_blur), mean, na.rm=T)
    #   cov_vals <- extract(cov_raster, SpatialPoints(ym_locations))
    #   covariate_value_list_unique[[j]] <- cov_vals
    # }
    
    
    #put covariate values into big vector
    covariate_vals_vec <- rep(NA, N_s * N_timelags)
    for(j in 1:N_ymu){
      covariate_vals_vec[year.month_index[[j]]] <- covariate_value_list_unique[[j]]
    }
    
    covariate_vals_list <- vector_split(covariate_vals_vec, N_parts = N_timelags)
    names(covariate_vals_list) <- paste0(covariates[i], "_", timelags[[i]])
    out_vals <- c(out_vals, covariate_vals_list)
  }
  
  return(out_vals)
}





# ########test
mad_raster <- raster(paste0(zdrive_path,
                            "Madagascar_incidence/CovariateData/5km/Monthly_variables/LST_day/LST_Day.2010.03.mean.5km.mean.tif"))
mad_coords <- coordinates(mad_raster)


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
#                             "10km",
#                             n_blur = 5,
#                             crop_extent=extent(mad_raster),
#                             timelags=list(c(0, 1), 0, 0, 0))
# set.seed(12)
# ptm <- proc.time()
# n_coords <- 10
# n_times <- 1
# years_use <- sample(2013:2016, n_times, replace=T)
# months_use <- sample(1:12, n_times, replace=T)
# mad_coords_use <- mad_coords[sample.int(dim(mad_coords)[1], n_coords), ]
# test2 <- monthly_covariates(mad_coords_use,
#                             years_use,
#                             months_use,
#                             c("Rain", "LST_day"),
#                             "5km",
#                             n_blur = 5,
#                             crop_extent=extent(mad_raster),
#                             timelags=list(c(0, 1), 0))
# print(proc.time() - ptm)
# 
# ptm <- proc.time()
# test2 <- monthly_covariates(mad_coords_use,
#                             years_use,
#                             months_use,
#                             c("Rain", "LST_day"),
#                             "5km",
#                             n_blur = 5,
#                             crop_extent=extent(mad_raster),
#                             timelags=list(c(0, 1), 0),
#                             run_parallel=TRUE,
#                             n_cores=2)
# print(proc.time() - ptm)


##test monthly individual
set.seed(12)
#ptm <- proc.time()
n_coords <- 10
#n_times <- 1
years_use <- sample(2013:2016, n_coords, replace=T)
months_use <- sample(1:12, n_coords, replace=T)
mad_coords_use <- mad_coords[sample.int(dim(mad_coords)[1], n_coords), ]
test2 <- monthly_covariates_individual(mad_coords_use,
                            years_use,
                            months_use,
                            c("Rain", "LST_day"),
                            "5km",
                            n_blur = 5,
                            crop_extent=extent(mad_raster),
                            timelags=list(c(0, 1), 0),
                            run_parallel = TRUE,
                            n_cores=2)
test3 <- monthly_covariates(mad_coords_use[1:3, ],
                            years_use[3],
                            months_use[3],
                            c("Rain", "LST_day"),
                            "5km",
                            n_blur = 5,
                            crop_extent=extent(mad_raster),
                            timelags=list(c(0, 1), 0),
                            run_parallel=TRUE,
                            n_cores=2)
