
setClass("Covariate",
         slots = c(
           Token_0 = "character",
           Token_1 = "character",
           Token_2 = "character",
           Token_3 = "character",
           Token_4 = "character",
           Token_5 = "character",
           covariate_index = "numeric",
           folder_path = "character",
           spatial_resolution = "character",
           temporal_resolution = "character",
           interpolate = "logical",
           format = "character"
         )
         )

#dictionaries
formats <- c("year-day", "year-month", "year")
Token_0_dict <- c("chirps-v2-0", 
                     "LST_Day_v6",
                     "LST_Night_v6",
                     "TCB_v6.",
                     "EVI_v6.")
covariate_names <- c("Rain",
                      "LST_Day",
                      "LST_Night",
                      "TCB",
                      "EVI")
names(Token_0_dict) <- covariate_names
folder_paths <- c("Z:/mastergrids/Other_Global_Covariates/Rainfall/CHIRPS",
                  "Z:/mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day",
                  "Z:/mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Night",
                  "Z:/mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/TCB_v6",
                  "Z:/mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/EVI_v6"
                  )
spatial_resolutions <- c("1km", "5km", "10km")
modis_resolution_spatial_summaries <- c("Data.tif", "mean.tif", NA)
spatial_summaries <- list(c("NN.tif", "NN.tif", "Data.tif"),
                          modis_resolution_spatial_summaries,
                          modis_resolution_spatial_summaries,
                          modis_resolution_spatial_summaries,
                          modis_resolution_spatial_summaries)
temporal_resolutions <- c("8-Daily", "Monthly", "Annual")
modis_temp_summary_defaults <- c("Data", "mean", "mean")
temporal_summary_defaults <- list(c(NA, "sum", "sum"),
                                  modis_temp_summary_defaults,
                                  modis_temp_summary_defaults,
                                  modis_temp_summary_defaults,
                                  modis_temp_summary_defaults)
                               

test <- new("Covariate", Token_0="1")

setGeneric("name", function(x) standardGeneric("name"))
setGeneric("name<-", function(x, value) standardGeneric("name<-"))

setMethod("name", "Covariate", function(x) names(Token_0_dict)[which(Token_0_dict==x@Token_0)])

setMethod("name<-", "Covariate", function(x, value) {
  x@Token_0 <- Token_0_dict[which(names(Token_0_dict) == value)]
  x
})


Covariate <- function(name,
                      spatial_resolution,
                      #temporal=TRUE,
                      temporal_resolution,
                      format = NULL,
                      interpolate=FALSE
){
  name <- match.arg(name, covariate_names)
  sptial_resolution <- match.arg(spatial_resolution, spatial_resolutions)
  temporal_resolution <- match.arg(temporal_resolution, temporal_resolutions)
  #set up which covariate it is
  sp_resolution_index <- which(spatial_resolutions == spatial_resolution)
  temp_resolution_index <- which(temporal_resolutions == temporal_resolution)
  covariate_index <- which(covariate_names == name)
  Token_0 <- Token_0_dict[[covariate_index]]
  Token_3 <- temporal_summary_defaults[[covariate_index]][temp_resolution_index]
  Token_5 <- spatial_summaries[[covariate_index]][sp_resolution_index]
  folder_path <- folder_paths[[covariate_index]]
  if(is.null(format)){
    format <- formats[temp_resolution_index]
  }else{
    format <- match.arg(format, formats)
  }
  
  #check if such a resolution exists for this covariate
  if(is.na(Token_5)) stop(paste0("resolution ", spatial_resolution, " does not exist for covariate ", name))
  if(is.na(Token_3)) stop(paste0("resolution ", temporal_resolution, " does not exist for covariate ", name))
  
  new("Covariate", 
      Token_0=Token_0, 
      Token_3=Token_3,
      Token_4=spatial_resolution,
      Token_5=Token_5,
      interpolate=interpolate,
      covariate_index=covariate_index, 
      folder_path=folder_path,
      spatial_resolution=spatial_resolution,
      temporal_resolution=temporal_resolution,
      format=format)
}

rain_5k_test <- Covariate("Rain", "5km", "Monthly")
#rain_10k_test <- Covariate("Rain", "10km", "8-Daily")
#temp_10k_test <- Covariate("LST_Day", "10km", "Annual")
temp_1k_test <- Covariate("LST_Day", "1km", "Annual")


setGeneric("makeFilePath", function(covariate, ...) standardGeneric("makeFilePath"))
setMethod("makeFilePath", "Covariate", function(covariate, year){
  #if temporal resolution this will work fine
  if(covariate@temporal_resolution == "Annual"){
    covariate@Token_1 <- as.character(year)
    covariate@Token_2 <- "Annual"
  }else{
    stop("full date not supplied")
  }
  return(pasteFilePath(covariate))
})

setMethod("makeFilePath", "Covariate", function(covariate, year, Token_2){
  #if temporal resolution is annual throw an error
  if(covariate@temporal_resolution == "Annual"){
    stop("date or month supplied to Annual covariate")
  }else{
    
  }
  return(pasteFilePath(covariate))
})


pasteFilePath <- function(covariate){
  return(paste0(
    paste(covariate@folder_path,
          covariate@spatial_resolution,
          covariate@temporal_resolution,
          "",
          sep="/"),
    paste(covariate@Token_0,
          covariate@Token_1,
          covariate@Token_2,
          covariate@Token_3,
          covariate@Token_4,
          covariate@Token_5,
           sep=".")
    )
  )
}


#makeFilePath(rain_5k_test, 2011)
print(makeFilePath(temp_1k_test, 2011))
#print(pasteFilePath(temp_1k_test))
