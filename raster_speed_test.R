#get list of rasters
set.seed(12)
library(raster)
folder_path <- "Z:/mastergrids/MODIS_Global/MOD11A2_v6_LST/LST_Day/5km/Monthly/"

file_end <- "mean.5km.mean.tif$"

file_list <- list.files(folder_path, pattern=file_end, full.names = TRUE)



#mad_raster <- raster("Z:/Madagascar_incidence/CovariateData/5km/Monthly_variables/LST_day/LST_Day.2010.03.mean.5km.mean.tif")
mad_raster <- raster("Z:/CHAI/Modelling/India/National/Elevation.tif")
mad_coords <- coordinates(mad_raster)
#choose some locations 
N_raster <- 10
N_locs <- 100
locs <- mad_coords[sample.int(dim(mad_coords)[1], N_locs), ]


####without blurring
####these three seem to be the same speed
###open each raster, crop and extract values
cov_vals_1 <- list()

test1_start <- proc.time()
for(i in 1:N_raster){
  print(i)
  cov_raster <- raster(file_list[i])
  cov_raster <- crop(cov_raster, extent(mad_raster))
  cov_vals_1[[i]] <- extract(cov_raster, SpatialPoints(locs))
}
test1_end <- proc.time()
print(test1_end - test1_start)


###open each raster, stack, crop and extract values
cov_vals_2 <- list()

test2_start <- proc.time()
raster_stack <- stack(file_list[1:N_raster])
raster_stack <- crop(raster_stack, extent(mad_raster))
cov_vals_2 <- extract(raster_stack, SpatialPoints(locs))
test2_end <- proc.time()
print(test2_end - test2_start)


###open first raster, get indexes, extract rest using those
cov_vals_3 <- list()

test3_start <- proc.time()
cov_raster <- raster(file_list[1])
cov_raster <- crop(cov_raster, extent(mad_raster))
cov_vals_w_index <- extract(cov_raster, SpatialPoints(locs), cellnumbers=TRUE)
cov_vals_3[[1]] <- cov_vals_w_index[, 2]

for(i in 2:N_raster){
  print(i)
  cov_raster <- raster(file_list[i])
  cov_raster <- crop(cov_raster, extent(mad_raster))
  cov_vals_3[[i]] <- values(cov_raster)[cov_vals_w_index[, 1]]
}
test3_end <- proc.time()
print(test3_end - test3_start)

plot.new()
plot(cov_vals_1[[2]], cov_vals_2[, 2])
plot.new()
plot(cov_vals_1[[2]], cov_vals_3[[2]])


####with blurring
#N_raster <- 2
# N_locs <- 1000
n_blur <- 5
# locs <- mad_coords[sample.int(dim(mad_coords)[1], N_locs), ]
###open each raster, crop, focal and extract values
cov_vals_1 <- list()

test4_start <- proc.time()
for(i in 1:N_raster){
  print(i)
  cov_raster <- raster(file_list[i])
  cov_raster <- crop(cov_raster, extent(mad_raster))
  cov_raster <- focal(cov_raster, w=matrix(1, n_blur, n_blur), mean, na.rm=T)
  cov_vals_1[[i]] <- extract(cov_raster, SpatialPoints(locs))
}
test4_end <- proc.time()
print(test4_end - test4_start)
