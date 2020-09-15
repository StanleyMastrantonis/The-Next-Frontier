############################################################################
#################R script for processing and modellling#####################
#################Human marine settlements###################################
library(raster)
library(rgdal)
library(gdalUtils)
library(snow)
library(ncdf4)
library(dplyr)
library(colorRamps)

raw_data_dir = "F:/City Model/Data dryad/Raw data"  #### set directory to raw data 
process_data_dir ='F:/City Model/Constraint & Opportunity Rasters'### set directory for processed rasters
output_data_dir = 'F:/City Model/Output' ### set directory for final model results
temp_dir = 'F:/City Model/Temp folder' ###temp working directory
  
options(scipen = 999)

##############functions##################################
########################################################

#please run the code for all these custom functions first-----------
find_func = function(seq, vector){
  y = which(abs(seq-vector) == min(abs(seq-vector)))
  y[1]
}

sigmoid_func = function(x, coef, d) {
  1 / (1 + exp((x-d)/coef))
}

sigmoid_func_2 = function(x, coef, d, max.add) {
  max.add / (1 + exp((x-d)/coef))
}

lin_fun = function(x,m,c){
  m*x+c
}

scale_func =  function(x) { (x-min(x))/(max(x)-min(x))}

scale_func_2 =  function(x) { 2*((x-min(x))/(max(x)-min(x)))-1}

gaus_func = function(x, a ,b, c ) {a*(2.71828^-(((x-b)^2)/(2*(c)^2)))}

mask_function = function(x,mask){
  ex = extent(mask)
  res = resample(x,mask ,method='bilinear')
  extend = extend(res,ex, value = NA)
  crop = crop(extend, mask, snap = 'near')
  mask = mask(crop, mask,maskvalue=TRUE)
  return(mask)
}

find_func = function(seq, vector){
  y = which(abs(seq-vector) == min(abs(seq-vector)))
  y[1]
}


##############################SECTION 1############################
###################################################################
#This section is for the prcoessing and transformation of the raw data
#You can skip to section 2 to use the processed data provided to run the
#Weighted sum model

##############################DEM and CRS##########################
##################################################################
setwd(raw_data_dir)
dem = raster('dem_moll.tif')
dem[dem > 0] = 0
CoordRS = '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
e = extent(dem)

dem_seq = seq(minValue(dem), maxValue(dem), 1)
dem_seq = seq(-2000, maxValue(dem), 1)
dem_sig = sigmoid_func_2(dem_seq,-100,-500,1)
dem_df = data.frame(dem_seq, dem_sig)
dem_df$col_seq = seq(1, nrow(dem_df),1)
dem_sub = subset(dem_df, dem_df$dem_seq >-2000)

pal = colorRampPalette(c('pink','yellow','dark green'))
pal2 = colorRampPalette(c('gray96','pink','yellow','red'))

   
plot(dem_sub$dem_seq, dem_sub$dem_sig, type ='n', lwd = 4
     ,xlab = 'Elevation (m)', ylab = 'Opportunity', main = 'Elevation Opportunity')

lines(dem_sub$dem_seq, dem_sub$dem_sig, lwd = 7, col = "black")
mapply(function(x, y, col) lines(x, y, col = col, lwd = 5),
       x = rbind.data.frame(dem_sub$dem_seq, dplyr::lag(dem_sub$dem_seq)), 
       y = rbind.data.frame(dem_sub$dem_sig, dplyr::lag(dem_sub$dem_sig)), 
       col = (pal(length(rev(dem_sub$col_seq)))))


dem_calc = calc(dem, function(x) {sigmoid_func_2(x,-100,-500,1)})
plot(dem_calc,xlim = c(-20010000,20010000), xaxt='n', yaxt='n', main = 'Elevation Opportunity',
     legend.width=1, legend.shrink=0.75,legend.args=list(text='Opportunity', side=4, font=2, line=2.5, cex=0.8)
     )

plot(dem_calc)

writeRaster(dem_calc,paste0(process_data_dir,'/dem.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)


#################SST Averages#########################
######################################################

nc_data = nc_open(list.files(path = getwd(), pattern='.*sst.*\\.nc$', 
                              all.files=TRUE, full.names=FALSE, ignore.case = TRUE), write = TRUE)

lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
t = ncvar_get(nc_data, "time")

sst_array = ncvar_get(nc_data, "sst")
sst_slice = sst_array[, , 1535:1547] 
dim(sst_slice)

stack_sst = stack(lapply(1:dim(sst_slice)[3], function(x){ 
  rotate(raster(t(sst_slice[,,x]), xmn=min(lon), xmx=max(lon), 
        ymn=min(lat), ymx=max(lat),crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
     )     
    }))

sst_mu = calc(stack_sst, mean, na.rm=T)

######################################
#important-------------------------
#if you dont have gdal you need to project the raster using the below code
#it will take some time, otherwise skip to the next segment to use gdal
#the rest the script uses gdal
#warnings are ok

sst_proj = projectRaster(sst_mu, crs = CoordRS )

####################################

writeRaster(sst_mu, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)

src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile= paste0(temp_dir, "/sst_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)


sst_proj = raster(paste0(temp_dir,"/sst_temp_proj.tif"))
plot(sst_proj)
beginCluster(7)
sst_mask = sst_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()
gc()
plot(sst_mask)

sst_seq = seq(minValue(sst_mask), maxValue(sst_mask), 0.1)
sst_gaus = gaus_func(sst_seq, 1, 18, 7)
sst_df = data.frame(sst_seq, sst_gaus)
peak = find_func(sst_df$sst_seq, 18)
sst_df$col_seq = c(seq(1, peak,1), (seq(peak,1,length.out = nrow(sst_df)-peak)))

plot(sst_df$sst_seq, sst_df$sst_gaus, type ='n', 
     xlab = 'Mean SST(°)', ylab = 'Opportunity', main = 'Temperature Opportunity')


mapply(function(x, y, col) lines(x, y, col = col, lwd = 5),
        x = rbind.data.frame(sst_df$sst_seq[1:peak], dplyr::lag(sst_df$sst_seq[1:peak])), 
        y = rbind.data.frame(sst_df$sst_gaus[1:peak], dplyr::lag(sst_df$sst_gaus[1:peak])), 
        col = (pal(peak)))

mapply(function(x, y, col) lines(x, y, col = col, lwd = 5),
       x = rbind.data.frame(sst_df$sst_seq[peak:nrow(sst_df)], dplyr::lag(sst_df$sst_seq[peak:nrow(sst_df)])), 
       y = rbind.data.frame(sst_df$sst_gaus[peak:nrow(sst_df)], dplyr::lag(sst_df$sst_gaus[peak:nrow(sst_df)])), 
       col = rev(pal(nrow(sst_df)-peak+1)))

sst_calc = calc(sst_mask, function(x) {gaus_func(x, 1, 18, 7)})

writeRaster(sst_calc,paste0(process_data_dir,'/sst.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)


nc_close(nc_data)
gc()
rm(sst_proj, sst_mu, sst_slice, sst_array, sst_calc)
##################ICE Averages#########################
#######################################################
nc_data = nc_open(list.files(path = getwd(), pattern='.*icec.*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE), write = TRUE)


lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
t = ncvar_get(nc_data, "time")

ice_array = ncvar_get(nc_data, "icec")
ice_slice = ice_array[, ,11:12 ] 


stack_ice = stack(lapply(1:dim(ice_slice)[3], function(x){ 
  rotate(raster(t(ice_slice[,,x]), xmn=min(lon), xmx=max(lon), 
         ymn=min(lat), ymx=max(lat),crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  )     
}))


ice_mu = calc(stack_ice,mean, na.rm=T)

writeRaster(ice_mu, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)
src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(temp_dir,"ice_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

ice_proj = raster(paste0(temp_dir,'/ice_temp_proj.tif'))

beginCluster(7)
ice_mask = ice_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()
gc()
plot(ice_mask)

ice_seq = seq(minValue(ice_mask), maxValue(ice_mask), 0.01)
ice_sig = sigmoid_func_2(ice_seq,-0.08,0.5,1)
ice_df = data.frame(ice_seq, ice_sig)
ice_df$col_seq = seq(1, nrow(ice_df),1)

plot(ice_df$ice_seq, ice_df$ice_sig, type ='n', lwd = 4
     ,xlab = 'Sea Ice Concentration', ylab = 'Constraint', main = 'Sea Ice Constraint')

mapply(function(x, y, col) lines(x, y, col = col, lwd = 5),
       x = rbind.data.frame(ice_df$ice_seq, dplyr::lag(ice_df$ice_seq)), 
       y = rbind.data.frame(ice_df$ice_sig, dplyr::lag(ice_df$ice_sig)), 
       col = (pal2(length(rev(ice_df$col_seq)))))


ice_calc = calc(ice_mask, function(x) {sigmoid_func_2(x,-0.08,0.5,1)})



writeRaster(ice_calc,paste0(process_data_dir,'/icec.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)

nc_close(nc_data)
gc()
rm(ice_proj, ice_mu, ice_slice, ice_array, ice_calc)

##################CHLA Averages#########################
#######################################################


chla_list = list.files(path = getwd(), pattern='.*chlor_a.*\\.nc$', 
           all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

nc_data = lapply(1:length(chla_list), function(x) {
                 nc_open(chla_list[x])})

chla_array_list = lapply(1:length(nc_data), function(x){
                        ncvar_get(nc_data[[x]], "chlor_a")
                        })

lon = ncvar_get(nc_data[[1]],"lon") 
lat = ncvar_get(nc_data[[1]], "lat", verbose = F)


stack_chla = stack(lapply(1:length(chla_array_list), function(x){ 
  raster(t(chla_array_list[[x]]), xmn=min(lon), xmx=max(lon), 
                ymn=min(lat), ymx=max(lat),crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
      
)}))


chla_mu = calc(stack_chla,mean, na.rm=T)


writeRaster(chla_mu, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)
src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(temp_dir,"/chla_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

chla_proj = raster(paste0(temp_dir,'/chla_temp_proj.tif'))

beginCluster(7)
chla_mask = chla_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(chla_mask)

chla_seq = seq(minValue(chla_mask), maxValue(chla_mask), 1)
chla_sig = sigmoid_func_2(chla_seq,-6,10,1)
chla_df = data.frame(chla_seq, chla_sig)
chla_df$col_seq = seq(1, nrow(chla_df),1)

plot(chla_df$chla_seq, chla_df$chla_sig, type ='n', lwd = 4
     ,xlab = 'Chl-a Concentration(mg/m^3)', ylab = 'Opportunity', main = 'Chl-a Opportunity')
mapply(function(x, y, col) lines(x, y, col = col, lwd = 5),
       x = rbind.data.frame(chla_df$chla_seq, dplyr::lag(chla_df$chla_seq)), 
       y = rbind.data.frame(chla_df$chla_sig, dplyr::lag(chla_df$chla_sig)), 
       col = rev(green2red(length(rev(chla_df$col_seq)))))


chla_calc = calc(chla_mask, function(x) {sigmoid_func_2(x,-6,10,1)})

writeRaster(chla_calc,  paste0(process_data_dir,'/chla_a.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)


nc_close(nc_data)
gc()
rm(chla_proj, chla_mu, chla_array_list,stack_chla, chla_calc, chla_mask )

##################biodiversity################
##############################################


nc_data = nc_open(list.files(path = getwd(), pattern='.*biod*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))

biod_array = ncvar_get(nc_data, "biod")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)

biod_ras = raster(t(biod_array), xmn=min(lon), xmx=max(lon), 
          ymn=min(lat), ymx=max(lat),
          crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
     
writeRaster(biod_ras, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)
src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(temp_dir,"/biod_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

biod_proj = raster(paste0(temp_dir,'/biod_temp_proj.tif'))

beginCluster(7)
biod_mask = biod_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()


plot(biod_mask)

biod_seq = seq(minValue(biod_mask), maxValue(biod_mask), 0.01)
biod_sig = sigmoid_func_2(biod_seq,-0.09,0.5,1)
biod_df = data.frame(biod_seq, biod_sig)
biod_df$col_seq = seq(1, nrow(biod_df),1)

plot(biod_df$biod_seq, biod_df$biod_sig, type ='l', lwd = 4
     ,xlab = 'Marine Biodiversity', ylab = 'Constraint', main = 'Biodiversity Constraint')

mapply(function(x, y, col) lines(x, y, col = col, lwd = 5),
       x = rbind.data.frame(biod_df$biod_seq, dplyr::lag(biod_df$biod_seq)), 
       y = rbind.data.frame(biod_df$biod_sig, dplyr::lag(biod_df$biod_sig)), 
       col = (green2red(length(rev(biod_df$col_seq)))))


biod_calc = calc(biod_mask, function(x) {sigmoid_func_2(x,-0.09,0.5,1)})

writeRaster(biod_calc, paste0(process_data_dir,'/biod.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)

nc_close(nc_data)
gc()
rm(biod_proj, biod_array , biod_calc, biod_ras)


####################shipping#####################
################################################

nc_data = nc_open(list.files(path = getwd(), pattern='.*shipping*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))
ship_array = ncvar_get(nc_data, "ship")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
dim(ship_array) 

ship = raster(t(ship_array), xmn=min(lon), xmx=max(lon), 
                       ymn=min(lat), ymx=max(lat),crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

gc()
writeRaster(ship, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)
gc()
src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(temp_dir,"/ship_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

ship_proj = raster(paste0(temp_dir,'/ship_temp_proj.tif'))


ship_mask = aggregate(ship_proj, c(2.2351554,2.2351554))

beginCluster(7)
ship_mask_ = ship_mask %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(ship_mask_)

ship_seq = seq(minValue(ship_mask_), maxValue(ship_mask_), 0.01)
ship_sig = sigmoid_func_2(ship_seq,-0.09,0.5,1)
ship_df = data.frame(ship_seq, ship_sig)
ship_df$col_seq = seq(1, nrow(ship_df),1)

plot(ship_df$ship_seq, ship_df$ship_sig, type ='l', lwd = 4
     ,xlab = 'Shipping Density', ylab = 'Constraint', main = 'Shipping Constraint')
mapply(function(x, y, col) lines(x, y, col = col, lwd = 5),
       x = rbind.data.frame(ship_df$ship_seq, dplyr::lag(ship_df$ship_seq)), 
       y = rbind.data.frame(ship_df$ship_sig, dplyr::lag(ship_df$ship_sig)), 
       col = rev(green2red(length(rev(ship_df$col_seq)))))


ship_calc = calc(ship_mask_, function(x) {sigmoid_func_2(x,-0.09,0.5,1)})
plot(ship_calc)
writeRaster(ship_calc, paste0(process_data_dir,'/shipping.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)



nc_close(nc_data)
gc()
rm(ship_proj,ship, ship_calc, ship_mask, ship_mask_ )

###############coast distance##################
##############################################

setwd("F:/City Model/Dist to coast")
nc_data = nc_open(list.files(path = getwd(), pattern='.*dist*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))

dist_array = ncvar_get(nc_data, "dist")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
dim(ship_array) 

dist = raster(t(dist_array), xmn=min(lon), xmx=max(lon), 
              ymn=min(lat), ymx=max(lat),crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

gc()
writeRaster(dist, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)
gc()
src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(temp_dir,"/dist_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

dist_proj = raster(paste0(temp_dir,'/dist_temp_proj.tif'))

gc()

beginCluster(7)
dist_mask = dist_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(dist_mask)

dist_seq = seq(minValue(dist_mask), maxValue(dist_mask), 1000)
dist_lin = lin_fun(dist_seq ,-1/4153133.50,1)
dist_df = data.frame(dist_seq, dist_lin)
dist_df$col_seq = seq(1, nrow(dist_df),1)

plot(dist_df$dist_seq, dist_df$dist_lin, type ='l', lwd = 4
     ,xlab = 'Distance to Coast (m)', ylab = 'Opportunity', main = 'Coastal Distance Opportunity')

mapply(function(x, y, col) lines(x, y, col = col, lwd = 5),
       x = rbind.data.frame(dist_df$dist_seq, dplyr::lag(dist_df$dist_seq)), 
       y = rbind.data.frame(dist_df$dist_lin, dplyr::lag(dist_df$dist_lin)), 
       col = (green2red(length(rev(dist_df$col_seq)))))


dist_calc = calc(dist_mask, function(x) {lin_fun(x,-1/4153133.50,1)})
plot(dist_calc)
writeRaster(dist_calc, paste0(process_data_dir,'/dist.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)



nc_close(nc_data)
gc()
rm(dist_proj,dist, dist_calc, dist_mask )


##################wind#########################
###############################################


nc_data = nc_open(list.files(path = getwd(), pattern='.*wind*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))
wind_array = ncvar_get(nc_data, "sfcWind")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
dim(wind_array) 

wind_slice = wind_array[, ,110:123 ] 

stack_wind = stack(lapply(1:dim(wind_slice)[3], function(x){ 
  rotate(raster(t(wind_slice[,,x]), xmn=min(lon), xmx=max(lon), 
                ymn=min(lat), ymx=max(lat),crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  )     
}))


wind_mu = stack_wind %>%
    calc(mean, na.rm=T)%>%
     flip(direction = 'y')

plot(wind_mu)
gc()
writeRaster(wind_mu, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)
gc()

src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(temp_dir,"/wind_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

wind_proj = raster(paste0(temp_dir,'/wind_temp_proj.tif'))

gc()

beginCluster(7)
wind_mask = wind_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(wind_mask)

wind_seq = seq(minValue(wind_mask), maxValue(wind_mask), 0.1)
wind_sig = sigmoid_func_2(wind_seq,-0.5,5,1)
wind_df = data.frame(wind_seq, wind_sig)
wind_df$col_seq = seq(1, nrow(wind_df),1)

plot(wind_df$wind_seq, wind_df$wind_sig, type ='l', lwd = 4
     ,xlab = 'Wind Speed(m/s)', ylab = 'Opportunity', main = 'Wind Speed Opportunity')

mapply(function(x, y, col) lines(x, y, col = col, lwd = 5),
       x = rbind.data.frame(wind_df$wind_seq, dplyr::lag(wind_df$wind_seq)), 
       y = rbind.data.frame(wind_df$wind_sig, dplyr::lag(wind_df$wind_sig)), 
       col = rev(green2red(length(rev(wind_df$col_seq)))))


wind_calc = calc(wind_mask, function(x) {sigmoid_func_2(x,-0.5,5,1)})
plot(wind_calc)
writeRaster(wind_calc, paste0(process_data_dir,'/wind.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)

nc_close(nc_data)
gc()
rm(wind_calc,wind_array, wind_mu, wind_proj, stack_wind )



#################solar insol#####################
#################################################

nc_data = nc_open(list.files(path = getwd(), pattern='.*solar*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))

solar_array = ncvar_get(nc_data, "solar")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
dim(solar_array) 


sol_ras = raster(t(solar_array), xmn=min(lon), xmx=max(lon), 
                  ymn=min(lat), ymx=max(lat),
                  crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))


writeRaster(sol_ras, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)
gc()

src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(temp_dir,"sol_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

sol_proj = raster(paste0(temp_dir,'/sol_temp_proj.tif'))

gc()

beginCluster(7)
sol_mask = sol_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(sol_mask)

sol_seq = seq(minValue(sol_mask), maxValue(sol_mask), 0.1)
sol_sig = sigmoid_func_2(sol_seq,-4,20,1)
sol_df = data.frame(sol_seq, sol_sig)
sol_df$col_seq = seq(1, nrow(sol_df),1)

plot(sol_df$sol_seq, sol_df$sol_sig, type ='l', lwd = 4
     ,xlab = 'Solar Insolation(w/m2)', ylab = 'Opportunity', main = 'Solar Opportunity')

mapply(function(x, y, col) lines(x, y, col = col, lwd = 5),
       x = rbind.data.frame(sol_df$sol_seq, dplyr::lag(sol_df$sol_seq)), 
       y = rbind.data.frame(sol_df$sol_sig, dplyr::lag(sol_df$sol_sig)), 
       col = rev(green2red(length(rev(sol_df$col_seq)))))


sol_calc = calc(sol_mask, function(x) {sigmoid_func_2(x,-4,20,1)})
plot(sol_calc)
writeRaster(sol_calc, paste0(process_data_dir,'/solar.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)

nc_close(nc_data)
gc()
rm(sol_calc,sol_ras, sol_proj )

##################fishery data#################
##############################################

nc_data = nc_open(list.files(path = getwd(), pattern='.*fish_sum*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))

fish_array = ncvar_get(nc_data, "landings")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
dim(fish_array) 

fish_ras = raster(t(fish_array), xmn=min(lon), xmx=max(lon), 
                 ymn=min(lat), ymx=max(lat),
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
plot(fish_ras)

writeRaster(fish_ras, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)

src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(temp_dir,"/fish_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

fish_proj = raster(paste0(temp_dir,'/fish_temp_proj.tif'))

gc()

beginCluster(7)
fish_mask = fish_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(fish_mask)

fish_seq = seq(minValue(fish_mask), maxValue(fish_mask), 1)
fish_sig = sigmoid_func_2(fish_seq,-50,20,1)
fish_df = data.frame(fish_seq, fish_sig)
fish_df$col_seq = seq(1, nrow(fish_df),1)

plot(fish_df$fish_seq, fish_df$fish_sig, type ='l', lwd = 4
     ,xlab = 'Landings (tonnes)', ylab = 'Constraint', main = 'Fishery COnstraint')

mapply(function(x, y, col) lines(x, y, col = col, lwd = 5),
       x = rbind.data.frame(fish_df$fish_seq, dplyr::lag(fish_df$fish_seq)), 
       y = rbind.data.frame(fish_df$fish_sig, dplyr::lag(fish_df$fish_sig)), 
       col = (green2red(length(rev(fish_df$col_seq)))))


fish_calc = calc(fish_mask, function(x) {sigmoid_func_2(x,-50,20,1)})

writeRaster(fish_calc, paste0(process_data_dir,'/fish.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)


nc_close(nc_data)
gc()
rm(fish_calc,fish_ras_list, fish_sum, fish_stack, fish_mask )


#################pollution#####################
###############################################


nc_data = nc_open(list.files(path = getwd(), pattern='.*pollution*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))
pol_array = ncvar_get(nc_data, "pollution")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
dim(pol_array) 


pol_ras = raster(t(pol_array), xmn=min(lon), xmx=max(lon), 
                 ymn=min(lat), ymx=max(lat),
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))


plot(pol_ras)

gc()
writeRaster(pol_ras, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)
gc()
src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset, paste0(temp_dir,"/pol_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

pol_proj = raster(paste0(temp_dir,'/pol_temp_proj.tif'))

gc()

beginCluster(7)
pol_mask = pol_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(pol_mask)

pol_seq = seq(minValue(pol_mask), maxValue(pol_mask), 0.01)
pol_sig = sigmoid_func_2(pol_seq,-0.2,0.5,1)
pol_df = data.frame(pol_seq, pol_sig)
pol_df$col_seq = seq(1, nrow(pol_df),1)

plot(pol_df$pol_seq, pol_df$pol_sig, type ='l', lwd = 4
     ,xlab = 'Pollution Density', ylab = 'Constraint', main = 'Pollution Constraint')
mapply(function(x, y, col) lines(x, y, col = col, lwd = 5),
       x = rbind.data.frame(pol_df$pol_seq, dplyr::lag(pol_df$pol_seq)), 
       y = rbind.data.frame(pol_df$pol_sig, dplyr::lag(pol_df$pol_sig)), 
       col = (green2red(length(rev(pol_df$col_seq)))))


pol_calc = calc(pol_mask, function(x) {sigmoid_func_2(x,-0.2,0.5,1)})

writeRaster(pol_calc, paste0(process_data_dir,'/pollution.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)

nc_close(nc_data)
gc()
rm(pol_calc,pol_ras, pol_proj, pol_array )


##################tsunami###################
###########################################

nc_data = nc_open(list.files(path = getwd(), pattern='.*tsunami*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))

tsu_array = ncvar_get(nc_data, "tsunami")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
dim(tsu_array) 


tsu_ras = raster(t(tsu_array), xmn=min(lon), xmx=max(lon), 
                 ymn=min(lat), ymx=max(lat),
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))


plot(tsu_ras)

gc()
writeRaster(tsu_ras, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)
gc()
src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(temp_dir,"/tsu_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

tsu_proj = raster( paste0(temp_dir,"/tsu_temp_proj.tif"))
plot(tsu_proj)
gc()

beginCluster(7)
tsu_mask = tsu_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(tsu_mask)

tsu_seq = seq(minValue(tsu_mask), maxValue(tsu_mask), 0.1)
tsu_sig = lin_fun(tsu_seq,-100,1)
tsu_df = data.frame(tsu_seq, tsu_sig)
tsu_df$col_seq = seq(1, nrow(tsu_df),1)

plot(tsu_df$tsu_seq, tsu_df$tsu_sig, type ='l', lwd = 4
     ,xlab = 'Distance to Tsunami Events (Degrees)', ylab = 'Constraint', main = 'Tsunami Constraint')
mapply(function(x, y, col) lines(x, y, col = col, lwd = 5),
       x = rbind.data.frame(tsu_df$tsu_seq, dplyr::lag(tsu_df$tsu_seq)), 
       y = rbind.data.frame(tsu_df$tsu_sig, dplyr::lag(tsu_df$tsu_sig)), 
       col = rev(green2red(length(rev(tsu_df$col_seq)))))


tsu_calc = calc(tsu_mask, function(x) {sigmoid_func_2(x,1,2,1)})
plot(tsu_calc)
writeRaster(tsu_calc, paste0(process_data_dir,'/tsunami.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)

nc_close(nc_data)
gc()
rm(tsu_calc,tsu_ras, tsu_proj, tsu_array )


################hard bottom###################
##############################################

nc_data = nc_open(list.files(path = getwd(), pattern='.*hard_bottom*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))

hb_array = ncvar_get(nc_data, "hard_bottom")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
dim(hb_array) 


hb_ras = raster(t(hb_array), xmn=min(lon), xmx=max(lon), 
                 ymn=min(lat), ymx=max(lat),
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))


plot(hb_ras)

gc()
writeRaster(hb_ras, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)
gc()
src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(temp_dir,"/hb_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

hb_proj = raster(paste0(temp_dir,"/hb_temp_proj.tif"))
plot(hb_proj)
gc()

beginCluster(7)
hb_mask = hb_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(hb_mask)

writeRaster(hb_mask, paste0(process_data_dir,'/hard_bottom.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)


################hard shelf######################
################################################


nc_data = nc_open(list.files(path = getwd(), pattern='.*hard_shelf*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))
hs_array = ncvar_get(nc_data, "hard_shelf")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
dim(hs_array) 


hs_ras = raster(t(hs_array), xmn=min(lon), xmx=max(lon), 
                ymn=min(lat), ymx=max(lat),
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))


plot(hs_ras)

gc()
writeRaster(hs_ras, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)
gc()
src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(temp_dir,"/hs_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

hs_proj = raster(paste0(temp_dir,"/hs_temp_proj.tif"))
plot(hs_proj)
gc()

beginCluster(7)
hs_mask = hs_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(hs_mask)

writeRaster(hs_mask, paste0(process_data_dir,'/hard_shelf.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)


################soft shelf##################
############################################


nc_data = nc_open(list.files(path = getwd(), pattern='.*soft_shelf*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))

ss_array = ncvar_get(nc_data, "soft_shelf")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
dim(ss_array) 


ss_ras = raster(t(ss_array), xmn=min(lon), xmx=max(lon), 
                ymn=min(lat), ymx=max(lat),
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))


plot(ss_ras)

gc()
writeRaster(ss_ras, paste0(temp_dir,'/temp.tif', 'GTiff'), prj = TRUE, overwrite = TRUE)
gc()
src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(tempdir,"/ss_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

ss_proj = raster( paste0(tempdir,"/ss_temp_proj.tif"))
plot(ss_proj)
gc()

beginCluster(7)
ss_mask = ss_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(ss_mask)

writeRaster(ss_mask, paste0(process_data_dir,'/soft_shelf.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)


##############seamounts################
######################################


nc_data = nc_open(list.files(path = getwd(), pattern='.*seamounts*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))
sm_array = ncvar_get(nc_data, "seamounts")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
dim(ss_array) 


sm_ras = raster(t(sm_array), xmn=min(lon), xmx=max(lon), 
                ymn=min(lat), ymx=max(lat),
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))


plot(sm_ras)

gc()
writeRaster(sm_ras, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)
gc()
src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(temp_dir,"/sm_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

sm_proj = raster(paste0(temp_dir,"/sm_temp_proj.tif"))
plot(sm_proj)
gc()

beginCluster(7)
sm_mask = sm_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(sm_mask)

writeRaster(sm_mask, paste0(process_data_dir,'/seamounts.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)


################hotspots######################
#############################################

nc_data = nc_open(list.files(path = getwd(), pattern='.*hotspots*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))

hotspot_array = ncvar_get(nc_data, "hotspots")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
dim(hotspot_array) 


hotspot_ras = raster(t(hotspot_array), xmn=min(lon), xmx=max(lon), 
                ymn=min(lat), ymx=max(lat),
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))


plot(hotspot_ras)

gc()
writeRaster(hotspot_ras, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)
gc()
src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(temp_dir,"/hotspot_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

hotspot_proj = raster(paste0(temp_dir,"/hotspot_temp_proj.tif"))
plot(hotspot_proj)
gc()

beginCluster(7)
hotspot_mask = hotspot_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(hotspot_mask)

writeRaster(hotspot_mask, paste0(process_data_dir,'/hotspots.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)


##################EEZ####################
########################################

nc_data = nc_open(list.files(path = getwd(), pattern='.*EEZ*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))
eez_array = ncvar_get(nc_data, "eez")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
dim(eez_array) 


eez_ras = raster(t(eez_array), xmn=min(lon), xmx=max(lon), 
                     ymn=min(lat), ymx=max(lat),
                     crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))


plot(eez_ras)

gc()
writeRaster(eez_ras, paste0(temp_dir,'temp.tif', 'GTiff'), prj = TRUE, overwrite = TRUE)
gc()
src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(tempdir,"/eez_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

eez_proj = raster(paste0(tempdir,"/eez_temp_proj.tif"))
plot(eez_proj)
gc()

beginCluster(7)
eez_mask = eez_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(eez_mask)

writeRaster(eez_mask, paste0(process_data_dir,'/eez.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)


#################MPA####################
########################################

nc_data = nc_open(list.files(path = getwd(), pattern='.*mpa*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))
mpa_array = ncvar_get(nc_data, "mpa")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
dim(mpa_array) 


mpa_ras = raster(t(mpa_array), xmn=min(lon), xmx=max(lon), 
                 ymn=min(lat), ymx=max(lat),
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))


plot(mpa_ras)

gc()
writeRaster(mpa_ras, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)
gc()
src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(temp_dir,"/mpa_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

mpa_proj = raster(paste0(temp_dir,"/mpa_temp_proj.tif"))
plot(mpa_proj)
gc()

beginCluster(7)
mpa_mask = mpa_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(mpa_mask)

writeRaster(mpa_mask, paste0(process_data_dir,'/mpa.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)


#################cables###################
##########################################


nc_data = nc_open(list.files(path = getwd(), pattern='.*cables*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))
cab_array = ncvar_get(nc_data, "cables")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
dim(cab_array) 


cab_ras = raster(t(cab_array), xmn=min(lon), xmx=max(lon), 
                 ymn=min(lat), ymx=max(lat),
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))


plot(cab_ras)

gc()
writeRaster(cab_ras, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)
gc()
src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(temp_dir,"/cab_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

cab_proj = raster(paste0(temp_dir,"/cab_temp_proj.tif"))
plot(mpa_proj)
gc()

beginCluster(7)
cab_mask = cab_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(cab_mask)

writeRaster(cab_mask, paste0(process_data_dir,'/cables.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)


####################oil rigs#####################
#################################################

nc_data = nc_open(list.files(path = getwd(), pattern='.*oil*\\.nc$', 
                             all.files=TRUE, full.names=FALSE, ignore.case = TRUE))
oil_array = ncvar_get(nc_data, "oil")
lon = ncvar_get(nc_data,"lon") 
lat = ncvar_get(nc_data, "lat", verbose = F)
dim(oil_array) 


oil_ras = raster(t(oil_array), xmn=min(lon), xmx=max(lon), 
                 ymn=min(lat), ymx=max(lat),
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))


plot(oil_ras)

gc()
writeRaster(oil_ras, paste0(temp_dir,'/temp.tif'), 'GTiff', prj = TRUE, overwrite = TRUE)
gc()
src_dataset = list.files(path = temp_dir, pattern='.*temp.*\\.tif$', 
                         all.files=TRUE, full.names=TRUE, ignore.case = TRUE)

gdalwarp(src_dataset,dstfile = paste0(temp_dir,"/oil_temp_proj.tif"),
         t_srs=CoordRS,output_Raster=TRUE,
         overwrite=TRUE,verbose=FALSE)

oil_proj = raster(paste0(temp_dir,"/oil_temp_proj.tif"))
plot(oil_proj)
gc()

beginCluster(7)
oil_mask = oil_proj %>% 
  resample(dem ,method='bilinear') %>%
  extend(e, value = NA)%>%
  mask(dem,maskvalue=TRUE)
endCluster()

plot(oil_mask)

writeRaster(oil_mask, paste0(process_data_dir,'/oil.tif'),
            'GTiff', prj = TRUE, overwrite = TRUE)



##################################################################
#####################SECTION 2####################################

###################################################################
###################################################################
###########################MODEL###################################

setwd(process_data_dir)

criteria = stack(lapply(as.list(list.files(path = getwd(), pattern='.*\\.tif$', 
                       all.files=TRUE, full.names=TRUE, ignore.case = TRUE)), raster))

criteria_df = data.frame(criteria = names(criteria))
weights = c(-0.15, -0.50 , 0.15, 0.35, 0.10, -0.50, -0.10, 0.15,0.15, -0.15, -0.20, -0.50, -0.30, 
            -0.10, 0.15, -0.25, -0.20, 0.15, 0.15, -0.15, 0.15 )

criteria_df$weights = weights

criteria_weights = stack(lapply(1:nrow(criteria_df), function(x){
                       criteria[[x]] *criteria_df$weights[x] 
                        }))
gc()

model_sum = calc(criteria_weights, sum, na.rm = TRUE)
model_sum[model_sum < -1] = -1 


writeRaster(model_sum, paste0(output_data_dir,"/model_result.tif") ,'GTiff', prj = TRUE, overwrite = TRUE)


########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################





