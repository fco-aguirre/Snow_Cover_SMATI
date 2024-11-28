#!/usr/local/bin/R

### Script downscaling

library(raster)
library(terra)
library(optparse)



#library(rgdal)
#library('rgdal',lib.loc = '/home/faguirre/paquetes/library')

# Implementa sobre un raster la funcion dem_z ()
Z_r = function(n){
  (n-cellStats(n, stat='mean', na.rm=TRUE))/cellStats(n, stat='sd', na.rm=TRUE)
}

## Downscaling function
# llamar todas las librerias necesarias para la funcion

library(nlme) # to gls models
library(gstat) # to variagram
library(atakrig) # to AREA_to_AREA Kriging


## Downscaling Implementing

## Check the input and output files!!

option_list = list(
  make_option(c("-s", "--source"), type="character", default="Brunswick", 
              help="folder name of the images source [default= %default]", metavar="character"),
  make_option(c("-y", "--year"), type="integer", default=NULL, 
              help="year to process", metavar="character")
  #make_option(c("-h", "--help"), type="character", default=NULL, 
  #help="include parameters", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

source_t <- opt$source

#source_t <- 'SFS_2024'
setwd('..')
getwd()
setwd(paste0('./DATA/', source_t))

year <- opt$year


# open files
#year <- opt$year
#link <- '~/Documents/Data/Data_modis/Order_files/Brunswick/'

#link_1 <- '/Users/fcoj_aguirre/Documents/Articulos en trabajo/Snow_Magallanes/Scripts/Snow_Cover_SMATI/Spatial_Downscaling' ## Revisar este link!!
#setwd(link_1)
#link_2 <- getwd()
#setwd('..')
#link_3 <- getwd()

#setwd('./Outputs/Order_files/Brunswick')
#link_4 <- getwd()

  
link <- getwd() # Set source folder
#setwd(link)

# set output location / default is in Outputs File
setwd('..')
setwd('..')

Out_dir_1 <- dir.create(paste0('./Outputs/', source_t))
Out_dir_2 <- dir.create(paste0('./Outputs/', source_t,'/Downscaling_files'))
Out_dir_3 <- dir.create(paste0('./Outputs/', source_t,'/Downscaling_files/',year))

setwd(paste0('./Outputs/', source_t))
link_2 <- getwd()

Out_link <- paste0(link_2,'/Downscaling_files/',year)
#Out_link_1 <- paste0(link,'/Downscaling_files')
#Out_link_2 <- paste0(link,'/Downscaling_files/',year)
#Out_dir_1 <- dir.create(Out_link_1)
#Out_dir <- dir.create(Out_link)
#library(tictoc)

#link_b <- paste0(link,'Reflectance_bands/',year)
#path_mod09ga_file <- paste0(link_f,'/MOD09GA/')

setwd(link)


## Read master files
#mod09ga_file <- readLines(paste0('Reflectance_bands/',year,'/MOD09GA/',year,'_MOD09GA'))
mod09ga_file <- readLines(paste0('Reflectance_bands/',year,'/MOD09GA/',year,'_MOD09GA'))
mod09gq_file <- readLines(paste0('Reflectance_bands/',year,'/MOD09GQ/',year,'_MOD09GQ'))
mod35_file <- readLines(paste0('Reflectance_bands/',year,'/MOD35/',year,'_MOD35'))

## Revisa las bandas y selecciona días

# MOD09GA file
j <- length(mod09ga_file)
day_f <- j/8

## Lectura del archivo
day_year <- vector(mode = "list", length = day_f)
#day_time <- vector(mode = "list", length = day_f)


band_1_250_p <- raster(paste0('Reflectance_bands/',year,'/MOD09GQ/',mod09gq_file[1]))
band_1_500_p <- raster(paste0('Reflectance_bands/',year,'/MOD09GA/',mod09ga_file[2]))

## Read parameters file
parameter_file <-readLines('Parameters_Data')

#shape <- readOGR('Cuencas/Cuencas_brunswick_UTM.shp') ## rgdal is discontinued
shape <- vect(paste0('./Watershed/',parameter_file[2])) # Here we used terra package

Band_water <- raster(paste0('./land_water_mask/',parameter_file[3]))
names(Band_water) <- 'Water'

dem <- raster(paste0('./DEM/',parameter_file[4])) # realizado con clip by extent Qgis

# Mascara de agua
Band_water[Band_water$Water == 1] <- NaN  # se cambio Na
Band_water[Band_water$Water == 0] <- 1

# Ajusta resolución mascara de agua, dado los diferentes metodos de procesamiento
Band_water_res_250 <- resample(Band_water, band_1_250_p, method='bilinear') #promedio
Band_water_res_500 <- resample(Band_water, band_1_500_p, method='bilinear')

# Ajusta resolución DEM
dem_250 <- resample(dem, band_1_250_p, method='bilinear')
dem_500 <- resample(dem, band_1_500_p, method='bilinear')

# Crea coberturas de Slope and Aspet y lo adiciona a la elevación
dem_250_t <- terrain(dem_250, opt=c('slope', 'aspect'), unit='degrees')
names(dem_250) <- 'elevation'
dem_250_f <- stack(dem_250,dem_250_t)

dem_500_t <- terrain(dem_500, opt=c('slope', 'aspect'), unit='degrees')
names(dem_500) <- 'elevation'
dem_500_f <- stack(dem_500,dem_500_t)

# Define ata-pred con valores NaNs!
# Check projection
crs_1 <- crs(shape) # extract projection in UTM
band_250_st <- stack(band_1_250_p,dem_250_f) 
band_250_nan <- projectRaster(band_250_st, res=250, crs = crs_1, method = 'ngb') ## the type of projection has been modified!!


#r1 <- band_250_nan
#r2 <- rast(band_250_nan)
#r3 <- brick(r2)

ext_b <- ext(shape) ## Extract the extension of shape

R.250_nan <- crop(rast(band_250_nan), ext_b)

Band_250_nan <- mask(R.250_nan, shape)
Band_250_nan <- brick(Band_250_nan) ## return to rasterlayer object!

names(Band_250_nan) <- c('band','elevation','slope', 'aspect')

nan_band <- Band_250_nan$band * NaN


## functions from modiscloud package:

# Convert a byte integer (0-255) to a list of 8 bits
byteint2bit <- function(intval, reverse=TRUE)
{
  require(sfsmisc)		# for digitsBase
  
  if (reverse == TRUE)
  {
    byte_in_binary = rev(t(digitsBase(intval, base= 2, 8)))
  } else {
    byte_in_binary = c(t(digitsBase(intval, base= 2, 8)))
  }
  return(byte_in_binary)
}


# Get the value of a particular bit in a byte

extract_bit <- function(intval, bitnum)
{
  require(sfsmisc)		# for digitsBase
  
  # Convert to binary, read correctly (MODIS reads right-to-left)
  binval = rev(t(digitsBase(intval, base= 2, 8)))
  
  bitval = binval[bitnum]
  return(bitval)
}



## Implement the.....
d1 <- 0
for (k in 1:day_f){
  print(k)
  day_b_t <- substr(mod09ga_file[d1*8 + 1], 14, 16)
  v_day = FALSE
  
  ga <- 0
  gq <- 0
  m5 <- 0
  ga_v <- 0
  gq_v <- 0
  #m5_v <- 0
  #check m_ga band
  for (im in 1:8){
    c_day <- substr(mod09ga_file[d1*8 + (1 + ga_v)], 14, 16)
    #print(c_day)
    if(c_day == day_b_t){
      ga <- ga +1
    }
    ga_v <- ga_v + 1
  }
  #check m_gq band
  for (im in 1:2){
    c_day <- substr(mod09gq_file[d1*2 + (1 + gq_v)], 14, 16)
    #print(c_day)
    if(c_day == day_b_t){
      gq <- gq +1
    }
    gq_v <- gq_v + 1
  }
  
  Mod_35_ord <- list()
  f <- 1
  for (m_35 in mod35_file){
    day_mod_35 <- substr(m_35, 15, 17)
    #print(day_mod_35)
    if (day_mod_35 == day_b_t){
      Mod_35_ord[[f]] <- m_35
      f <- f + 1
      m5 <- m5 + 1
    }
  }
  
  if (ga==8 & gq==2 & m5==6){
    v_day = TRUE
  } else{
    v_day = FALSE
  }
  
  if (v_day == FALSE){
    
    ## Alta resolución GLS_ATA
    high_res_GLS_ATA <- stack(nan_band$band, nan_band$band, nan_band$band, nan_band$band, 
                              nan_band$band, nan_band$band, nan_band$band, Band_250_nan$elevation, Band_250_nan$slope, 
                              Band_250_nan$aspect, nan_band$band, nan_band$band) # con dem y Zenith, y cloud cover
    
    #high_res_GLS_ATA <- stack(exp(Band_250$band_1),exp(Band_250$band_2),p_b3_250_f_2, p_b4_250_f_2, 
    #p_b5_250_f_2, p_b6_250_f_2, p_b7_250_f_2,Band_250$Dem_250,Band_250$slope,Band_250$aspect,Band_250$zenith, Band_250$cloud_mask) # con dem y Zenith, y cloud cover
    names(high_res_GLS_ATA) <- c('band_1','band_2','band_3', 'band_4', 'band_5', 'band_6', 'band_7','elevation',
                                 'slope','aspect','zenith','cloud_mask')
    
    day <- day_b_t
    
    
  }else{
    #tic.clearlog()
    
    #tic(d1 + 1)
    #Sys.sleep(1)
    
    Zenith_500_d <- raster(paste0('Reflectance_bands/',year,'/MOD09GA/',mod09ga_file[d1*8 + 1]))
    Band_1_500_d <- raster(paste0('Reflectance_bands/',year,'/MOD09GA/',mod09ga_file[d1*8 + 2]))
    Band_2_500_d <- raster(paste0('Reflectance_bands/',year,'/MOD09GA/',mod09ga_file[d1*8 + 3]))
    Band_3_500_d <- raster(paste0('Reflectance_bands/',year,'/MOD09GA/',mod09ga_file[d1*8 + 4]))
    Band_4_500_d <- raster(paste0('Reflectance_bands/',year,'/MOD09GA/',mod09ga_file[d1*8 + 5]))
    Band_5_500_d <- raster(paste0('Reflectance_bands/',year,'/MOD09GA/',mod09ga_file[d1*8 + 6]))
    Band_6_500_d <- raster(paste0('Reflectance_bands/',year,'/MOD09GA/',mod09ga_file[d1*8 + 7]))
    Band_7_500_d <- raster(paste0('Reflectance_bands/',year,'/MOD09GA/',mod09ga_file[d1*8 + 8]))
    
    Band_1_250_d <- raster(paste0('Reflectance_bands/',year,'/MOD09GQ/',mod09gq_file[d1*2 + 1]))
    Band_2_250_d <- raster(paste0('Reflectance_bands/',year,'/MOD09GQ/',mod09gq_file[d1*2 + 2]))
    
    Band_cloud_1000_d <- raster(paste0('Reflectance_bands/',year,'/MOD35/',Mod_35_ord[[1]]))
    
    ## Aplicar la funcion
    #tic()
    #downs_func(Band_1_250_d,Band_2_250_d,Zenith_500_d,Band_1_500_d,Band_2_500_d,Band_3_500_d,Band_4_500_d,Band_5_500_d,
    #Band_6_500_d,Band_7_500_d,Band_water_res_250,Band_water_res_500,dem_250_f,dem_500_f,Band_cloud_1000_d, day_f, shape)
    #toc()
    
    print(day_b_t)
    day_year[d1+1] <- day_b_t
    
    ## Script
    
    # Setea nombre de los rasters
    Zenith_500 <- Zenith_500_d
    Band_1_500 <- Band_1_500_d
    Band_2_500 <- Band_2_500_d
    Band_3_500 <- Band_3_500_d
    Band_4_500 <- Band_4_500_d
    Band_5_500 <- Band_5_500_d
    Band_6_500 <- Band_6_500_d
    Band_7_500 <- Band_7_500_d
    
    Band_1_250 <- Band_1_250_d
    Band_2_250 <- Band_2_250_d
    
    Band_cloud_1000 <- Band_cloud_1000_d
    
    Band_water_1 <- Band_water_res_250
    Band_water_2 <- Band_water_res_500
    
    dem_1 <- dem_250_f
    dem_2 <- dem_500_f
    
    day <- day_b_t
    
    shape_pol <- shape
    
    ## Setear nombre de las bandas
    names(Band_1_250)[1] <- 'Band_1_250_r'
    names(Band_2_250)[1] <- 'Band_2_250_r'
    names(Zenith_500)[1] <- 'Zenith_250_r'
    names(Band_1_500)[1] <- 'Band_1_500_r'
    names(Band_2_500)[1] <- 'Band_2_500_r'
    names(Band_3_500)[1] <- 'Band_3_500_r'
    names(Band_4_500)[1] <- 'Band_4_500_r'
    names(Band_5_500)[1] <- 'Band_5_500_r'
    names(Band_6_500)[1] <- 'Band_6_500_r'
    names(Band_7_500)[1] <- 'Band_7_500_r'
    names(Band_water_1)[1] <- 'Water_250'
    names(Band_water_2)[1] <- 'Water_500'
    
    ## Setear valores NA
    Band_1_250_v <- clamp(Band_1_250, lower=-100, upper=16000, useValues=FALSE)
    Band_1_250_v[is.na(Band_1_250_v$Band_1_250_r)] <- NaN  # se cambio Na
    
    Band_2_250_v <- clamp(Band_2_250, lower=-100, upper=16000, useValues=FALSE)
    Band_2_250_v[is.na(Band_2_250_v$Band_2_250_r)] <- NaN  # se cambio Na
    
    #Zenith_500_v <- clamp(Zenith_500, lower=0, upper=18000, useValues=FALSE)
    #Zenith_500_v[is.na(Zenith_500_v$Zenith_250_r)] <- NaN  # se cambio Na
    
    Band_1_500_v <- clamp(Band_1_500, lower=-100, upper=16000, useValues=FALSE)
    Band_1_500_v[is.na(Band_1_500_v$Band_1_500_r)] <- NaN # se cambio Na
    
    Band_2_500_v <- clamp(Band_2_500, lower=-100, upper=16000, useValues=FALSE)
    Band_2_500_v[is.na(Band_2_500_v$Band_2_500_r)] <- NaN # se cambio Na
    
    Band_3_500_v <- clamp(Band_3_500, lower=-100, upper=16000, useValues=FALSE)
    Band_3_500_v[is.na(Band_3_500_v$Band_3_500_r)] <- NaN # se cambio Na
    
    Band_4_500_v <- clamp(Band_4_500, lower=-100, upper=16000, useValues=FALSE)
    Band_4_500_v[is.na(Band_4_500_v$Band_4_500_r)] <- NaN # se cambio Na
    
    Band_5_500_v <- clamp(Band_5_500, lower=-100, upper=16000, useValues=FALSE)
    Band_5_500_v[is.na(Band_5_500_v$Band_5_500_r)] <- NaN # se cambio Na
    
    Band_6_500_v <- clamp(Band_6_500, lower=-100, upper=16000, useValues=FALSE)
    Band_6_500_v[is.na(Band_6_500_v$Band_6_500_r)] <- NaN # se cambio Na
    
    Band_7_500_v <- clamp(Band_7_500, lower=-100, upper=16000, useValues=FALSE)
    Band_7_500_v[is.na(Band_7_500_v$Band_7_500_r)] <- NaN # se cambio Na
    
    # Ajusta resolución mascara de agua, dado los diferentes metodos de procesamiento
    Zenith_250 <- resample(Zenith_500, Band_1_250, method='bilinear') #promedio
    
    Cloud_cover_250 <- resample(Band_cloud_1000, Band_1_250, method='ngb')
    
    # Mascara de agua
    #Band_water[Band_water$Water == 1] <- NaN  # se cambio Na
    #Band_water[Band_water$Water == 0] <- 1
    
    # Ajusta resolución mascara de agua, dado los diferentes metodos de procesamiento
    #Band_water_res_250 <- resample(Band_water, Band_1_250, method='bilinear') #promedio
    #Band_water_res_500 <- resample(Band_water, Band_1_500, method='bilinear')
    
    # Ajusta resolución DEM
    #dem_250 <- resample(dem, Band_1_250, method='bilinear')
    #dem_500 <- resample(dem, Band_1_500, method='bilinear')
    
    # Mascara de NA
    Band_1_250_m <- Band_1_250_v
    values(Band_1_250_m)[values(Band_1_250_m) >= -100] = 1
    
    Band_2_250_m <- Band_2_250_v
    values(Band_2_250_m)[values(Band_2_250_m) >= -100] = 1
    
    Band_1_500_m <- Band_1_500_v
    values(Band_1_500_m)[values(Band_1_500_m) >= -100] = 1
    
    Band_2_500_m <- Band_2_500_v
    values(Band_2_500_m)[values(Band_2_500_m) >= -100] = 1
    
    Band_3_500_m <- Band_3_500_v
    values(Band_3_500_m)[values(Band_3_500_m) >= -100] = 1
    
    Band_3_500_m <- Band_3_500_v
    values(Band_3_500_m)[values(Band_3_500_m) >= -100] = 1
    
    Band_4_500_m <- Band_4_500_v
    values(Band_4_500_m)[values(Band_4_500_m) >= -100] = 1
    
    Band_5_500_m <- Band_5_500_v
    values(Band_5_500_m)[values(Band_5_500_m) >= -100] = 1
    
    Band_6_500_m <- Band_6_500_v
    values(Band_6_500_m)[values(Band_6_500_m) >= -100] = 1
    
    Band_7_500_m <- Band_7_500_v
    values(Band_7_500_m)[values(Band_7_500_m) >= -100] = 1
    
    mask_f_250 <- Band_water_1 * Band_1_250_m * Band_2_250_m  # Mascara 250 metros
    
    mask_f_500 <- Band_water_2 * Band_1_500_m * Band_2_500_m * Band_3_500_m * Band_4_500_m * 
      Band_5_500_m * Band_6_500_m * Band_7_500_m  # Mascara 500 metros
    
    ## Se aplica la mascara a los diferentes archivos
    Band_1_250_r <- Band_1_250 * mask_f_250
    Band_2_250_r <- Band_2_250 * mask_f_250
    Zenith_250_r <- Zenith_250 * mask_f_250
    Band_1_500_r <- Band_1_500 * mask_f_500
    Band_2_500_r <- Band_2_500 * mask_f_500
    Band_3_500_r <- Band_3_500 * mask_f_500
    Band_4_500_r <- Band_4_500 * mask_f_500
    Band_5_500_r <- Band_5_500 * mask_f_500
    Band_6_500_r <- Band_6_500 * mask_f_500
    Band_7_500_r <- Band_7_500 * mask_f_500
    
    dem_250_r <- dem_1 * mask_f_250
    names(dem_250_r) <- c('elevation','slope','aspect')
    dem_500_r <- dem_2 * mask_f_500
    names(dem_500_r) <- c('elevation','slope','aspect')
    
    
    
    ## Check band valid (no NaNs)
    b1_500_m <- cellStats(Band_1_500_r, 'mean')
    b2_500_m <- cellStats(Band_2_500_r, 'mean')
    b3_500_m <- cellStats(Band_3_500_r, 'mean')
    b4_500_m <- cellStats(Band_4_500_r, 'mean')
    b5_500_m <- cellStats(Band_5_500_r, 'mean')
    b6_500_m <- cellStats(Band_6_500_r, 'mean')
    b7_500_m <- cellStats(Band_7_500_r, 'mean')
    
    b1_500_sd <- cellStats(Band_1_500_r, 'sd')
    b2_500_sd <- cellStats(Band_2_500_r, 'sd')
    b3_500_sd <- cellStats(Band_3_500_r, 'sd')
    b4_500_sd <- cellStats(Band_4_500_r, 'sd')
    b5_500_sd <- cellStats(Band_5_500_r, 'sd')
    b6_500_sd <- cellStats(Band_6_500_r, 'sd')
    b7_500_sd <- cellStats(Band_7_500_r, 'sd')
    
    if (is.nan(b1_500_m)|is.nan(b2_500_m)|is.nan(b3_500_m)|is.nan(b4_500_m)|is.nan(b5_500_m)|
        is.nan(b6_500_m)|is.nan(b7_500_m)|b1_500_sd < 10|b2_500_sd < 10|b3_500_sd < 10|b4_500_sd < 10|
        b5_500_sd < 10|b6_500_sd < 10|b7_500_sd < 10){
      
      ## Alta resolución GLS_ATA
      high_res_GLS_ATA <- stack(nan_band$band, nan_band$band, nan_band$band, nan_band$band, 
                                nan_band$band, nan_band$band, nan_band$band, Band_250_nan$elevation, Band_250_nan$slope, 
                                Band_250_nan$aspect, nan_band$band, nan_band$band) # con dem y Zenith, y cloud cover
      
      #high_res_GLS_ATA <- stack(exp(Band_250$band_1),exp(Band_250$band_2),p_b3_250_f_2, p_b4_250_f_2, 
      #p_b5_250_f_2, p_b6_250_f_2, p_b7_250_f_2,Band_250$Dem_250,Band_250$slope,Band_250$aspect,Band_250$zenith, Band_250$cloud_mask) # con dem y Zenith, y cloud cover
      names(high_res_GLS_ATA) <- c('band_1','band_2','band_3', 'band_4', 'band_5', 'band_6', 'band_7','elevation',
                                   'slope','aspect','zenith','cloud_mask')
      #print(1)
      
      #band_250_s <- stack(Band_1_250_r, Band_2_250_r, dem_250_r, mask_f_250, Zenith_250_r, Cloud_cover_250)
      #names(band_250_s) <- c('band_1', 'band_2','Dem_250','slope','aspect', 'Mask', 'zenith','cloud')
      
      #band_500_s <- stack(Band_1_500_r, Band_2_500_r, Band_3_500_r, Band_4_500_r, Band_5_500_r, Band_6_500_r,
                          #Band_7_500_r, dem_500_r, mask_f_500)
      #names(band_500_s) <- c('band_1', 'band_2', 'band_3', 'band_4', 'band_5', 'band_6', 'band_7','Dem_500','slope','aspect', 'Mask')
      
      #reproject to UTM
      #band_250_s_UTM <- projectRaster(band_250_s, res=250, crs=CRS("+init=epsg:32719"), method = 'ngb')
      #band_500_s_UTM <- projectRaster(band_500_s, res=500, crs=CRS("+init=epsg:32719"), method = 'ngb')
      
      # data
      #R.500 <- crop(band_500_s_UTM, extent(shape_pol))
      #Band_500 <- mask(R.500, shape_pol)
      
      #R.250 <- crop(band_250_s_UTM, extent(shape_pol))
      #Band_250 <- mask(R.250, shape_pol)
      
      #nan_b1 <- Band_250$band_1 * NaN
      
      #high_res_GLS_ATA <- stack(nan_b1, nan_b1, nan_b1, nan_b1, 
                                #nan_b1, nan_b1, nan_b1, Band_250$Dem_250, Band_250$slope, 
                                #Band_250$aspect, Band_250$zenith, nan_b1) # con dem y Zenith, y cloud cover
      
      #high_res_GLS_ATA <- stack(exp(Band_250$band_1),exp(Band_250$band_2),p_b3_250_f_2, p_b4_250_f_2, 
      #p_b5_250_f_2, p_b6_250_f_2, p_b7_250_f_2,Band_250$Dem_250,Band_250$slope,Band_250$aspect,Band_250$zenith, Band_250$cloud_mask) # con dem y Zenith, y cloud cover
      #names(high_res_GLS_ATA) <- c('band_1','band_2','band_3', 'band_4', 'band_5', 'band_6', 'band_7','elevation',
                                   #'slope','aspect','zenith','cloud_mask')
    } else{
      ## Aumentar 200 para no tener problemas de falsas convergencias
      Band_1_250_r <- Band_1_250_r + 200
      Band_2_250_r <- Band_1_250_r + 200
      
      Band_1_500_r <- Band_1_500_r + 200
      Band_2_500_r <- Band_2_500_r + 200
      #hist(Band_2_500_r)
      Band_3_500_r <- Band_3_500_r + 200
      #hist(Band_3_500_r)
      Band_4_500_r <- Band_4_500_r + 200
      #hist(Band_3_500_r)
      Band_5_500_r <- Band_5_500_r + 200
      #hist(Band_3_500_r)
      Band_6_500_r <- Band_6_500_r + 200
      #hist(Band_3_500_r)
      Band_7_500_r <- Band_7_500_r + 200
      
      ## Normalizar por log
      Band_1_250_rl <- log(Band_1_250_r)
      Band_2_250_rl <- log(Band_2_250_r)
      
      Band_1_500_rl <- log(Band_1_500_r)
      Band_2_500_rl <- log(Band_2_500_r)
      Band_3_500_rl <- log(Band_3_500_r)
      Band_4_500_rl <- log(Band_4_500_r)
      Band_5_500_rl <- log(Band_5_500_r)
      Band_6_500_rl <- log(Band_6_500_r)
      Band_7_500_rl <- log(Band_7_500_r)
      
      # Ordena los archivos para la regresion
      band_250_s <- stack(Band_1_250_rl, Band_2_250_rl, dem_250_r, mask_f_250, Zenith_250_r, Cloud_cover_250)
      names(band_250_s) <- c('band_1', 'band_2','Dem_250','slope','aspect', 'Mask', 'zenith','cloud')
      
      band_500_s <- stack(Band_1_500_rl, Band_2_500_rl, Band_3_500_rl, Band_4_500_rl, Band_5_500_rl, Band_6_500_rl,
                          Band_7_500_rl, dem_500_r, mask_f_500)
      names(band_500_s) <- c('band_1', 'band_2', 'band_3', 'band_4', 'band_5', 'band_6', 'band_7','Dem_500','slope','aspect', 'Mask')
      
      #reproject to UTM
      band_250_s_UTM <- projectRaster(band_250_s, res=250, crs = crs_1, method = 'ngb') # was changed from crs=CRS("EPSG:32719")
      band_500_s_UTM <- projectRaster(band_500_s, res=500, crs = crs_1, method = 'ngb')
      
      # data
      R.500 <- crop(rast(band_500_s_UTM), ext(shape_pol))
      Band_500 <- mask(R.500, shape_pol)
      Band_500 <- brick(Band_500)
      
      R.250 <- crop(rast(band_250_s_UTM), ext(shape_pol)) # Extend by ext, and add rast function
      Band_250 <- mask(R.250, shape_pol)
      Band_250 <- brick(Band_250)
      
      Band_500$Dem_500_Z <- Z_r(Band_500$Dem_500)
      Band_250$Dem_500_Z <- Z_r(Band_250$Dem_250)
      
      data_500.p <- as(Band_500, 'SpatialPointsDataFrame')
      data_500.p.frame <- as.data.frame(na.omit(data_500.p))
      data_250.p <- as(Band_250, 'SpatialPointsDataFrame')
      data_250.p.frame <- as.data.frame(na.omit(data_250.p))
      
      ## GLS models
      fb3.2 = formula(band_3 ~ band_1)
      fb4 = formula(band_4 ~ band_1)
      fb5.2 = formula(band_5 ~ band_1 + Dem_500_Z)
      fb6 = formula(band_6 ~ band_1 + Dem_500_Z)
      fb7 = formula(band_7 ~ band_1 + Dem_500_Z)
      
      form_lis <- list(fb3.2, fb4, fb5.2, fb6, fb7) 
      
      # Linear regresion
      M3_lm <- gls(fb3.2, data = data_500.p)
      M4_lm <- gls(fb4, data = data_500.p)
      M5_lm <- gls(fb5.2, data = data_500.p)
      M6_lm <- gls(fb6, data = data_500.p)
      M7_lm <- gls(fb7, data = data_500.p)
      
      lm_m_list <- list(M3_lm, M4_lm, M5_lm, M6_lm, M7_lm)
      
      # variance structure
      vf1 <- varFixed(~ band_1)
      vf2 <- varIdent(form = ~band_1)
      
      vf3 <- varPower(form = ~band_1)
      vf4 <- varExp(form = ~band_1)
      vf5 <- varConstPower(form = ~band_1)
      
      var_list <- list(vf1, vf2, vf3, vf4, vf5)
      
      # models selection
      
      m <- 1
      gls_mod_f <- list()
      for (b in form_lis){
        list_gls <- list()
        list_gls[[1]] <- lm_m_list[[m]]
        l <- 2
        for (vl in var_list){
          check_gls <- tryCatch( expr = { gls(b,  weights = vl, data = data_500.p)
            TRUE
          },
          error = function(e){
            FALSE
          }
          )
          if (check_gls == FALSE){
          } else {
            M_gls <- gls(b,  weights = vl, data = data_500.p)
            #gls_b6_m[2,2] = TRUE
            list_gls[[l]] <- M_gls
            l <- l + 1
          }
          #print(check_gls.5)
        }
        min_r <- 1
        pos_min <- 1
        min_gls <- AIC(list_gls[[1]])
        for (g in list_gls){
          if (AIC(g) < min_gls){
            min_gls <- AIC(g)
            pos_min <- min_r
          }
          min_r <- min_r + 1
        }
        gls_mod_f[[m]] <- list_gls[[pos_min]]
        m <- m+1
      }
      
      #band3
      M3.2_gls.3 <- gls_mod_f[[1]]
      
      #band4
      M4_gls.4 <- gls_mod_f[[2]]
      
      #band5
      M5_gls.5 <- gls_mod_f[[3]]
      
      #band_6
      M6_gls.5 <- gls_mod_f[[4]]
      
      #band_7
      M7_gls.5 <- gls_mod_f[[5]]
      
      #residuos at 500 meter resolution
      p_l_2 <- predict(Band_500, M3.2_gls.3)
      p_2 <- exp(p_l_2)
      resid_b3_2 <- exp(Band_500$band_3) - p_2
      
      p_l_b4_2 <- predict(Band_500, M4_gls.4)
      p_b4_2 <- exp(p_l_b4_2)
      resid_b4_2 <- exp(Band_500$band_4) - p_b4_2
      
      p_l_b5_2 <- predict(Band_500, M5_gls.5)
      p_b5_2 <- exp(p_l_b5_2)
      resid_b5_2 <- exp(Band_500$band_5) - p_b5_2
      
      p_l_b6_2 <- predict(Band_500, M6_gls.5)
      p_b6_2 <- exp(p_l_b6_2)
      resid_b6_2 <- exp(Band_500$band_6) - p_b6_2
      
      p_l_b7_2 <- predict(Band_500, M7_gls.5)
      p_b7_2 <- exp(p_l_b7_2)
      resid_b7_2 <- exp(Band_500$band_7) - p_b7_2
      
      vario <- c('Exp', 'Gau', 'Mat', 'Sph', 'Bes', 'Exc', 'Cir', 'Pen', 'Log', 'Nug', 'Lin', 'Wav', 'Spl')
      
      ## Variograma para los residuos
      Res_b3.500_point_2 <- as(resid_b3_2, 'SpatialPointsDataFrame')
      V.r3_2 <- variogram(layer ~1, data =Res_b3.500_point_2)
      
      # Selección de variograma b3
      SS_Err_1 <- 100000000000
      for (v in vario){
        v.fit <- fit.variogram(V.r3_2, vgm(v))
        SS_Err_2 <- attr(v.fit,"SSErr")
        #print(v)
        #print(SS_Err_2)
        if (is.na(attr(v.fit,"SSErr")) != TRUE){
          if (SS_Err_2 < SS_Err_1){
            SS_Err_1 <- SS_Err_2
            v_sel_3 <- v
          }
        }
      }
      v.fit_b3 <- fit.variogram(V.r3_2, vgm(v_sel_3))
      
      Res_b4.500_point_2 <- as(resid_b4_2, 'SpatialPointsDataFrame')
      V.r4_2 <- variogram(layer ~1, data =Res_b4.500_point_2)
      # Selección de variograma b4
      SS_Err_1 <- 100000000000
      for (v in vario){
        v.fit <- fit.variogram(V.r4_2, vgm(v))
        SS_Err_2 <- attr(v.fit,"SSErr")
        #print(v)
        #print(SS_Err_2)
        if (is.na(attr(v.fit,"SSErr")) != TRUE){
          if (SS_Err_2 < SS_Err_1){
            SS_Err_1 <- SS_Err_2
            v_sel_4 <- v
          }
        }
      }
      v.fit_b4 <- fit.variogram(V.r4_2, vgm(v_sel_4))
      
      # Selección de variograma b5
      Res_b5.500_point_2 <- as(resid_b5_2, 'SpatialPointsDataFrame')
      V.r5_2 <- variogram(layer ~1, data =Res_b5.500_point_2)
      # Selección de variograma b5
      SS_Err_1 <- 100000000000
      for (v in vario){
        v.fit <- fit.variogram(V.r5_2, vgm(v))
        SS_Err_2 <- attr(v.fit,"SSErr")
        if (is.na(attr(v.fit,"SSErr")) != TRUE){
          if (SS_Err_2 < SS_Err_1){
            SS_Err_1 <- SS_Err_2
            v_sel_5 <- v
          }
        }
      }
      v.fit_b5 <- fit.variogram(V.r5_2, vgm(v_sel_5))
      
      Res_b6.500_point_2 <- as(resid_b6_2, 'SpatialPointsDataFrame')
      V.r6_2 <- variogram(layer ~1, data =Res_b6.500_point_2)
      # Selección de variograma b6
      SS_Err_1 <- 100000000000
      for (v in vario){
        v.fit <- fit.variogram(V.r6_2, vgm(v))
        SS_Err_2 <- attr(v.fit,"SSErr")
        if (is.na(attr(v.fit,"SSErr")) != TRUE){
          if (SS_Err_2 < SS_Err_1){
            SS_Err_1 <- SS_Err_2
            v_sel_6 <- v
          }
        }
      }
      v.fit_b6 <- fit.variogram(V.r6_2, vgm(v_sel_6))
      
      Res_b7.500_point_2 <- as(resid_b7_2, 'SpatialPointsDataFrame')
      V.r7_2 <- variogram(layer ~1, data =Res_b7.500_point_2)
      # Selección de variograma b7
      SS_Err_1 <- 100000000000
      for (v in vario){
        v.fit <- fit.variogram(V.r7_2, vgm(v))
        SS_Err_2 <- attr(v.fit,"SSErr")
        if (is.na(attr(v.fit,"SSErr")) != TRUE){
          if (SS_Err_2 < SS_Err_1){
            SS_Err_1 <- SS_Err_2
            v_sel_7 <- v
          }
        }
      }
      v.fit_b7 <- fit.variogram(V.r7_2, vgm(v_sel_7))
      
      ## Atakrig
      
      # raster format Rasterlayer (raster) to SpatRaster (Terra):
      #Band_250_sr <- as(Band_250$band_1, "SpatRaster") # is working!
      Band_250_sr <- rast(Band_250$band_1) # working
      grid.pred_250 <- discretizeRaster(Band_250_sr, 200, type = "all") # grilla para la predicción (raster)
      
      res_500_b3.d_2_sr <- rast(resid_b3_2$layer)
      res_500_b3.d_2 <- discretizeRaster(res_500_b3.d_2_sr, 200) # este es el que quiero interpolar
      
      res_500_b4.d_2_sr <- rast(resid_b4_2$layer)
      res_500_b4.d_2 <- discretizeRaster(res_500_b4.d_2_sr, 200)
      
      res_500_b5.d_2_sr <- rast(resid_b5_2$layer)
      res_500_b5.d_2 <- discretizeRaster(res_500_b5.d_2_sr, 200)
      
      res_500_b6.d_2_sr <- rast(resid_b6_2$layer)
      res_500_b6.d_2 <- discretizeRaster(res_500_b6.d_2_sr, 200)
      
      res_500_b7.d_2_sr <- rast(resid_b7_2$layer)
      res_500_b7.d_2 <- discretizeRaster(res_500_b7.d_2_sr, 200)
      
      ## Control de error: Error in { : 
      #task 1 failed - "dims [product 16] do not match the length of object [1]"
      
      v_fit_list <- list(v.fit_b3, v.fit_b4, v.fit_b5, v.fit_b6, v.fit_b7)
      res_500_list <- list(res_500_b3.d_2, res_500_b4.d_2, res_500_b5.d_2, res_500_b6.d_2, res_500_b7.d_2)
      band_e_lis <- list(Band_500$band_3, Band_500$band_4, Band_500$band_5, Band_500$band_6, Band_500$band_7)
      
      pred_list <- list()
      
      ataStartCluster()
      
      ak <- 1
      #check_ata_list <- list()
      #prueba <- ataKriging(res_500_list[ak], grid.pred_250, v_fit_list[ak], showProgress=T)
      for (f in v_fit_list){
        check_ata <- tryCatch( expr = { pred.ataok_prob <- ataKriging(res_500_list[[ak]], grid.pred_250, v_fit_list[[ak]], showProgress=T)
        check_ata_v <- TRUE
        pred.ataok <- pred.ataok_prob
        },
        error = function(e){
          check_ata_v <- FALSE
        }
        )
        print(check_ata_v)
        print(ak)
        #print(check_ata)
        if (check_ata_v == FALSE){
          band_e <- resample(exp(band_e_lis[[ak]]), Band_250$band_1, method='bilinear')
          band_e_f <- band_e * Band_250$Mask
          pred_list[[ak]] <- band_e_f
        } else {
          ## Generar el raster final
          #pred.ataok <- ataKriging(res_500_list[[ak]], grid.pred_250, v_fit_list[[ak]], showProgress=T)
          res.r_250_2 <- rasterFromXYZ(pred.ataok[,2:4], crs = crs_1) # was changed from crs = "+init=epsg:32719"
          res.r_250_f_2 <- res.r_250_2 * Band_250$Mask
          
          p_l_250_2 <- predict(Band_250, gls_mod_f[[ak]])
          p_b_250_2 <- exp(p_l_250_2)
          p_b_250_f_2 <- p_b_250_2 + res.r_250_f_2
          
          pred_list[[ak]] <- p_b_250_f_2
        }
        ak <- ak + 1
      }
      
      ataStopCluster()
      
      
      
      #pred.ataok_rb3_2 <- ataKriging(res_500_b3.d_2, grid.pred_250, v.fit_b3, showProgress=T)
      #pred.ataok_rb4_2 <- ataKriging(res_500_b4.d_2, grid.pred_250, v.fit_b4, showProgress=T)
      #pred.ataok_rb5_2 <- ataKriging(res_500_b5.d_2, grid.pred_250, v.fit_b5, showProgress=T)
      #pred.ataok_rb6_2 <- ataKriging(res_500_b6.d_2, grid.pred_250, v.fit_b6, showProgress=T)
      #pred.ataok_rb7_2 <- ataKriging(res_500_b7.d_2, grid.pred_250, v.fit_b7, showProgress=T)
      #ataStopCluster()
      
      #res_b3.r_250_2 <- rasterFromXYZ(pred.ataok_rb3_2[,2:4], crs = "+init=epsg:32719")
      #res_b4.r_250_2 <- rasterFromXYZ(pred.ataok_rb4_2[,2:4], crs = "+init=epsg:32719")
      #res_b5.r_250_2 <- rasterFromXYZ(pred.ataok_rb5_2[,2:4], crs = "+init=epsg:32719")
      #res_b6.r_250_2 <- rasterFromXYZ(pred.ataok_rb6_2[,2:4], crs = "+init=epsg:32719")
      #res_b7.r_250_2 <- rasterFromXYZ(pred.ataok_rb7_2[,2:4], crs = "+init=epsg:32719")
      
      #res_b3.r_250_f_2 <- res_b3.r_250_2 * Band_250$Mask
      #res_b4.r_250_f_2 <- res_b4.r_250_2 * Band_250$Mask
      #res_b5.r_250_f_2 <- res_b5.r_250_2 * Band_250$Mask
      #res_b6.r_250_f_2 <- res_b6.r_250_2 * Band_250$Mask
      #res_b7.r_250_f_2 <- res_b7.r_250_2 * Band_250$Mask
      
      ## Aplicar modelo a 250 GLS
      #p_l_250_2 <- predict(Band_250, M3.2_gls.3)
      #p_b3_250_2 <- exp(p_l_250_2)
      #p_b3_250_f_2 <- p_b3_250_2 + res_b3.r_250_f_2
      
      #p_l_250_b4_2 <- predict(Band_250, M4_gls.4)
      #p_b4_250_2 <- exp(p_l_250_b4_2)
      #p_b4_250_f_2 <- p_b4_250_2 + res_b4.r_250_f_2
      
      #p_l_250_b5_2 <- predict(Band_250, M5_gls.5)
      #p_b5_250_2 <- exp(p_l_250_b5_2)
      #p_b5_250_f_2 <- p_b5_250_2 + res_b5.r_250_f_2
      
      #p_l_250_b6_2 <- predict(Band_250, M6_gls.5)
      #p_b6_250_2 <- exp(p_l_250_b6_2)
      #p_b6_250_f_2 <- p_b6_250_2 + res_b6.r_250_f_2
      
      #p_l_250_b7_2 <- predict(Band_250, M7_gls.5)
      #p_b7_250_2 <- exp(p_l_250_b7_2)
      #p_b7_250_f_2 <- p_b7_250_2 + res_b7.r_250_f_2
      
      ## Cloud mask
      cloud_mask.p <- as(Band_250$cloud, 'SpatialPointsDataFrame')
      cloud_mask.p.frame <- as.data.frame(cloud_mask.p)
      colnames(cloud_mask.p.frame) <- c('value','x' ,'y')
      
      #extrae bit 0
      CM_bit_0 = function(n){
        extract_bit(n, bitnum=1)
      }
      
      #extrae bit 1
      CM_bit_1 = function(n){
        extract_bit(n, bitnum=2)
      }
      
      #extrae bit 2
      CM_bit_2 = function(n){
        extract_bit(n, bitnum=3)
      }
      
      for (i in 1:nrow(cloud_mask.p.frame)){
        cloud_mask.p.frame[i,4] <- CM_bit_0(cloud_mask.p.frame[i,1])
        cloud_mask.p.frame[i,5] <- CM_bit_1(cloud_mask.p.frame[i,1])
        cloud_mask.p.frame[i,6] <- CM_bit_2(cloud_mask.p.frame[i,1])
        #print(CM_bit_2(cloud_mask.p.frame[i,1]))
      }
      for (i in 1:nrow(cloud_mask.p.frame)){
        if(cloud_mask.p.frame[i,4] == 1 & cloud_mask.p.frame[i,5] == 0 & cloud_mask.p.frame[i,6] == 0){
          cloud_mask.p.frame[i,7]<- 0
        } else cloud_mask.p.frame[i,7]<- 1
      }
      #arma raster cloud cover
      coordinates(cloud_mask.p.frame) <- ~x+y
      proj4string(cloud_mask.p.frame) <- crs_1 # UTM huso 19S WGS84 # was changed from CRS("+init=epsg:32719")
      cloud_mask_r <- rasterFromXYZ(cloud_mask.p.frame[,5])
      
      Band_250$cloud_mask <- resample(cloud_mask_r, Band_250$band_1, method='ngb')
      
      ## Alta resolución GLS_ATA
      high_res_GLS_ATA <- stack(exp(Band_250$band_1),exp(Band_250$band_2),pred_list[[1]], pred_list[[2]], 
                                pred_list[[3]], pred_list[[4]], pred_list[[5]], Band_250$Dem_250, Band_250$slope, 
                                Band_250$aspect, Band_250$zenith, Band_250$cloud_mask) # con dem y Zenith, y cloud cover
      
      #high_res_GLS_ATA <- stack(exp(Band_250$band_1),exp(Band_250$band_2),p_b3_250_f_2, p_b4_250_f_2, 
      #p_b5_250_f_2, p_b6_250_f_2, p_b7_250_f_2,Band_250$Dem_250,Band_250$slope,Band_250$aspect,Band_250$zenith, Band_250$cloud_mask) # con dem y Zenith, y cloud cover
      names(high_res_GLS_ATA) <- c('band_1','band_2','band_3', 'band_4', 'band_5', 'band_6', 'band_7','elevation',
                                   'slope','aspect','zenith','cloud_mask')
    }
    
  }
  
  # 7 bands at 250 meters (incluye:dem(elevation, slope,aspect),zenith, cloud cover)
  writeRaster(high_res_GLS_ATA, paste0(Out_link,'/',year,'_',day,'_pred_GLS-ATAK.tif'), 
              format="GTiff", overwrite=TRUE)
  #toc()
  #day_time
  
  #plot(high_res_GLS_ATA)
  
  day_year[d1+1] <- paste0(year,'_',day,'_pred_GLS-ATAK.tif')
  #day_time[d1+1] <- toc()
  
  
  #toc(log = TRUE, quiet = TRUE)
  
  d1 <- d1 + 1
  
}

#log.txt <- tic.log(format = TRUE)
#log.lst <- tic.log(format = FALSE)
#tic.clearlog()
#timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))

# Escribe archivo

#day_year <- c('2015_217_pred_GLS-ATAK.tif','2015_268_pred_GLS-ATAK.tif')
#timings <- c('305', '306')

fileConn<-file(paste0(Out_link,'/',year,'_downcaling.txt'))
writeLines(as.character(day_year), fileConn)
close(fileConn)

#ata_1 <- brick(paste0(Out_link,'/',year,'_',209,'_pred_GLS-ATAK.tif'))

#ata_2 <- brick(paste0(Out_link,'/',year,'_',210,'_pred_GLS-ATAK.tif'))
#ata_3 <-stack(ata_1,ata_2)



