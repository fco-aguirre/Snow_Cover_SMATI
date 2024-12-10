## Revisar formula Albedo!!

## Calculo de albedo

# Values to broadband albedo Solar 
Sw_A_30 <- 0.0765
Sw_B_30 <- 0.2205
pend_Sw_A <- -3.9*10^(-4)
pend_Sw_B <- 1.77*10^(-4)

## revisar limite de snow fraction para realizar el calculo del albedo

for(a in 1:nrow(Snow_st.p.frame)){
  if(is.nan(Snow_st.p.frame[a,9])){
    Snow_st.p.frame[a,11] <- NaN
  }else{
    if(Snow_st.p.frame[a,9] >= 0.01){
      Sw_A_ang <- pend_Sw_A*(Snow_st.p.frame[a,6]-30) + Sw_A_30
      Sw_B_ang <- pend_Sw_B*(Snow_st.p.frame[a,6]-30) + Sw_B_30
      
      Snow_st.p.frame[a,11] <- 1 - (Sw_A_ang * (Snow_st.p.frame[a,2]^Sw_B_ang))
    }
    else{
      Snow_st.p.frame[a,11] <- NaN
    }
  }
}

# Arma el raster
coordinates(Snow_st.p.frame) <- ~x+y
proj4string(Snow_st.p.frame) <- crs_1 # UTM huso 19S WGS84 !new format!

Snow_f_r <- rasterFromXYZ(Snow_st.p.frame[,1])
Snow_gz_r <- rasterFromXYZ(Snow_st.p.frame[,2])
NDSI_r <- rasterFromXYZ(Snow_st.p.frame[,3])
MADI_r <- rasterFromXYZ(Snow_st.p.frame[,4])
Snow_M <- rasterFromXYZ(Snow_st.p.frame[,7])
Snow_albedo <- rasterFromXYZ(Snow_st.p.frame[,9])

Snow_umb_r <- stack(Snow_f_r, Snow_gz_r, NDSI_r, MADI_r, Snow_M, Snow_albedo)
crs(Snow_umb_r) <- crs_1
plot(Snow_umb_r)

# Print sub_product_2
writeRaster(Snow_umb_r, paste0(Out_link,'/',day_b_d,'_mesma_albedo_p2.tif'), overwrite=TRUE)
