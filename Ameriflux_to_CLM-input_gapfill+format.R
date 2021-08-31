library(REddyProc)
library(dplyr)
library(ggplot2)
library(lubridate)
library(ncdf4)
library(gridExtra)

# tell R not to use scientific notation
options(scipen = 999)

# bring in ameriflux data
rcls <- read.csv("C:/Users/Derek/Google Drive/RCrk/RC Flux Tower data/RC-Rls/from_Ameriflux/AMF_US-Rls_BASE_HH_3-5.csv", as.is = T)
flux_file_name <- "AMF_US-Rls_BASE_HH_3-5"

# Change -9999 to NA
rcls[rcls == -9999] <- NA

# Remove empty rows at beginning for data
rcls <- rcls[12125:nrow(rcls),]

# Fix time format
rcls$TIMESTAMP_START <- as.POSIXct(as.character(rcls$TIMESTAMP_START),format="%Y%m%d%H%M", tz="Etc/GMT-7")
rcls$TIMESTAMP_END <- as.POSIXct(as.character(rcls$TIMESTAMP_END),format="%Y%m%d%H%M", tz="Etc/GMT-7")
rcls$DateTime <- rcls$TIMESTAMP_END

#Start at first full day
#End at last full day
set_first_day = as.POSIXct("201409110030",format="%Y%m%d%H%M", tz="Etc/GMT-7") 
set_last_day =  as.POSIXct("201809300000",format="%Y%m%d%H%M", tz="Etc/GMT-7")  
rcls <- rcls %>% filter(between(DateTime, set_first_day,set_last_day))

##############################################################################################
### Plot to check sunrise time
sunrise_plt_strt <- as.POSIXct("201604010000",format="%Y%m%d%H%M", tz="Etc/GMT-7") 
sunrise_plt_end <- as.POSIXct("201604040000",format="%Y%m%d%H%M", tz="Etc/GMT-7")  

sunrise_plt_df <- rcls
sunrise_plt_df$DateTime <- sunrise_plt_df$DateTime# + hours(7)
sunrise_plt_df <- sunrise_plt_df %>% filter(between(DateTime, sunrise_plt_strt, sunrise_plt_end))

time_plt_raw <- ggplot(sunrise_plt_df, aes(x=DateTime, y=SW_IN)) + geom_point() + geom_line() +ggtitle("Rls raw - 2016") +
  xlim(sunrise_plt_strt, sunrise_plt_end)
time_plt_raw 

time_plt_up7 <- ggplot(sunrise_plt_df, aes(x=DateTime+hours(7), y=SW_IN)) + geom_point() + geom_line() + ggtitle("Rls shift to UTC - 2016") +
  xlim(sunrise_plt_strt, sunrise_plt_end)
time_plt_up7
#############################################################################################


# create year column
rcls$Year <- as.numeric(strftime(rcls$DateTime, format = "%Y"))

# create day of year column
rcls$DoY <- as.numeric(strftime(rcls$DateTime, format = "%j"))

# create decimal time column
rcls$Hour <- as.numeric(strftime(rcls$DateTime, format = "%H")) + as.numeric(strftime(rcls$DateTime, format = "%M"))/60

#calc net radiation
rcls$radNet = as.numeric(as.character((rcls$SW_IN - rcls$SW_OUT) + (rcls$LW_IN - rcls$LW_OUT)))

# build required flux data dataframe
rcls_flux_df <- rcls %>% select(DateTime, Year, DoY, Hour, NEE_PI, LE, H, TA, TS_PI_1, RH, USTAR, P_PI_F, PA, WS, SW_IN, LW_IN, radNet) 

# fix bound on RH, can't be over 100%
rcls_flux_df <- rcls_flux_df%>% mutate(RH=replace(RH, RH > 100, 100))

# calc vapor pressure deficit
rcls_flux_df$VPD = fCalcVPDfromRHandTair(rH = rcls_flux_df$RH, Tair = rcls_flux_df$TA)
  #tail(sort(unlist(rcls_flux_df$VPD, use.names = FALSE)), 7) #check highest VPD values

# calc PRECTmms 
rcls_flux_df$PRECTmms = rcls_flux_df$P_PI_F / (30*60) # mm/s 

# set column names
colnames(rcls_flux_df) <- c("DateTime", "Year", "DoY", "Hour", "NEE", "LE", "H", "Tair", "Tsoil", "rH", "Ustar", "P", "PA", "U", "Rg", "FLDS", "radNet", "VPD", "PRECTmms")

#check time intervals
# for(i in 2:nrow(rcls_flux_df)) {
#   print(i)
#   interval <- as.numeric(rcls_flux_df$DateTime[i] - rcls_flux_df$DateTime[i-1])
#   if(interval != 30) {
#     print(interval)
#     print(i)
#     break
#   }
# }

EddyProc.C <- sEddyProc$new('Rcls', rcls_flux_df, c("NEE", "LE", "H", "Tair", "Tsoil", "rH", "Ustar", "P", "PA", "U", "Rg", "FLDS", "radNet", "VPD", "PRECTmms"))

# gap filling (without prior Ustar filtering)
# +++ Fill gaps in variables with MDS gap filling algorithm (without prior ustar filtering)
#  Note, this also takes a long time to complete!

EddyProc.C$sMDSGapFill('NEE', FillAll = TRUE) #Fill all values to estimate flux uncertainties
EddyProc.C$sMDSGapFill('LE', FillAll = TRUE)
EddyProc.C$sMDSGapFill('H', FillAll = TRUE)
EddyProc.C$sMDSGapFill('Ustar', FillAll = TRUE)
EddyProc.C$sMDSGapFill('Tair', FillAll = TRUE) #DP switched from FALSE to TRUE  
EddyProc.C$sMDSGapFill('VPD', FillAll = FALSE) 
EddyProc.C$sMDSGapFill('rH', FillAll = TRUE) #DP switched from FALSE to TRUE  
EddyProc.C$sMDSGapFill('U', FillAll = FALSE) # wind
EddyProc.C$sMDSGapFill('PRECTmms', FillAll = TRUE) #DP switched from FALSE to TRUE  
EddyProc.C$sMDSGapFill('P', FillAll = FALSE) 
EddyProc.C$sMDSGapFill('FLDS', FillAll = TRUE) #DP switched from FALSE to TRUE  
EddyProc.C$sMDSGapFill('Rg', FillAll = TRUE) #DP switched from FALSE to TRUE  
EddyProc.C$sMDSGapFill('radNet', FillAll = FALSE) 
EddyProc.C$sMDSGapFill('Tsoil', FillAll = FALSE)  
EddyProc.C$sMDSGapFill('PA', FillAll = TRUE) #DP switched from FALSE to TRUE  


# set location of site
#latSite = 43.1439 #43.1439 for Rls
latSite = 42.8795811518324

#set lonSite for clm
#lonSiteE = -116.7356 + 360
lonSiteE = 243.75

#lonSiteW = -116.7356 #-116.7356 for Rls, +360 to convert from W to E
lonSiteW = lonSiteE - 360

EddyProc.C$sSetLocationInfo(LatDeg = latSite, LongDeg = lonSiteW, TimeZoneHour = -7)

# Nighttime-based partitioning of measured net ecosystem fluxes into gross primary production (GPP) and ecosystem respiration (Reco)
EddyProc.C$sMRFluxPartition()

#+++ Export gap filled and partitioned data to standard data frame
FilledEddyData.F <- EddyProc.C$sExportResults()

# Grab just the filled data products
dataClm <- FilledEddyData.F[,grep(pattern = "_f$", x = names(FilledEddyData.F))]

# Grab the POSIX timestamp
dataClm$DateTime <- rcls_flux_df$DateTime #- lubridate::minutes(30) # putting back to original position

# fix names
names(dataClm) <- gsub("_f", "", names(dataClm))

# A fingerprint-plot of the source half-hourly shows already several gaps. 
# A fingerprint-plot is a color-coded image of the half-hourly fluxes by daytime on the x and and day of the year on the y axis.
#EddyProc.C$sPlotFingerprintY('PA', Year = 2017)

# Convert degC to K for temperature
dataClm$Tair <- dataClm$Tair + 273.15
attributes(obj = dataClm$Tair)$units <- "K"

# Convert kPa to Pa for pressure
dataClm$PA <- dataClm$PA * 1000.0
attributes(obj = dataClm$PA)$units <- "Pa"

# Create tower height measurement field
dataClm$ZBOT <- rep(2.09,nrow(dataClm))  #RCrk systems were mounted to towers between
                                         #1.7 and 2.5 m above the plant canopy; heights
                                         #were 2.05, 2.09, 3.5, and 2.5 m above the ground
                                         #surface for the WBS, LoS, PFS, and MBS site,
                                         #respectively.




##########################################################################################
# Check for NA in dataframe
for(i in 1:(ncol(dataClm)-3)){
  if(i==1){
    print("Count of NA values by column:")
  }
  print(paste0(colnames(dataClm)[i], ": ", sum(is.na(dataClm[,i])))) 
}

#QA plotter
plot.df <- dataClm 

### Debug cache dataCLM in dataClm2
dataClm2 <- dataClm
#dataClm <- dataClm2

#ggplot(plot.df, aes(x=DateTime, y=PA)) + geom_line()

#from flux tower data
#ggplot(rcls_flux_df, aes(x=DateTime, y=PA)) + geom_line()

#########################################################################################
# Select only the continuous data set (no gaps)
###############################################################
# Convert time from MT to GMT == 7 hours forward
#dataClm$DateTime[1]
dataClm$DateTime <- dataClm$DateTime + hours(7)
#dataClm$DateTime[1]

#set period to select from data
set_start_date = as.POSIXct("201502190000",format="%Y%m%d%H%M", tz="Etc/GMT-7") 
set_end_date =  as.POSIXct("201809010000",format="%Y%m%d%H%M", tz="Etc/GMT-7")  
dataClm <- dataClm %>% filter(between(DateTime, set_start_date,set_end_date))

#remove leap year days (Feb 29th)
dataClm <- dataClm[!(format(dataClm$DateTime,"%m") == "02" & format(dataClm$DateTime, "%d") == "29"), ,drop = FALSE]

# Year month combination for data filtering
dataClm$yearMon <- paste0(year(dataClm$DateTime), "-", 
                          sprintf("%02d", month(dataClm$DateTime)))

# Quick check of data for gaps
#ggplot(dataClm, aes(x=DateTime, y=PA)) + geom_line()


##############################################################################
# QC output
##############################################################################

#str(dataClm$DateTime)
#str(dataClm2$DateTime)

setYearMon <- unique(dataClm$yearMon)

for (j in setYearMon) {

  dataClm.mon <- dataClm[dataClm$yearMon == j,]
  
  print(paste0(j, ", nrow: ", as.character(nrow(dataClm.mon)), ", start: ", as.character(min(dataClm.mon$DateTime))))
}


#as.Date(dataClm$DateTime[1],format="%Y%m%d") 


##############################################################################
# Write output to CLM
##############################################################################
DirOut <- "C:/Users/Derek/Google Drive/RCrk/RC Flux Tower data/RC-Rls/gap-filled_monthly_netCDF"
tower_name <- "rcls"

#Define missing value fill
mv <- -9999  

#Set of year/month combinations for netCDF output
setYearMon <- unique(dataClm$yearMon)


# DEBUG: Check month start days 
#############################################
# for (m in setYearMon) {
#   print(m)
#   Data.mon <- dataClm[dataClm$yearMon == m,]
#   print(paste("days since",as.Date(Data.mon$DateTime[1],format="%Y%m%d"), "00:00:00"))
# }
# 
# df_check <- dataClm[dataClm$yearMon == "2016-03",]
# df_check$DateTime[1]
# as.Date(df_check$DateTime[1],format="%Y%m%d", tz = "Etc/GMT-7")
# print(paste("days since",as.Date(df_check$DateTime[1],format="%Y%m%d", tz = "Etc/GMT-7"), "00:00:00"))

#################################################

verbose = TRUE
for (m in setYearMon) {
  
  #DEBUG
  #m <- setYearMon[10] #for testing
  
  Data.mon <- dataClm[dataClm$yearMon == m,]
  timeStep <- seq(0,nrow(Data.mon)-1,1)
  time     <- timeStep/48
  #endStep  <- startStep + nsteps[m]-1
  
  if (verbose) {
    print(paste(m,"Data date =", as.Date(Data.mon$DateTime[1],format="%Y%m%d", tz = "Etc/GMT-7"), "00:00:00"))
    names(Data.mon)
  }
  
  #NetCDF output filename
  fileOutNcdf <- paste(DirOut,"/",m,".nc", sep = "")
  if (verbose) {
    print(fileOutNcdf)
  }
  
  # define the netcdf coordinate variables (name, units, type)
  lat  <- ncdf4::ncdim_def("lat","degrees_north", as.double(latSite), create_dimvar=TRUE)
  lon <- ncdf4::ncdim_def("lon","degrees_east", as.double(lonSiteE), create_dimvar=TRUE)
  
  #Variables to output to netCDF
  time <- ncdf4::ncdim_def("time", paste("days since",as.Date(Data.mon$DateTime[1],format="%Y%m%d", tz = "Etc/GMT-7"), "00:00:00"),
                           vals=as.double(time),unlim=FALSE, create_dimvar=TRUE,
                           calendar = "gregorian")
  LATIXY  <- ncdf4::ncvar_def("LATIXY", "degrees N", list(lat), mv,
                              longname="latitude", prec="double")
  LONGXY  <- ncdf4::ncvar_def("LONGXY", "degrees E", list(lon), mv,
                              longname="longitude", prec="double")
  FLDS  <- ncdf4::ncvar_def("FLDS", "W/m^2", list(lon,lat,time), mv,
                            longname="incident longwave (FLDS)", prec="double")
  FSDS  <- ncdf4::ncvar_def("FSDS", "W/m^2", list(lon,lat,time), mv,
                            longname="incident shortwave (FSDS)", prec="double")
  PRECTmms <- ncdf4::ncvar_def("PRECTmms", "mm/s", list(lon,lat,time), mv,
                               longname="precipitation (PRECTmms)", prec="double")
  PSRF  <- ncdf4::ncvar_def("PSRF", "Pa", list(lon,lat,time), mv,
                            longname="pressure at the lowest atmospheric level (PSRF)", prec="double")
  RH    <- ncdf4::ncvar_def("RH", "%", list(lon,lat,time), mv,
                            longname="relative humidity at lowest atm level (RH)", prec="double")
  TBOT  <- ncdf4::ncvar_def("TBOT", "K", list(lon,lat,time), mv,
                            longname="temperature at lowest atm level (TBOT)", prec="double")
  WIND  <- ncdf4::ncvar_def("WIND", "m/s", list(lon,lat,time), mv,
                            longname="wind at lowest atm level (WIND)", prec="double")
  ZBOT  <- ncdf4::ncvar_def("ZBOT", "m", list(lon,lat,time), mv,
                            longname="observational height", prec="double")
  NEE <- ncdf4::ncvar_def("NEE", "umolm-2s-1", list(lon,lat,time), mv,
                          longname="net ecosystem exchange", prec="double")
  FSH  <- ncdf4::ncvar_def("FSH", "Wm-2", list(lon,lat,time), mv,
                           longname="sensible heat flux", prec="double")
  EFLX_LH_TOT  <- ncdf4::ncvar_def("EFLX_LH_TOT", "Wm-2", list(lon,lat,time), mv,
                                   longname="latent heat flux", prec="double")
  GPP <- ncdf4::ncvar_def("GPP", "umolm-2s-1", list(lon,lat,time), mv,
                          longname="gross primary productivity", prec="double")
  Rnet  <- ncdf4::ncvar_def("Rnet", "W/m^2", list(lon,lat,time), mv,
                            longname="net radiation", prec="double")
  
  #Create the output file
  ncnew <- ncdf4::nc_create(fileOutNcdf, list(LATIXY,LONGXY,FLDS,FSDS,PRECTmms,RH,PSRF,TBOT,WIND,ZBOT,FSH,EFLX_LH_TOT,NEE,GPP,Rnet))
  
  
  # Write some values to this variable on disk.
  ncdf4::ncvar_put(ncnew, LATIXY, latSite)
  ncdf4::ncvar_put(ncnew, LONGXY, lonSiteE)
  ncdf4::ncvar_put(ncnew, FLDS, Data.mon$FLDS) 
  ncdf4::ncvar_put(ncnew, FSDS, Data.mon$Rg)
  ncdf4::ncvar_put(ncnew, RH,   Data.mon$rH)
  ncdf4::ncvar_put(ncnew, PRECTmms, Data.mon$PRECTmms)
  ncdf4::ncvar_put(ncnew, PSRF, Data.mon$PA)
  ncdf4::ncvar_put(ncnew, TBOT, Data.mon$Tair)
  ncdf4::ncvar_put(ncnew, WIND, Data.mon$U)
  ncdf4::ncvar_put(ncnew, ZBOT, Data.mon$ZBOT)
  ncdf4::ncvar_put(ncnew, NEE, Data.mon$NEE)
  ncdf4::ncvar_put(ncnew, FSH, Data.mon$H)
  ncdf4::ncvar_put(ncnew, EFLX_LH_TOT, Data.mon$LE)
  ncdf4::ncvar_put(ncnew, GPP, Data.mon$GPP)
  ncdf4::ncvar_put(ncnew, Rnet, Data.mon$radNet)
  #add attributes
  # ncdf4::ncatt_put(ncnew, time,"calendar", "gregorian" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, FLDS,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, FSDS,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, RH  ,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, PRECTmms,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, PSRF,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, TBOT,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, WIND,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, ZBOT,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, NEE,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, FSH,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, EFLX_LH_TOT,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, GPP,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, Rnet,"mode","time-dependent" ,prec=NA,verbose=FALSE,definemode=FALSE )
  #ncdf4::ncatt_put(ncnew, 0, "veg_community_type", veg_community_list[veg_community],prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, 0, "created_on",date(),prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, 0, "created_by","Derek Pierson",prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, 0, "flux_tower_site_name",tower_name ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, 0, "created_from",flux_file_name ,prec=NA,verbose=FALSE,definemode=FALSE )
  ncdf4::ncatt_put(ncnew, 0, "created_with", "prepare_forcings_for_clm.R",prec=NA,verbose=FALSE,definemode=FALSE )
  
  
  #Close Netcdf file connection
  ncdf4::nc_close(ncnew)
  #Add step
  #startStep <- endStep + 1
  #Remove not needed variables
  remove(time, timeStep, fileOutNcdf, ncnew, Data.mon,
         FLDS,FSDS,RH,PRECTmms,PSRF,TBOT,WIND,ZBOT)
} #End of monthloop


############################################################################################
#save dataframe as .csv
write.csv(dataClm, "C:/Users/Derek/Google Drive/RCrk/RC Flux Tower data/RC-Rls/gap-filled_csv_rds/rcls_dataClm.csv")

#save as RDS
saveRDS(dataClm, file = "C:/Users/Derek/Google Drive/RCrk/RC Flux Tower data/RC-Rls/gap-filled_csv_rds/rcls_dataClm.rds")


# Open and plot netCDF file to check time conversion
#############################################################
nc <- nc_open("C:/Users/Derek/Google Drive/RCrk/RC Flux Tower data/RC-Rls/gap-filled_monthly_netCDF/2016-04.nc")
nc_time <- ncvar_get(nc, "time")
nc_fsds <- ncvar_get(nc, "FSDS")
nc_df <- data.frame(time=nc_time, FSDS = nc_fsds)
nc_test_plt <- ggplot(nc_df[1:144,], aes(x=time, y=FSDS)) + geom_point() + geom_line() + ggtitle("netCDF data: April 2016")

#plot should match!
grid.arrange(time_plt_raw, time_plt_up7, nc_test_plt, ncol=1) 



