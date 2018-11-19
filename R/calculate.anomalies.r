################################################################
# calculate.anomalies.r 

################################################################

#---------------------------------------------------------------
#' File Time Series Function
#'
#' Reads the netcdf file time dimensions as a PCICt time
#' object and returns a PCICt time series
#' @param nc ncdf4 object of opened netcdf file.
#' @keywords time
#' @importFrom PCICt as.PCICt
#' @import ncdf4
#' @export

pcict_time <- function(nc) {
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)  
  dates <- origin.pcict + time.values*86400
}
#---------------------------------------------------------------

#---------------------------------------------------------------
#' Find Interval Dates Function
#'
#' Find the start and end dates for a time series given an 
#' interval with a start and end year.
#' @param dates PCICt time series.
#' @param yst Start year
#' @param yen End year
#' @keywords interval
#' @importFrom utils  head tail
#' @export

find_interval_dates <- function(dates,yst,yen) {
  years <- format(dates,'%Y')
  st.ix <- head(grep(yst,years),1)
  en.ix <- tail(grep(yen,years),1)
  rv <- list(start=st.ix,end=en.ix,count=en.ix-st.ix+1)
  return(rv)
}
#---------------------------------------------------------------


#---------------------------------------------------------------
#' Extend Area Edge for Interpolation Function
#'
#' Repeats values at the edge of the file area to add a buffer
#' for the interpolation step. Only necessary for edges with
#' missing data.
#' @param data Input anomaly data
#' @param n Defaults to 3.
#' @keywords buffer
#' @export

edge_buffer <- function(data,n=3) {
  
  ncol <- dim(data)[2]
  for (l in 1:n) { ##Repeat n times to add buffer
    for (j in 1:(ncol-1)) {
      ix <- is.na(data[,j,1])
      data[ix,j,] <- data[ix,(j+1),]    
    }
    for (j in ncol:2) {
      ix <- is.na(data[,j,1])
      data[ix,j,] <- data[ix,(j-1),]    
    }

    nrow <- dim(data)[1]
    for (k in 1:(nrow-1)) {
      kx <- is.na(data[k,,1])
      data[k,kx,] <- data[(k+1),kx,]    
    }
    for (k in nrow:2) {
      kx <- is.na(data[k,,1])
      data[k,kx,] <- data[(k-1),kx,]    
    }   
  }
  return(data)
}
#---------------------------------------------------------------

#---------------------------------------------------------------
#' Create Anomaly File Function
#'
#' Calculates monthly anomalies of daily values relative to
#' a climatological period. 
#' @param  var.name Variable name for the netcdf file.
#' @param  series.file File name for input netcdf file.
#' @param  read.dir Location for the input file.
#' @param  anomaly.file File name for anomaly netcdf file.
#' @param  write.dir Location for the anomaly file.
#' @param  monthly.clim Monthly climatologies for the baseline.
#' @param  buffer.edge Fill missing data at edges DEFAULT=NULL
#' @keywords anomalies
#' @import ncdf4
#' @export

create_anomalies <- function(var.name,series.file,read.dir,anomaly.file,write.dir,monthly.clim,buffer.edge) {

  ## Create a copy of the daily file to receive the anomalies
  file.copy(from=paste0(read.dir,series.file),
            to=paste0(write.dir,anomaly.file),overwrite=T)
  Sys.sleep(5)
  if (!file.exists(paste0(read.dir,anomaly.file)))
    warning('Anomaly file was not created')
   
  nc <- nc_open(paste0(write.dir,series.file))
  anc <- nc_open(paste0(write.dir,anomaly.file),write=TRUE)
 
  lat <- ncvar_get(nc,'lat')  
  n.lat <- length(lat)
  l.seq <- seq(1,n.lat,by=1) ###

  lon <- ncvar_get(nc,'lon')
  n.lon <- length(lon)

  var.dates <- pcict_time(nc)
  n.time <- length(var.dates)
  monthly.fac <- as.factor(format(var.dates,'%m'))

  ##Iterate over latitude bands to lower memory usage
  for (ltx in l.seq) {
    var.data <- ncvar_get(nc,var.name,start=c(1,ltx,1),count=c(-1,1,-1))    

    ##Load full time series and take anomalies from this
    var.anoms <- var.data*0
    for(mn in 1:12) {
      var.ix <- which(monthly.fac==sprintf('%02d',mn))
      mlen <- length(var.ix)      
      monthly.mean <- t(sapply(monthly.clim[,ltx,mn],rep,mlen))
      if (var.name=='pr') {
        var.anoms[,var.ix] <- var.data[,var.ix]/monthly.mean
      }
      if (grepl('tas',var.name)) {
          var.anoms[,var.ix] <- var.data[,var.ix] - monthly.mean
      }
    }
    if (var.name=='pr') {
      var.anoms[is.na(var.anoms)] <- 0
    }
    if (!is.null(buffer.edge)) {
      var.anoms <- edge_buffer(var.anoms,n=buffer.edge)
    }
    ncvar_put(anc,varid=var.name,vals=var.anoms,start=c(1,ltx,1),count=c(n.lon,1,n.time)) 
  }
  nc_close(nc)
  nc_close(anc)
  gc()
}

#---------------------------------------------------------------
#' Make Anomalies Function
#'
#' Calculates monthly anomalies of daily values relative to
#' a climatological period. 
#' @param  var.name Variable name for the netcdf file.
#' @param  series.file File name for input netcdf file.
#' @param  read.dir Location for the input file.
#' @param  anomaly.file File name for anomaly netcdf file.
#' @param  write.dir Location for the anomaly file.
#' @param  yst Climatology starting year
#' @param  yen Climatology ending year
#' @param  buffer.edge Fill missing data at edges DEFAULT=NULL
#' @keywords anomalies
#' @import ncdf4
#' @export

make_anomalies <- function(var.name,
                           series.file,read.dir,
                           anomaly.file,write.dir,
                           yst,yen,buffer.edge=NULL) {

  nc <- nc_open(series.file)
  var.dates <- pcict_time(nc)
  interval.bounds <- find_interval_dates(var.dates,yst,yen)
  var.data <- ncvar_get(nc,var.name)
  clim.data <- ncvar_get(nc,var.name,start=c(1,1,interval.bounds$start),
                                    count=c(-1,-1,interval.bounds$count))
  clim.dates <- var.dates[interval.bounds$start:interval.bounds$end]

  ##------------------------------------  
  monthly.fac <- as.factor(format(clim.dates,'%m'))

  ##If using precipitation data apply monthly totals otherwise use monthly means
  if (var.name=='pr') {
    monthly.ts.fac <- as.factor(format(clim.dates,'%Y-%m'))
    mon.facs <- as.factor(format(as.Date(paste(levels(monthly.ts.fac),'-01',sep='')),'%m'))
    clim.data[clim.data <=0] <- NA    
    clim.total <- apply(clim.data,c(1,2),function(x,fac){tapply(x,fac,sum,na.rm=T)},monthly.ts.fac)
    monthly.clim <-  aperm(apply(clim.total,c(2,3),function(x,fac){tapply(x,fac,mean,na.rm=T)},mon.facs),c(2,3,1))
  } else {
    monthly.clim <- aperm(apply(clim.data,c(1,2),function(x,fac){tapply(x,fac,mean,na.rm=T)},monthly.fac),c(2,3,1))
  }
  nc_close(nc)

  ##------------------------------------  
  create_anomalies(var.name,series.file,read.dir,anomaly.file,write.dir,monthly.clim,buffer.edge)

}

#---------------------------------------------------------------

