# note this assumes we removed the 0930 date from the raw soil moisture files per Hoang on 5/27/2022...

library(terra)
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)
library(hydrostats)
library(lubridate)

process_year <- function(year, band) # band 1 = 0.1 m, band 2 = 0.4 m, band 3 = 1 m, etc
{
  basedir = sprintf('../../SOIL MOISTURE/%d/',year)
  
  files = dir(basedir)
  print(files)

  # select a given depth
  rasters = lapply(file.path(basedir,files), function(x) { rast(x)[[band]] })
  
  # stack up rasters
  rasters = do.call("c", rasters)
  names(rasters) = files
  
  cat('.')
  
  return(rasters)
}  

years = 2012:2019
rasters_all <- lapply(years, process_year, band=1)
names(rasters_all) = make.names(years)


rasters_all_q01 = do.call("c",lapply(1:length(rasters_all), function(x) { quantile(rasters_all[[x]],0.01) }))
names(rasters_all_q01) = make.names(years)

# example time series extraction
#as.array(rasters_all[[1]])[10,7,] %>% plot(type='l')

# create time series on daily basis for hydrostats package for a given pixel i,j
make_ts_from_pixel <- function(raster_stack_for_year, i,j)
{
  vals = as.array(raster_stack_for_year)[i,j,]
  times = parse_date_time(gsub("\\.tif","",gsub("SM_","",names(rasters_all[[1]]))), "ymd")
  vals = vals[1:365]
  times = times[1:365]
  
  #print(paste(i, j, length(vals), length(times)))
  
  return(data.frame(Date=times,Q=vals))
}

# see https://cran.r-project.org/web/packages/hydrostats/hydrostats.pdf
make_drought_run_map <- function(rasters_year, threshold_sm)
{
  # set up empty rasters
  r_med.low.spell.duration = rasters_year[[1]]
  r_med.low.spell.duration[] = NA
  
  r_max.low.duration = rasters_year[[1]]
  r_max.low.duration[] = NA
  
  r_low.spell.freq = rasters_year[[1]]
  r_low.spell.freq[] = NA
  
  # iterate through pixels
  for (i in 1:nrow(rasters_year[[1]]))
  {
    for (j in 1:ncol(rasters_year[[1]]))
    {
      ts_this = make_ts_from_pixel(rasters_year,i=i,j=j)
      ls_this = low.spells(ts_this, threshold=threshold_sm, plot=FALSE)
      
      r_med.low.spell.duration[i,j] = ls_this$med.low.spell.duration
      r_max.low.duration[i,j] = ls_this$max.low.duration
      r_low.spell.freq[i,j] = ls_this$low.spell.freq
    }
    print(i/nrow(rasters_year[[1]]))
  }
  r_out = c(r_med.low.spell.duration, r_max.low.duration, r_low.spell.freq)
  names(r_out) = c("med.low.spell.duration", "max.low.duration", "low.spell.freq")
    
  r_out = mask(r_out, rasters_year[[1]])
  
  return(r_out)
}

# calculate the runs
runs_all_years = lapply(rasters_all, make_drought_run_map, threshold_sm = 0.2) # so this is 30% SM at 0.1 m depth

# break out by variable
runs_med.low.spell.duration = rast(lapply(runs_all_years, function(x) {x[["med.low.spell.duration"]]}))
names(runs_med.low.spell.duration) = make.names(years)

runs_max.low.duration = rast(lapply(runs_all_years, function(x) {x[["max.low.duration"]]}))
names(runs_max.low.duration) = make.names(years)

runs_low.spell.freq = rast(lapply(runs_all_years, function(x) {x[["low.spell.freq"]]}))
names(runs_low.spell.freq) = make.names(years)



# write out rasters
plot_raster <- function(var)
{
  if(!file.exists('metrics'))
  {
    dir.create('metrics')
  }
  vn = deparse(substitute(var))
  
  png(width=1000,height=800,file=sprintf("metrics/metric_annual_0.1m_%s.png",vn))
  par(oma=c(0,0,3,0))
  plot(var)
  title(vn,outer=TRUE)
  dev.off()
  
  writeRaster(var,file=sprintf("metrics/metric_annual_0.1m_%s.tif",vn),overwrite=TRUE)
}

plot_raster(rasters_all_q01)
plot_raster(runs_med.low.spell.duration)
plot_raster(runs_max.low.duration)
plot_raster(runs_low.spell.freq)