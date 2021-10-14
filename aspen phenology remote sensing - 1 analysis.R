library(terra)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggpubr)
library(ranger)
library(pdp)
library(viridis)
library(tidyr)
library(caret)

try(dir.create('rasters_for_plotting'))

r_phenology_2016 = rast('../13SCD/MSLSP_13SCD_2016.nc')
r_phenology_2017 = rast('../13SCD/MSLSP_13SCD_2017.nc')
r_phenology_2018 = rast('../13SCD/MSLSP_13SCD_2018.nc')
r_phenology_2019 = rast('../13SCD/MSLSP_13SCD_2019.nc')

# get variable names
vars_all = read.csv('../13SCD/MSLSP_Layers_V0.csv') %>% 
  mutate(short_name = make.names(short_name)) %>%
  select(long_name, short_name) %>%
  rbind(data.frame(long_name="Growing season length (50% greenup - 50% greendown)",short_name="GSL.50"))

vars_phenology = c("OGI", "OGMn", "GSL.50", "EVImax")#c("OGI", "X50PCGI", "OGMx", "OGD", "X50PCGD", "OGMn", "EVImax", "GSL.15", "GSL.50", "GSL.90")


# load topography rasters
r_elev = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/layers/dtm_mosaic_min_phase_me.tif')
r_cos_aspect = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/layers/r_cos_aspect.tif')
r_slope = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/layers/r_slope.tif')
r_height = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/layers/tch_mosaic_min_phase_me.tif')

# load cytotype
r_cytotype_masked = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/r_cytotype_masked.tif') 
# 0 = triploid, 1=diploid

# make aspen cover
r_aspen_cover = !is.na(r_cytotype_masked)


r_neon = c(r_elev, r_cos_aspect, r_slope, r_height, r_cytotype_masked, r_aspen_cover)
names(r_neon) = c("elevation","cos_aspect","slope","height_canopy","is_diploid", "aspen_cover")




# switch aspen data to 30 m resolution and clip to aspen
r_neon_projected = project(r_neon, r_phenology_2016)
# trim neon data to relevant area
r_neon_trimmed = crop(r_neon_projected, trim(r_neon_projected["is_diploid"]))

# check cover predictions
#r_neon_trimmed["aspen_cover"] %>% hist(ylim=c(0,1e4))
# fix strange outlier issues after projection
#r_neon_trimmed["cos_aspect"] = clamp(r_neon_trimmed["cos_aspect"],lower=-1,upper=1)
#r_neon_trimmed["slope"] = clamp(r_neon_trimmed["slope"],lower=0,upper=90)

# clip phenology to aspen
r_phenology_2016_masked = mask(r_phenology_2016, r_neon_projected["is_diploid"])
r_phenology_2017_masked = mask(r_phenology_2017, r_neon_projected["is_diploid"])
r_phenology_2018_masked = mask(r_phenology_2018, r_neon_projected["is_diploid"])
r_phenology_2019_masked = mask(r_phenology_2019, r_neon_projected["is_diploid"])


# trim phenology data to relevant area
r_phenology_2016_trimmed = crop(r_phenology_2016_masked, trim(r_neon_projected["is_diploid"]))
r_phenology_2017_trimmed = crop(r_phenology_2017_masked, trim(r_neon_projected["is_diploid"]))
r_phenology_2018_trimmed = crop(r_phenology_2018_masked, trim(r_neon_projected["is_diploid"]))
r_phenology_2019_trimmed = crop(r_phenology_2019_masked, trim(r_neon_projected["is_diploid"]))


# make time series plots
extract_ts_raster <- function(yvar)
{
  if (yvar=="GSL.50")
  {
    r_phenology_2016_trimmed = c(r_phenology_2016_trimmed, r_phenology_2016_trimmed[["50PCGD"]] - r_phenology_2016_trimmed[["50PCGI"]])
    r_phenology_2017_trimmed = c(r_phenology_2017_trimmed, GSL.50 = r_phenology_2017_trimmed[["50PCGD"]] - r_phenology_2017_trimmed[["50PCGI"]])
    r_phenology_2018_trimmed = c(r_phenology_2018_trimmed, GSL.50 = r_phenology_2018_trimmed[["50PCGD"]] - r_phenology_2018_trimmed[["50PCGI"]])
    r_phenology_2019_trimmed = c(r_phenology_2019_trimmed, GSL.50 = r_phenology_2019_trimmed[["50PCGD"]] - r_phenology_2019_trimmed[["50PCGI"]])
    names(r_phenology_2016_trimmed)[nlyr(r_phenology_2016_trimmed)] = "GSL.50"
    names(r_phenology_2017_trimmed)[nlyr(r_phenology_2017_trimmed)] = "GSL.50"
    names(r_phenology_2018_trimmed)[nlyr(r_phenology_2018_trimmed)] = "GSL.50"
    names(r_phenology_2019_trimmed)[nlyr(r_phenology_2019_trimmed)] = "GSL.50"
  }  
  r_yvar = c(r_phenology_2016_trimmed[[yvar]], r_phenology_2017_trimmed[[yvar]], r_phenology_2018_trimmed[[yvar]], r_phenology_2019_trimmed[[yvar]]) 
  
  names(r_yvar) = make.names(2016:2019)
  
  
  writeRaster(r_yvar, file=sprintf('rasters_for_plotting/r_pheno_%s.tif',yvar))
  

  return(r_yvar)
}

# make plots of each phenology predictor by time
r_pheno_all_2016_2019 = lapply(vars_phenology, extract_ts_raster)



# load in temp metrics
fn_tmax = dir(path='../gridmet',pattern="*.nc",full.names = TRUE)
r_tmax_list = lapply(fn_tmax, rast)
r_tmax_q99 = lapply(r_tmax_list, quantile, 0.99, na.rm=TRUE)
r_tmax_q99 = rast(r_tmax_q99)
r_tmax_q99_projected = project(r_tmax_q99, r_phenology_2016)
r_tmax_q99_masked = mask(r_tmax_q99_projected, r_neon_projected["is_diploid"])
r_tmax_q99_trimmed = crop(r_tmax_q99_masked, trim(r_neon_projected["is_diploid"]))
names(r_tmax_q99_trimmed) = make.names(2012:2019)

# load in soil moisture metrics
r_sm_q01 = rast('metrics/metric_annual_0.1m_rasters_all_q01.tif')
r_sm_q01_projected = project(r_sm_q01, r_phenology_2016)
r_sm_q01_masked = mask(r_sm_q01_projected, r_neon_projected["is_diploid"])
r_sm_q01_trimmed = crop(r_sm_q01_masked, trim(r_neon_projected["is_diploid"]))

r_sm_runs_med_dur = rast('metrics/metric_annual_0.1m_runs_med.low.spell.duration.tif')
r_sm_runs_med_dur_projected = project(r_sm_runs_med_dur, r_phenology_2016)
r_sm_runs_med_dur_projected[is.infinite(r_sm_runs_med_dur_projected)] = 0
r_sm_runs_med_dur_masked = mask(r_sm_runs_med_dur_projected, r_neon_projected["is_diploid"])
r_sm_runs_med_dur_trimmed = crop(r_sm_runs_med_dur_masked, trim(r_neon_projected["is_diploid"]))

r_sm_runs_max_dur = rast('metrics/metric_annual_0.1m_runs_max.low.duration.tif')
r_sm_runs_max_dur_projected = project(r_sm_runs_max_dur, r_phenology_2016)
r_sm_runs_max_dur_projected[is.infinite(r_sm_runs_max_dur_projected)] = 0
r_sm_runs_max_dur_masked = mask(r_sm_runs_max_dur_projected, r_neon_projected["is_diploid"])
r_sm_runs_max_dur_trimmed = crop(r_sm_runs_max_dur_masked, trim(r_neon_projected["is_diploid"]))



# load in snowmelt metrics
r_snowmelt = rast(c('../Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2012_V2.tif',
                    '../Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2013_V2.tif',
                    '../Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2014_V2.tif',
                    '../Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2015_V2.tif',
                    '../Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2016_V2.tif',
                    '../Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2017_V2.tif',
                    '../Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2018_V2.tif'))

r_snowmelt_2019 = rast('../Snowmelt_timing_unpublished/Snowmelt_Timing_h09v05_2019.tif')
r_snowmelt_2019 = terra::extend(r_snowmelt_2019, ext(r_snowmelt))
r_snowmelt = c(r_snowmelt, r_snowmelt_2019)

names(r_snowmelt) = make.names(2012:2019)
r_snowmelt_projected = project(r_snowmelt, r_phenology_2016)
r_snowmelt_masked = mask(r_snowmelt_projected, r_neon_projected["is_diploid"])
r_snowmelt_trimmed = crop(r_snowmelt_masked, trim(r_neon_projected["is_diploid"]))



# make output rasters for plotting
r_snowmelt_trimmed_for_output = r_snowmelt_trimmed
names(r_snowmelt_trimmed_for_output) = paste("snowmelt",make.names(2012:2019),sep=".")

r_tmax_q99_trimmed_for_output = r_tmax_q99_trimmed
names(r_tmax_q99_trimmed_for_output) = paste("tmax_q99",make.names(2012:2019),sep=".")

r_sm_q01_trimmed_for_output = r_sm_q01_trimmed
names(r_sm_q01_trimmed_for_output) = paste("sm_q01",make.names(2012:2019),sep=".")


r_env_all = c(r_neon_trimmed, 
              r_snowmelt_trimmed_for_output,
              r_tmax_q99_trimmed_for_output,
              r_sm_q01_trimmed_for_output)
  
for (i in 1:nlyr(r_env_all))
{
  writeRaster(r_env_all[[i]], file=sprintf('rasters_for_plotting/r_pred_%s.tif',names(r_env_all)[i]))
}





get_raster_for_years <- function(r, year, num_lags, varname) # assumes names are X2015, X2016, ...
{
  years_this = make.names((year-num_lags+1):year)
  r_out = r[[years_this]]
  print(names(r_out)) #to verify output
  names(r_out) = sprintf("%s.t.minus.%d", varname, (num_lags-1):0)
  return(r_out)
}


process_data <- function(r, year, num_lags_moisture=3, num_lags_snowmelt=3, num_lags_temp = 3, aspen_cover_min = 0.5, qa_max = 2)
{
  r_soilmoisture_q01 = get_raster_for_years(r_sm_q01_trimmed, varname="soilmoisture0.1.q01",year=year, num_lags=num_lags_moisture)
  r_soilmoisture_run.med = get_raster_for_years(r_sm_runs_med_dur_trimmed, varname="soilmoisture0.1.run.med",year=year, num_lags=num_lags_moisture)
  r_snowmelt_date = get_raster_for_years(r_snowmelt_trimmed, varname="snowmelt.date",year=year, num_lags=num_lags_snowmelt)
  r_maxtemp = get_raster_for_years(r_tmax_q99_trimmed, varname="tmax.q99",year=year, num_lags=num_lags_temp)
  
  r_joined = c(r,r_neon_trimmed, r_maxtemp, r_snowmelt_date, r_soilmoisture_q01, r_soilmoisture_run.med)
  data = as.data.frame(r_joined[])
  
  xy = xyFromCell(r_joined, 1:ncell(r_joined))
  
  data = cbind(data, xy)
  
  names(data) = make.names(names(data))
  
  data = data %>% 
    # only pixels that are likely to be unshadowed during mid-morning sentinel-2 overpass - check if realistic
    filter(!is.na(elevation) & aspen_cover > aspen_cover_min) %>%  #slope < 60 & height_canopy > 2 & cos_aspect < 0.9 & 
    # see # page 12 of MSLSP_User_Guide_V1.pdf - QA 1, high quality; QA 2, moderate quality, 3 poor quality with successful fill
    filter(gupQA <= qa_max & gdownQA <= qa_max) %>% 
    # add a year column
    mutate(year=factor(year,ordered=FALSE)) %>%
    # make a growing season column
    mutate(GSL.15 = OGMn - OGI) %>%
    mutate(GSL.50 = X50PCGD - X50PCGI) %>%
    mutate(GSL.90 = OGD - OGMx) %>%
    select(year, x, y, all_of(vars_phenology), all_of(names(r_neon)),all_of(names(r_maxtemp)), all_of(names(r_snowmelt_date)), all_of(names(r_soilmoisture_q01)), all_of(names(r_soilmoisture_run.med))) %>%
    mutate(cytotype=factor(is_diploid>0.5,levels=c(FALSE,TRUE),labels=c("triploid","diploid"),ordered=FALSE)) %>%
    select(-is_diploid)
  
  return(data)
}

# put all datasets together after QC
df_all_2016 = process_data(r_phenology_2016_trimmed, 2016)
df_all_2017 = process_data(r_phenology_2017_trimmed, 2017)
df_all_2018 = process_data(r_phenology_2018_trimmed, 2018)
df_all_2019 = process_data(r_phenology_2019_trimmed, 2019)

df_all = rbind(df_all_2016, df_all_2017, df_all_2018, df_all_2019)

write.csv(df_all, 'df_all.csv',row.names = FALSE)




