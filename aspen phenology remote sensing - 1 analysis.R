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
library(raster)
library(ncdf4)

try(dir.create('outputs'))

# this is a miserable hack due to some ncdf reading issues when i upgraded to R 4.2 and the packages changed.
# should work just as well with rast(fn) but apparently not...
read_ncdf <- function(fn, names_drop)
{
  nc_this = nc_open(fn)
  varnames_this = attributes(nc_this$var)$names
  print(varnames_this)
  varnames_this = setdiff(varnames_this, names_drop)
  r_all = lapply(varnames_this, function(v) { raster(fn,varname=v) })
  r_all = stack(r_all)
  names(r_all) = varnames_this
  
  fn_final = sprintf("~/Downloads/%s_%.6f.grd",basename(fn),runif(1))
  writeRaster(r_all,fn_final, format="raster", overwrite=TRUE,options=c("COMPRESS=LZW"))
  r_all = rast(fn_final)
  
  return(r_all)
}

# read in phenology and convert EVImax to fractions
names_drop_pheno = c("transverse_mercator","NumCycles","Peak","EVIamp","EVIarea","OGI_2","50PCGI_2","OGMx_2","Peak_2","OGD_2","50PCGD_2","OGMn_2","EVImax_2","EVIamp_2","EVIarea_2","numObs","gupQA_2","gdownQA_2")

r_phenology_2016 = read_ncdf('../../13SCD/MSLSP_13SCD_2016.nc', names_drop=names_drop_pheno)
r_phenology_2016[["EVImax"]] = r_phenology_2016[["EVImax"]] / 10000
r_phenology_2017 = read_ncdf('../../13SCD/MSLSP_13SCD_2017.nc', names_drop=names_drop_pheno)
r_phenology_2017[["EVImax"]] = r_phenology_2017[["EVImax"]] / 10000
r_phenology_2018 = read_ncdf('../../13SCD/MSLSP_13SCD_2018.nc', names_drop=names_drop_pheno)
r_phenology_2018[["EVImax"]] = r_phenology_2018[["EVImax"]] / 10000
r_phenology_2019 = read_ncdf('../../13SCD/MSLSP_13SCD_2019.nc', names_drop=names_drop_pheno)
r_phenology_2019[["EVImax"]] = r_phenology_2019[["EVImax"]] / 10000

# get variable names
vars_all = read.csv('../../13SCD/MSLSP_Layers_V0.csv') %>% 
  mutate(short_name = make.names(short_name)) %>%
  dplyr::select(long_name, short_name) %>%
  rbind(data.frame(long_name="Growing season length (50% greenup - 50% greendown)",short_name="GSL.50"))

vars_phenology = c("OGI", "OGMn", "GSL.50", "EVImax","gupQA","gdownQA")#c("OGI", "X50PCGI", "OGMx", "OGD", "X50PCGD", "OGMn", "EVImax", "GSL.15", "GSL.50", "GSL.90")


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




# switch aspen data to phenology projection 
r_neon_projected = project(r_neon, r_phenology_2016)
# trim neon data to aspen area
r_neon_trimmed = crop(r_neon_projected, trim(r_neon_projected["is_diploid"]))


# trim phenology data to relevant area

r_phenology_2016_trimmed = crop(r_phenology_2016, r_neon_trimmed)
r_phenology_2017_trimmed = crop(r_phenology_2017, r_neon_trimmed)
r_phenology_2018_trimmed = crop(r_phenology_2018, r_neon_trimmed)
r_phenology_2019_trimmed = crop(r_phenology_2019, r_neon_trimmed)


# clip phenology to aspen
r_phenology_2016_trimmed_masked = mask(r_phenology_2016_trimmed, r_neon_trimmed["is_diploid"])
r_phenology_2017_trimmed_masked = mask(r_phenology_2017_trimmed, r_neon_trimmed["is_diploid"])
r_phenology_2018_trimmed_masked = mask(r_phenology_2018_trimmed, r_neon_trimmed["is_diploid"])
r_phenology_2019_trimmed_masked = mask(r_phenology_2019_trimmed, r_neon_trimmed["is_diploid"])


# deal with phenology rasters
process_pheno <- function(r_pheno)
{
  names(r_pheno) = make.names(names(r_pheno))
  r_pheno[["GSL.50"]] = r_pheno[["X50PCGD"]] - r_pheno[["X50PCGI"]]
  r_pheno = r_pheno[[vars_phenology]]
  return(r_pheno)
}

r_pheno_all_trimmed = lapply(list(r_phenology_2016_trimmed, r_phenology_2017_trimmed, r_phenology_2018_trimmed, r_phenology_2019_trimmed),
                     process_pheno)
names(r_pheno_all_trimmed) = make.names(2016:2019)

r_pheno_all_trimmed_masked = lapply(list(r_phenology_2016_trimmed_masked, r_phenology_2017_trimmed_masked, r_phenology_2018_trimmed_masked, r_phenology_2019_trimmed_masked),
                             process_pheno)
names(r_pheno_all_trimmed_masked) = make.names(2016:2019)




# load in temp metrics
fn_tmax = dir(path='../../gridmet',pattern="*.nc",full.names = TRUE)
r_tmax_list = lapply(fn_tmax, stack) # avoid reading ncdf with rast
# make things smaller to go faster
r_tmax_list = lapply(r_tmax_list, function(x) {
  crop(x, y=c(-109,-105,37,41))
  })

# weird hack due to base raster package use
r_tmax_q99 = lapply(r_tmax_list, function(x) {calc(stack(x), fun=function(y) {quantile(y, 0.99, na.rm=TRUE)})})
r_tmax_q99 = rast(lapply(r_tmax_q99,rast))
r_tmax_q99_projected = project(r_tmax_q99, r_phenology_2016)
r_tmax_q99_trimmed = crop(r_tmax_q99_projected, trim(r_neon_trimmed["is_diploid"]))
r_tmax_q99_trimmed = r_tmax_q99_trimmed / 10 # convert to deg C
names(r_tmax_q99_trimmed) = make.names(2012:2019)
# also make masked version
r_tmax_q99_trimmed_masked = mask(r_tmax_q99_trimmed, r_neon_trimmed["is_diploid"])

# load in soil moisture metrics
r_sm_q01 = rast('metrics/metric_annual_0.1m_rasters_all_q01.tif')
r_sm_q01_projected = project(r_sm_q01, r_phenology_2016)
r_sm_q01_trimmed = crop(r_sm_q01_projected, trim(r_neon_trimmed["is_diploid"]))
r_sm_q01_trimmed_masked = mask(r_sm_q01_trimmed, r_neon_trimmed["is_diploid"])

r_sm_runs_med_dur = rast('metrics/metric_annual_0.1m_runs_med.low.spell.duration.tif')
r_sm_runs_med_dur_projected = project(r_sm_runs_med_dur, r_phenology_2016)
r_sm_runs_med_dur_projected[is.na(r_sm_runs_med_dur_projected)] = 0
r_sm_runs_med_dur_trimmed = crop(r_sm_runs_med_dur_projected, r_neon_trimmed["is_diploid"])
r_sm_runs_med_dur_trimmed_masked = mask(r_sm_runs_med_dur_trimmed, r_neon_trimmed["is_diploid"])



# load in snowmelt metrics
r_snowmelt = rast(c('../../Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2012_V2.tif',
                    '../../Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2013_V2.tif',
                    '../../Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2014_V2.tif',
                    '../../Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2015_V2.tif',
                    '../../Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2016_V2.tif',
                    '../../Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2017_V2.tif',
                    '../../Snowmelt_timing_maps_V2_1712/data/Snowmelt_Timing_North_America_2018_V2.tif'))

r_snowmelt_2019 = rast('../../Snowmelt_timing_unpublished/Snowmelt_Timing_h09v05_2019.tif')
r_snowmelt_2019 = terra::extend(r_snowmelt_2019, ext(r_snowmelt))
r_snowmelt = c(r_snowmelt, r_snowmelt_2019)

names(r_snowmelt) = make.names(2012:2019)
r_snowmelt_projected = project(r_snowmelt, r_phenology_2016)
r_snowmelt_trimmed = crop(r_snowmelt_projected, trim(r_neon_trimmed["is_diploid"]))
r_snowmelt_trimmed_masked = mask(r_snowmelt_trimmed, r_neon_trimmed["is_diploid"])


# make output rasters for plotting
r_snowmelt_trimmed_for_output = r_snowmelt_trimmed
names(r_snowmelt_trimmed_for_output) = paste("snowmelt",make.names(2012:2019),sep=".")

r_tmax_q99_trimmed_for_output = r_tmax_q99_trimmed
names(r_tmax_q99_trimmed_for_output) = paste("tmax_q99",make.names(2012:2019),sep=".")

r_sm_q01_trimmed_for_output = r_sm_q01_trimmed
names(r_sm_q01_trimmed_for_output) = paste("sm_q01",make.names(2012:2019),sep=".")

r_sm_runs_med_dur_trimmed_for_output = r_sm_runs_med_dur_trimmed
names(r_sm_runs_med_dur_trimmed_for_output) = paste("sm_runs_med_dur",make.names(2012:2019),sep=".")

# write out predictor rasters
r_env_all = c(r_neon_trimmed, 
              r_snowmelt_trimmed_for_output,
              r_tmax_q99_trimmed_for_output,
              r_sm_q01_trimmed_for_output,
              r_sm_runs_med_dur_trimmed_for_output)

for (i in 1:nlyr(r_env_all))
{
  writeRaster(r_env_all[[i]], file=sprintf('outputs/r_pred_%s.tif',names(r_env_all)[i]),overwrite=TRUE)
}

# write out phenology rasters, trimmed
for (i in 1:length(r_pheno_all_trimmed))
{
  for (j in 1:nlyr(r_pheno_all_trimmed[[i]]))
  {
    writeRaster(r_pheno_all_trimmed[[i]][[j]], 
                file=sprintf('outputs/r_pheno_all_trimmed_%s_%s.tif',names(r_pheno_all_trimmed[[i]])[[j]],names(r_pheno_all_trimmed)[i]),overwrite=TRUE)
  }
}

# write out phenology rasters, trimmed and masked
for (i in 1:length(r_pheno_all_trimmed_masked))
{
  for (j in 1:nlyr(r_pheno_all_trimmed_masked[[i]]))
  {
    writeRaster(r_pheno_all_trimmed_masked[[i]][[j]], 
                file=sprintf('outputs/r_pheno_all_trimmed_masked_%s_%s.tif',names(r_pheno_all_trimmed_masked[[i]])[[j]],names(r_pheno_all_trimmed_masked)[i]),overwrite=TRUE)
  }
}





get_raster_for_years <- function(r, year, num_lags, varname) # assumes names are X2015, X2016, ...
{
  years_this = paste(varname,make.names((year-num_lags+1):year),sep=".")

  r_out = r[[years_this]]
  print(names(r_out)) #to verify output
  names(r_out) = sprintf("%s.t.minus.%d", varname, (num_lags-1):0)
  return(r_out)
}


process_data <- function(r, year, 
                         num_lags_moisture = 4, 
                         num_lags_snowmelt = 4, 
                         num_lags_temp = 4, 
                         aspen_cover_min = 0.25, 
                         qa_max = 2)
{
  r_sm_q01 = get_raster_for_years(r_sm_q01_trimmed_for_output, varname="sm_q01",year=year, num_lags=num_lags_moisture)
  r_sm_runs_med_dur = get_raster_for_years(r_sm_runs_med_dur_trimmed_for_output, varname="sm_runs_med_dur",year=year, num_lags=num_lags_moisture)
  r_snowmelt_date = get_raster_for_years(r_snowmelt_trimmed_for_output, varname="snowmelt",year=year, num_lags=num_lags_snowmelt)
  r_maxtemp = get_raster_for_years(r_tmax_q99_trimmed_for_output, varname="tmax_q99",year=year, num_lags=num_lags_temp)
  
  r_joined = c(r,r_neon_trimmed, r_maxtemp, r_snowmelt_date, r_sm_q01, r_sm_runs_med_dur)
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
    mutate(GSL.50 = X50PCGD - X50PCGI) %>%
    dplyr::select(year, x, y, aspen_cover, 
                  all_of(vars_phenology), 
                  all_of(names(r_neon)),
                  all_of(names(r_maxtemp)), 
                  all_of(names(r_snowmelt_date)), 
                  all_of(names(r_sm_q01)),
                  all_of(names(r_sm_runs_med_dur))) %>%
    rename(cytotype_fraction_diploid=is_diploid)
  
  return(data)
}
 
# put all datasets together after QC, masking to aspen cover
df_all_2016 = process_data(r_phenology_2016_trimmed_masked, 2016)
df_all_2017 = process_data(r_phenology_2017_trimmed_masked, 2017)
df_all_2018 = process_data(r_phenology_2018_trimmed_masked, 2018)
df_all_2019 = process_data(r_phenology_2019_trimmed_masked, 2019)

df_all = rbind(df_all_2016, df_all_2017, df_all_2018, df_all_2019)

# make counts
df_all %>% group_by(year) %>% summarize(n())

write.csv(df_all, 'outputs/df_all_aspen_cover_0.25.csv',row.names = FALSE)




