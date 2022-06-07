library(dplyr)
library(terra)
library(raster)
library(ncdf4)

simplify_geology <- function(unit)
{
  if (unit %in% c("Qal","Qg","Qu","Qlu","Qdf"))
  {
    #return("Quaternary alluvium or outwash gravel or deltas or undifferentiated")
    #return("Quaternary sorted")
    return("Quaternary deposit")
  }
  else if (unit %in% c("Ql","Qlf","Qd","Qdu"))
  {
    #return("Quaternary landslide or debris flow deposit")
    #return("Quaternary unsorted")
    return("Quaternary deposit")
  }
  else if (unit %in% c("Qm","Qmy"))
  {
    #return("Quaternary glacial deposit")
    #return("Quaternary unsorted")
    return("Quaternary deposit")
  }
  else if (unit %in% c("Qt","Qr"))
  {
    #return("Quaternary talus or rock glacier")
    return("Quaternary talus / rock glacier")
  }
  else if (unit %in% c("Tt","Tg","Tf","Tp","Tcd","Tbx","Td","qm"))
  {
    return("Igneous intrusive")
  }
  else if (unit %in% c("Kmu","Kml","Km"))
  {
    return("Shale or limestone")  
  }
  else if (unit %in% c("Kmf"))
  {
    #return("Fort Hays limestone") 
    return("Shale or limestone") # because we have no diploids on limestone - messes with mixed model
  }
  else if (unit %in% c("Kd","Je","Kmv","Kmvo"))
  {
    #return("Sandstone") 
    return("Sandstone or conglomerate or siltstone")
  }
  else if (unit %in% c("Jm","PPm","PPM")) #PPM is probably a typo
  {
    #return("Conglomerate or siltstone")
    return("Sandstone or conglomerate or siltstone")
  }
  else
  {
    return("")
  }
}

df_site_raw = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/aspen data site-level processed 30 Mar 2020.csv')

df_site = df_site_raw %>%
  dplyr::select(Site_Code, x=X.UTM, y=Y.UTM, Elevation, Cos.aspect, Slope, Summer.Insolation, DBH.mean, Canopy_openness, Soil.type, Geologic.Unit, Cytotype = Ploidy_level)

df_sex = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/aspen sex markers/aspen_sex_aug_11_2021.csv') %>%
  dplyr::select(Site_Code, geneticSexID)

df_site = df_site %>%
  left_join(df_sex, by="Site_Code")

# simplify the geology
df_site$Rock_Unit = factor(sapply(df_site$Geologic.Unit,simplify_geology))

df_site = df_site %>% 
  dplyr::select(-Geologic.Unit)






files_pheno_raw = dir(path='outputs',pattern='r_pheno',full.names = TRUE)
# keep the unmasked entries as the ground based plots don't necessarily overlap the 
# remotely sensed and clipped files
files_pheno = files_pheno_raw[!grepl("masked",files_pheno_raw)]

select_files_for_year <- function(files, year)
{
  files_this = files[grep(year,files)]
  
  #files_this = files_this[!grepl("QA",files_this)]
  
  r_this = rast(files_this)
  
  return(r_this)
}

r_pheno_2016 = select_files_for_year(files_pheno, "2016")
r_pheno_2017 = select_files_for_year(files_pheno, "2017")
r_pheno_2018 = select_files_for_year(files_pheno, "2018")
r_pheno_2019 = select_files_for_year(files_pheno, "2019")

extract_pheno_data <- function(r_pheno_this, year)
{
  e = terra::extract(x=r_pheno_this, y=df_site %>% dplyr::select(x,y))
  
  e$OGI[e$gupQA > 4] = NA
  e$OGMn[e$gdownQA > 4] = NA
  e$GSL.50[e$gdownQA > 4 | e$gupQA > 4] = NA
  # no QC needed for EVImax
  
  e = e %>% 
    dplyr::select(-gupQA, -gdownQA, -ID)
  
  names(e) = paste("pheno",year,names(e),sep=".")
  
  return(e)
}

df_pheno_2016 = extract_pheno_data(r_pheno_this=r_pheno_2016, year=2016)
df_pheno_2017 = extract_pheno_data(r_pheno_this=r_pheno_2017, year=2017)
df_pheno_2018 = extract_pheno_data(r_pheno_this=r_pheno_2018, year=2018)
df_pheno_2019 = extract_pheno_data(r_pheno_this=r_pheno_2019, year=2019)



df_site_combined = cbind(df_site, df_pheno_2016, df_pheno_2017, df_pheno_2018, df_pheno_2019)



# add in other rasters
r_sm_q01 = rast(dir('outputs',pattern='r_pred_sm_q01',full.names = TRUE))
r_sm_runs_med_dur = rast(dir('outputs',pattern='r_pred_sm_runs_med_dur',full.names = TRUE))
r_snowmelt = rast(dir('outputs',pattern='r_pred_snowmelt',full.names = TRUE))
r_tmax_q99 = rast(dir('outputs',pattern='r_pred_tmax_q99',full.names = TRUE))
r_cover = rast('outputs/r_pred_aspen_cover.tif')

r_all_non_pheno = c(r_sm_q01,
                    r_sm_runs_med_dur,
                    r_snowmelt,
                    r_tmax_q99,
                    r_cover)

df_other = terra::extract(x=r_all_non_pheno, y=df_site %>% dplyr::select(x,y)) %>%
  dplyr::select(-ID)


df_final = cbind(df_site_combined, df_other)

# write out file
write.csv(df_final, file='outputs/data_for_heritability.csv', row.names=FALSE)
