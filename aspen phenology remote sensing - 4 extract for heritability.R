library(dplyr)
library(terra)

r_phenology_2016 = rast('../../13SCD/MSLSP_13SCD_2016.nc')
r_phenology_2017 = rast('../../13SCD/MSLSP_13SCD_2017.nc')
r_phenology_2018 = rast('../../13SCD/MSLSP_13SCD_2018.nc')
r_phenology_2019 = rast('../../13SCD/MSLSP_13SCD_2019.nc')

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
  select(Site_Code, x=X.UTM, y=Y.UTM, Elevation, Cos.aspect, Slope, Summer.Insolation, DBH.mean, Canopy_openness, Soil.type, Geologic.Unit, Cytotype = Ploidy_level)

df_sex = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/aspen sex markers/aspen_sex_aug_11_2021.csv') %>%
  select(Site_Code, geneticSexID)

df_site = df_site %>%
  left_join(df_sex, by="Site_Code")

# simplify the geology
df_site$Rock_Unit = factor(sapply(df_site$Geologic.Unit,simplify_geology))

df_site = df_site %>% select(-Geologic.Unit)


extract_pheno_data <- function(r, year)
{
  e_pheno = terra::extract(x=r, y=df_site %>% select(x,y)) %>%
    select("OGI","50PCGI","50PCGD","OGMn","EVImax","gupQA","gdownQA")
  
  e_pheno$OGI[e_pheno$gupQA > 4] = NA
  e_pheno$`50PCGI`[e_pheno$gupQA > 4] = NA
  e_pheno$OGMn[e_pheno$gdownQA > 4] = NA
  e_pheno$`50PCGD`[e_pheno$gdownQA > 4] = NA
  e_pheno$GSL = e_pheno$`50PCGD` - e_pheno$`50PCGI`
  
  names(e_pheno) = paste("pheno",year,names(e_pheno),sep=".")
  
  return(e_pheno)
}

df_2016 = extract_pheno_data(r_phenology_2016, 2016)
df_2017 = extract_pheno_data(r_phenology_2017, 2017)
df_2018 = extract_pheno_data(r_phenology_2018, 2018)
df_2019 = extract_pheno_data(r_phenology_2019, 2019)


df_site_combined = cbind(df_site, df_2016, df_2017, df_2018, df_2019)

df_site_combined_trimmed = df_site_combined %>%
  select(-contains("QA"),-contains("50"))

# get fractional cover
r_cytotype_masked = rast('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/r_cytotype_masked.tif') 
r_aspen_cover = !is.na(r_cytotype_masked)
r_aspen_cover_projected = project(r_aspen_cover, r_phenology_2016)
r_aspen_cover_trimmed = crop(r_aspen_cover_projected, r_cytotype_masked)

df_site_combined_trimmed$fraction_aspen = terra::extract(r_aspen_cover_trimmed, y=df_site_combined_trimmed %>% select(x,y))[,2]


# write out file
write.csv(df_site_combined_trimmed, file='outputs/data_for_heritability.csv', row.names=FALSE)
