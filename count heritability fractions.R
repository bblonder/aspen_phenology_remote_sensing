data = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/aspen phenology/aspen phenology remote sensing analysis/aspen_phenology_remote_sensing/outputs/data_for_heritability.csv')

library(dplyr)
values_na = data %>% select(contains("pheno")) %>% as.matrix %>% as.numeric %>% is.na

length(which(!values_na))/length(values)

data2 = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/aspen phenology/aspen phenology remote sensing analysis/aspen_phenology_remote_sensing/outputs/df_all_aspen_cover_0.25.csv')
data2 %>% filter(aspen_cover>= 0.5) %>% nrow
