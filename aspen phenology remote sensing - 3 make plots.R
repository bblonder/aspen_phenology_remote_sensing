library(dplyr)
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
library(sf)
library(MuMIn)
library(MASS)
library(ggExtra)
library(ggspatial)
library(ggsn)
library(RStoolbox)
library(maps)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(raster)
library(reshape)

#load('outputs/workspace script 2.Rdata')
if(!file.exists('figures'))
{
  dir.create('figures')
}

file_pheno = dir('outputs',pattern="*pheno*",full.names = TRUE)
rasters_pheno = lapply(file_pheno,rast)

nice_names_pheno = c(OGI='Greenup date (doy)', 
                     OGMn='Greendown date (doy)', 
                     GSL.50='Growing season length (days)', 
                     EVImax='Maximum greenness (fraction)')


info_pheno = data.frame(file=file_pheno, 
                        year=gsub("X","",gsub("\\.tif","",sapply(strsplit(file_pheno,'_'),tail,1))),
                        var=sapply(strsplit(file_pheno,'_'),function(x) {x[3]})
                        )
info_pheno$nice_name = nice_names_pheno[info_pheno$var]

# show time series predictors
file_tmax = dir('outputs',pattern="*tmax*",full.names = TRUE)
file_tmax = file_tmax[-c(1)] # start at 2013
rasters_tmax = do.call("c",lapply(file_tmax,rast))

file_sm_q01 = dir('outputs',pattern="*sm_q01*",full.names = TRUE)
file_sm_q01 = file_sm_q01[-c(1)] # start at 2013
rasters_sm_q01 = do.call("c",lapply(file_sm_q01,rast))

file_sm_runs_med_dur = dir('outputs',pattern="*sm_runs_med_dur*",full.names = TRUE)
file_sm_runs_med_dur = file_sm_runs_med_dur[-c(1)] # start at 2013
rasters_sm_runs_med_dur = do.call("c",lapply(file_sm_runs_med_dur,rast))

file_snow = dir('outputs',pattern="*snowmelt*",full.names = TRUE)
file_snow = file_snow[-c(1)] # start at 2013
rasters_snow = do.call("c",lapply(file_snow,rast))





plot_climate <- function(rasters, years=c(2013:2019), varname_nice, colors_diverging=FALSE, qclip=TRUE, na.value='black')
{
  g_final = ggarrange(plotlist = lapply(1:nlyr(rasters), function(yearid) {
    
    if (qclip==TRUE)
    {
      qvals = rasters[] %>% na.omit %>% as.numeric %>% quantile(c(0.02,0.98))
    }
    else
    {
      qvals = rasters[] %>% na.omit %>% as.numeric %>% quantile(c(0.0,1.0))
    }
    
    if(colors_diverging==TRUE)
    {
      qvals = max(abs(qvals)) * c(-1, 1)
    }
    
    g = ggR(raster(rasters[[yearid]]),geom_raster = TRUE,maxpixels=1e7) + 
      theme_bw() +
      theme(panel.grid = element_blank()) +
      coord_equal() + 
      labs(fill=varname_nice) + 
      theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
      xlab("") + 
      ylab("") +
      ggtitle(years[yearid]) +
      theme(plot.margin=margin(0.01, 0.01, 0.01, 0.01, "cm")) +
      theme(panel.background = element_rect(fill = "black"))
    
    if (colors_diverging==FALSE)
    {
      g = g + 
        scale_fill_viridis(na.value=na.value,limits=qvals) 
    }
    else
    {
      g = g +
        scale_fill_gradientn(na.value=na.value,limits=qvals,colours = c('pink','red','orange','lightgray','green','blue','cyan'))
        #scale_fill_distiller(na.value=na.value,limits=qvals,palette='RdBu')
    }
    
    if (yearid==nlyr(rasters))
    {
      g = g + annotation_scale(location = "bl", height = unit(0.1, "cm"),bar_cols=c('darkgray','white'),text_col='white')# +
        #annotation_north_arrow(location = "br",
        #                       style = north_arrow_minimal,
        #                       height = unit(0.5, "cm"))
    }
    return(g)
  }),common.legend = TRUE,legend='top',nrow = 1,ncol=nlyr(rasters))
  return(g_final)
}

r_tmax = plot_climate(rasters_tmax, varname_nice = "Maximum temperature (°C)")
r_sm_q01 = plot_climate(rasters_sm_q01, varname_nice = "Minimum soil moisture (fraction)")
r_sm_runs_max_dur = plot_climate(rasters_sm_runs_med_dur, varname_nice = "Median drought duration (days)")
r_snow = plot_climate(rasters_snow, varname_nice = "Snowmelt date (doy)")

g_env = ggarrange(r_snow, r_tmax, r_sm_q01, r_sm_runs_max_dur, nrow=4,ncol=1,
                    labels='auto',
                    align='hv')
ggsave(g_env, file='figures/g_env.png',width=9,height=9)





# plot phenology
r_pheno_all = lapply(unique(info_pheno$var), function(var) {
  indices = which(info_pheno$var==var)
  
  meanval = mean(unlist(sapply(rasters_pheno[indices], spatSample, size=10000)),na.rm=T)
  print(meanval)
  
  name_printed = gsub("\\(", sprintf("anomaly\n4-year mean=%.2f ", meanval), info_pheno$nice_name[indices[1]])
  name_printed = gsub("\\)", "", name_printed)
  name_printed = gsub("doy", "days", name_printed)
  name_printed = gsub("fraction", "", name_printed)
  
  plots = plot_climate(rast(rasters_pheno[indices]) - meanval,years=2016:2019, 
                       varname_nice=name_printed,
                      colors_diverging = TRUE)
  })

#ggsave(r_pheno_all[[3]],file='figures/g_pheno_greenup_date.png',width=9,height=3.5)
g_pheno_other = ggarrange(r_pheno_all[[1]], r_pheno_all[[2]], r_pheno_all[[3]], r_pheno_all[[4]], nrow=4,ncol=1,
                  labels='auto',
                  align='hv')
ggsave(g_pheno_other, file='figures/g_pheno_other.png',width=9,height=12)





# plot topography
raster_cos_aspect = rast('outputs/r_pred_cos_aspect.tif')
raster_slope = rast('outputs/r_pred_slope.tif')
raster_elevation  = rast('outputs/r_pred_elevation.tif')
raster_height  = terra::mask(rast('outputs/r_pred_height_canopy.tif'), rasters_snow[[1]])


r_t1 = plot_climate(raster_cos_aspect,years=NULL,varname_nice="Cosine aspect",qclip = FALSE)
r_t2 = plot_climate(raster_slope,years=NULL,,varname_nice="Slope (°)",qclip = FALSE)
r_t3 = plot_climate(raster_elevation,years=NULL,,varname_nice="Elevation (m)",qclip = FALSE)
r_t4 = plot_climate(raster_height,years=NULL,,varname_nice="Canopy height (m)",qclip = FALSE)


raster_cytotype = rast('outputs/r_pred_is_diploid.tif')
r_cytotype = plot_climate(raster_cytotype, years=NULL, varname_nice="Cytotype (fraction diploid)", qclip=FALSE)


g_topo = ggarrange(r_t1, r_t2, r_t3, r_t4, r_cytotype, align='hv',nrow=2,ncol=3, legend='top',labels='auto')
ggsave(g_topo, file='figures/g_topo_genetics.png',width=9,height=6)










# make a map
s_ws = shapefile('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/spectra analysis neon aop/cytotype analysis/watershed.shp')
s_ws_df = fortify(s_ws)

df_plots = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/SI/aspen data tree-level processed 30 Mar 2020.csv')

g_map_zoom = ggplot() + 
  geom_polygon(data=s_ws_df,aes(x=long,y=lat,group=group),color='purple',fill='white') +
  geom_point(data=df_plots,aes(x=X.UTM,y=Y.UTM),alpha=0.03) +
  theme_bw() + 
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
  xlab("") + ylab("") +
  coord_equal() +
  annotation_scale(location = "bl", height = unit(0.1, "cm")) +
  annotation_north_arrow(location = "br", 
                         style = north_arrow_minimal, 
                         height = unit(0.5, "cm"))



s_ws_ll = spTransform(s_ws,  "+init=epsg:4121 +proj=longlat +ellps=GRS80")
s_ws_ll_df = fortify(s_ws_ll)

USA <- map_data("state")

aspen = shapefile('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2020/aspen site selection for roots/Little range map/poputrem.shp')
crs(aspen) = "+init=epsg:4121 +proj=longlat +ellps=GRS80"
aspen_df = broom::tidy(aspen)

ws_df = broom::tidy(ws_boundaries)

data_cities = us.cities %>% 
  filter(country.etc %in% "CO") %>%
  mutate(name=gsub(" CO","",name)) %>%
  filter(name %in% c("Denver"))

g_map_big <- ggplot() + geom_polygon(data = USA, 
                                 aes(x=long, y = lat, group = group), 
                                 fill = 'white', 
                                 color='gray') +
  geom_polygon(data = aspen_df %>% filter(hole==FALSE), aes(x=long, y=lat, group=group),alpha=0.9,color=NA,fill='lightgray') +
  geom_polygon(data = aspen_df %>% filter(hole==TRUE), aes(x=long, y=lat, group=group),color=NA,fill='white') +
  geom_text_repel(data=data_cities,aes(x=long,y=lat,label=name),
                   size = 3,
                   nudge_x      = 0.15,
                   direction    = "y",
                   hjust        = 0,
                   alpha=0.5) +
  geom_polygon(data = s_ws_ll_df, aes(x=long, y=lat, group=group),fill='purple',col='purple',size=0.01) +
  geom_point(data=data_cities,aes(x=long,y=lat),size=0.2) +
  theme_bw() +
  coord_equal(ylim=c(34.5,42),xlim = c(-109,-102)) +
  xlab("Longitude (°)") +
  ylab("Latitude (°)")

g_map = ggarrange(g_map_big, g_map_zoom,align='hv',labels = 'auto')
ggsave(g_map, file='figures/g_map.png',width=8,height=5)





# get some distributions
df_all = read.csv('outputs/df_all_aspen_cover_0.25.csv')

df_all_for_plotting = df_all %>%
  filter(aspen_cover >= 0.5 & cytotype_fraction_diploid >= 0) %>%
  #group_by(year, cytotype) %>%
  dplyr::select(year, cytotype_fraction_diploid, OGI, OGMn, GSL.50,EVImax) %>%
  melt(id.vars=c("year","cytotype_fraction_diploid")) %>%
  mutate(variable=as.character(variable))

# make a table
df_all_for_plotting %>% group_by(variable) %>%
  summarize(q10=quantile(value,0.025),q50=quantile(value,0.50),q90=quantile(value,0.975))


g_range_cyto = ggplot(df_all_for_plotting,aes(x=factor(year),y=value,
                               color=cut(cytotype_fraction_diploid,breaks=seq(0,1,by=0.25)),
                               fill=cut(cytotype_fraction_diploid,breaks=seq(0,1,by=0.25)))) + 
  geom_violin(draw_quantiles=c(0.1,0.5,0.9),alpha=0.5) + 
  facet_wrap(~variable,scales='free',labeller = as_labeller(nice_names_pheno)) +
  theme_bw() +
  scale_color_viridis(name='Cytotype diploid fraction',discrete=TRUE, na.translate=FALSE) +
  scale_fill_viridis(name='Cytotype diploid fraction',discrete=TRUE, na.translate=FALSE) +
  ylab('Value') +
  xlab('Year')





df_aspen_sex = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/aspen sex markers/aspen_sex_aug_11_2021.csv')
df_aspen_site = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/aspen data site-level processed 30 Mar 2020.csv')

df_aspen_joined = df_aspen_site %>% dplyr::select(Site_Code, X.UTM, Y.UTM, Cytotype=Ploidy_level) %>%
  left_join(df_aspen_sex %>% dplyr::select(Site_Code, geneticSexID))  

paths_pheno = dir(path = 'outputs', pattern='r_pheno*',full.names = TRUE)
names_vars = gsub('\\.tif','',gsub("r_pheno_","",basename(paths_pheno)))
rasters_pheno = rast(paths_pheno)
names(rasters_pheno) = names_vars

raster_cover = rast('outputs/r_pred_aspen_cover.tif')



vals_pheno = extract(c(raster_cover,rasters_pheno), df_aspen_joined %>% dplyr::select(X.UTM,Y.UTM))

df_aspen_joined_final = cbind(df_aspen_joined, vals_pheno) %>% dplyr::select(-ID,-X.UTM,-Y.UTM)

df_aspen_joined_final_melted = melt(df_aspen_joined_final, 
                                    id.vars=c('Site_Code','Cytotype','geneticSexID','aspen_cover')) %>%
  mutate(var_pheno = sapply(strsplit(as.character(variable),'_'),head,1)) %>%
  mutate(year = gsub("X","",sapply(strsplit(as.character(variable),'_'),tail,1))) %>%
  arrange(Site_Code,var_pheno,year)

# get final dataset with viable levels of cover
df_for_plotting = df_aspen_joined_final_melted %>% filter(aspen_cover >= 0.25) %>% na.omit

# count # of sites at 25% cover
df_for_plotting$Site_Code %>% unique %>% length
nrow(df_for_plotting)



nice_names_pheno_ground_based = c(OGI='Greenup date (doy)', 
                                  OGMn='Greendown date (doy)', 
                                  GSL.50='Growing season length (days)', 
                                  EVImax='Maximum greenness (fraction)',
                                  `2016`=2016,
                                  `2017`=2017,
                                  `2018`=2018,
                                  `2019`=2019)


g_range_sex = ggplot(df_for_plotting, aes(color=Cytotype,fill=Cytotype,linetype=geneticSexID,x=year,y=value)) +
  facet_wrap(~var_pheno,scales='free',labeller = as_labeller(nice_names_pheno_ground_based)) +
  geom_violin(draw_quantiles=c(0.5),alpha=0.5) +
  theme_bw() +
  scale_color_viridis_d(name='Cytotype',direction=-1) +
  scale_fill_viridis_d(name='Cytotype',direction=-1) +
  scale_linetype_manual(values=c('solid','dotted'),name='Sex') +
  ylab('Value') +
  xlab('Year')



ggsave(ggarrange(g_range_cyto, g_range_sex, nrow=2,ncol=1,labels='auto',align='hv'), 
       file='figures/g_range.png',width=7,height=8)




# draw the PDPs

nice_names_preds = c(`year`="Year",
                     `cytotype_fraction_diploid`='Cytotype (fraction diploid)',
                     `elevation`='Elevation (m)', 
                     `cos_aspect`='Cosine aspect',  
                     `slope`='Slope (°)', 
                     `height_canopy`='Canopy height (m)',
                     `tmax_q99.t.minus.3`='Maximum temperature (t-3) (°C)',
                     `tmax_q99.t.minus.2`='Maximum temperature (t-2) (°C)', 
                     `tmax_q99.t.minus.1`='Maximum temperature (t-1) (°C)', 
                     `tmax_q99.t.minus.0`='Maximum temperature (t-0) (°C)', 
                     `snowmelt.t.minus.3`='Snowmelt date (t-3) (doy)',
                     `snowmelt.t.minus.2`='Snowmelt date (t-2) (doy)', 
                     `snowmelt.t.minus.1`='Snowmelt date (t-1) (doy)', 
                     `snowmelt.t.minus.0`='Snowmelt date (t-0) (doy)', 
                     `sm_runs_med_dur.t.minus.3`='Median drought duration (t-3) (days)', 
                     `sm_runs_med_dur.t.minus.2`='Median drought duration (t-2) (days)', 
                     `sm_runs_med_dur.t.minus.1`='Median drought duration (t-1) (days)', 
                     `sm_runs_med_dur.t.minus.0`='Median drought duration (t-0) (days)',
                     `sm_q01.t.minus.3`='Minimum soil moisture (t-3) (g/g)', 
                     `sm_q01.t.minus.2`='Minimum soil moisture (t-2) (g/g)', 
                     `sm_q01.t.minus.1`='Minimum soil moisture (t-1) (g/g)', 
                     `sm_q01.t.minus.0`='Minimum soil moisture (t-0) (g/g)'
                     )

# summarize each phenology response variable
rf_pdp_summaries = read.csv('outputs/rf_pdp_summaries.csv')

names_short = unlist(lapply(strwrap(c(nice_names_pheno, nice_names_preds), width=15, simplify=FALSE), paste, 
                           collapse="\n"))
names(names_short) = names(c(nice_names_pheno, nice_names_preds))


rf_importances = read.csv('outputs/rf_importances.csv')
rf_relative_importances = rf_importances %>% 
  group_by(yvar) %>% 
  mutate(importance = importance/max(importance)) %>% 
  ungroup %>%
  group_by(yvar, variable) %>% 
  summarize(relimp=mean(importance)) %>%
  ungroup %>%
  mutate(combo=paste(yvar, variable)) %>%
  dplyr::select(yvar, combo, relimp)

rf_relative_importances_top_n = rf_relative_importances %>%
  mutate(contains_year = grepl("year",combo)) %>%
  filter(contains_year==FALSE) %>%
  ungroup %>%
  group_by(yvar) %>%
  arrange(yvar, relimp) %>% 
  slice(tail(row_number(), 5))

rf_pdp_summaries_joined = rf_pdp_summaries %>% 
  mutate(combo=paste(yvar, xvar)) %>%
  left_join(rf_relative_importances %>% dplyr::select(-yvar), by='combo')

g_pdps = ggplot(rf_pdp_summaries_joined, aes(x=x,
                             y=yhat.q50,
                             ymin=yhat.q05,
                             ymax=yhat.q95,
                             color=factor(year),
                             group=factor(year),
                             alpha=relimp)) +
  geom_line(size=2) +
  geom_errorbar() +
  theme_bw() +
  facet_grid(yvar~xvar,scales='free',labeller = as_labeller(names_short)) +
  scale_color_brewer(palette='Spectral',name='Year') +
  scale_alpha(name='Mean relative\nvariable importance') +
  xlab("Predictor") +
  ylab("Response") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(g_pdps, file='figures/g_pdps.png',width=22,height=8)



rf_pdp_summaries_joined_top_n = rf_pdp_summaries %>% 
  ungroup %>%
  mutate(combo=paste(yvar, xvar)) %>%
  right_join(rf_relative_importances_top_n %>% 
               ungroup %>% 
               dplyr::select(-yvar) %>% as.data.frame, by='combo')

g_pdp_summaries_top_n_raw = lapply(rf_pdp_summaries_joined_top_n %>% group_by(yvar) %>% group_split, function(df) {
  xvar_ordered = df %>% 
    mutate(relimp_order=as.numeric(factor(relimp))) %>% 
    dplyr::select(xvar, relimp, relimp_order) %>% 
    unique %>% 
    arrange(relimp_order)
  print(xvar_ordered)
  
  df_ordered = df %>% 
    mutate(xvar = factor(xvar, levels=rev(xvar_ordered$xvar), ordered=TRUE))
  
  ggplot(df_ordered, aes(x,
                 y=yhat.q50,
                 ymin=yhat.q05,
                 ymax=yhat.q95,
                 color=factor(year),
                 group=factor(year))) +
    geom_line(size=2,alpha=0.75) +
    geom_errorbar() +
    theme_bw() +
    facet_wrap(~xvar,scales='free_x',labeller = as_labeller(names_short),nrow=1) +   
    scale_color_brewer(palette='Spectral',name='Year') +
    #scale_alpha(name='Mean relative\nvariable importance') +
    xlab("Predictor") +
    ylab("Response") +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(nice_names_pheno[df$yvar[1]])
})

g_pdp_summaries_top_n_all = ggarrange(plotlist=g_pdp_summaries_top_n_raw, 
                                      ncol=1, common.legend = TRUE,
                                      legend='bottom',
                                      labels='auto')

ggsave(g_pdp_summaries_top_n_all,file='figures/g_pdp_summaries_top_n_all.png',width=9,height=10)









g_imp = ggplot(rf_importances, aes(x=variable,y=importance)) +
  facet_wrap(~yvar,scales='free_x',labeller = as_labeller(nice_names_pheno),nrow=1,ncol=4) + 
  geom_boxplot() +
  theme_bw() +
  xlab("Predictor") +
  ylab("Variable importance") +
  coord_flip() +
  scale_x_discrete(labels=nice_names_preds)

ggsave(g_imp, file='figures/g_imp.png',width=11,height=5)


rf_r2 = read.csv('outputs/rf_r2.csv')

# summarize
rf_r2 %>% group_by(yvar) %>% summarize(mean(r2),sd(r2))

g_r2 = ggplot(rf_r2, aes(x=yvar,y=r2)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Response variable") +
  ylab(expression(paste("Random forest model ",R^2))) +
  ylim(0,1) +
  scale_x_discrete(labels=nice_names_pheno) +
  coord_flip()
ggsave(g_r2, file='figures/g_r2.png',width=6,height=2)
