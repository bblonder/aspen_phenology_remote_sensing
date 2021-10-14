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
library(dplyr)
library(ggplot2)
library(MuMIn)
library(ggpubr)
library(MASS)
library(visreg)
library(ggExtra)
library(ggspatial)
library(ggsn)
library(RStoolbox)
library(wesanderson)
library(maps)
library(ggrepel)
library(rgdal)
library(mgcv)
library(terra)
library(RColorBrewer)
library(viridis)
library(pals)
library(mgcViz)
library(raster)
library(reshape)
library(yarrr)

file_pheno = dir('rasters_for_plotting',pattern="*pheno*",full.names = TRUE)
rasters_pheno = lapply(file_pheno,rast)

plot_pheno <- function(phenoid, varname, varname_nice)
{
  g_final = ggarrange(plotlist = lapply(1:4, function(yearid) {
    years = c(2016:2019)
    
    qvals = rasters_pheno[[phenoid]][] %>% na.omit %>% as.numeric %>% quantile(c(0.02,0.98))
    
    r_this = raster::stack(rasters_pheno[[phenoid]])
    
    g = ggR(r_this[[yearid]],geom_raster = TRUE) + 
      scale_fill_viridis(na.value='gray',limits=qvals) + 
      theme_bw() +
      coord_equal() + 
      labs(fill=varname_nice) + 
      theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
      xlab("") + 
      ylab("") +
      ggtitle(years[yearid])
    
    if (yearid==4)
    {
      g = g + annotation_scale(location = "bl", height = unit(0.1, "cm")) +
        annotation_north_arrow(location = "br",
                               style = north_arrow_minimal,
                               height = unit(0.5, "cm"))
    }
    return(g)
  }),common.legend = TRUE,legend='top',nrow = 1,ncol=4)
  return(g_final)
}

g_pheno_1 = plot_pheno(1, 'EVImax','Maximum greenness')
g_pheno_2 = plot_pheno(2, 'GSL.50','Growing season length (days)')
g_pheno_3 = plot_pheno(3, 'OGI','Greenup date (doy)')
g_pheno_4 = plot_pheno(4, 'OGMn','Greendown date (doy)')

g_pheno = ggarrange(g_pheno_1, g_pheno_2, g_pheno_3, g_pheno_4,nrow=4,ncol=1,
                    labels='auto',
                    align='hv')
ggsave(g_pheno, file='g_pheno.png',width=8,height=12)


# show time series predictors
file_tmax = dir('rasters_for_plotting',pattern="*tmax*",full.names = TRUE)
file_tmax = file_tmax[-c(1,2)] # start at 2014
rasters_tmax = do.call("c",lapply(file_tmax,rast))

file_sm = dir('rasters_for_plotting',pattern="*sm_q01*",full.names = TRUE)
file_sm = file_sm[-c(1,2)] # start at 2014
rasters_sm = do.call("c",lapply(file_sm,rast))

file_snow = dir('rasters_for_plotting',pattern="*snowmelt*",full.names = TRUE)
file_snow = file_snow[-c(1,2)] # start at 2014
rasters_snow = do.call("c",lapply(file_snow,rast))





plot_climate <- function(rasters, years=c(2014:2019), varname, varname_nice)
{
  g_final = ggarrange(plotlist = lapply(1:nlyr(rasters), function(yearid) {
    
    qvals = rasters[] %>% na.omit %>% as.numeric %>% quantile(c(0.02,0.98))
    
    g = ggR(raster(rasters[[yearid]]),geom_raster = TRUE) + 
      scale_fill_viridis(na.value='gray',limits=qvals, option = 'plasma') + 
      theme_bw() +
      coord_equal() + 
      labs(fill=varname_nice) + 
      theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
      xlab("") + 
      ylab("") +
      ggtitle(years[yearid])
    
    if (yearid==nlyr(rasters))
    {
      g = g + annotation_scale(location = "bl", height = unit(0.1, "cm")) +
        annotation_north_arrow(location = "br",
                               style = north_arrow_minimal,
                               height = unit(0.5, "cm"))
    }
    return(g)
  }),common.legend = TRUE,legend='top',nrow = 1,ncol=nlyr(rasters))
  return(g_final)
}

r_tmax = plot_climate(rasters_tmax / 10, varname="tmax_q99",varname_nice = "Maximum temperature (°C)")
r_sm = plot_climate(rasters_sm, varname="sm_q01",varname_nice = "Minimum soil moisture (fraction)")
r_snow = plot_climate(rasters_snow, varname="snowmelt",varname_nice = "Snowmelt date (doy)")

g_env = ggarrange(r_tmax, r_sm, r_snow, nrow=3,ncol=1,
                    labels='auto',
                    align='hv')
ggsave(g_env, file='g_env.png',width=8,height=8)




# plot topography
raster_cos_aspect = rast('rasters_for_plotting/r_pred_cos_aspect.tif')
raster_slope = rast('rasters_for_plotting/r_pred_slope.tif')
raster_elevation  = rast('rasters_for_plotting/r_pred_elevation.tif')
raster_height  = terra::mask(rast('rasters_for_plotting/r_pred_height_canopy.tif'), rasters_snow[[1]])


plot_topo <- function(r, varname, varname_nice, add_annotation=FALSE, quantile_cut = FALSE)
{
  if (quantile_cut==TRUE)
  {
    qvals = r[] %>% na.omit %>% as.numeric %>% quantile(c(0.02,0.98))
  }
  else
  {
    qvals = range(r[] %>% na.omit %>% as.numeric)
  }
    
  g = ggR(raster(r),geom_raster = TRUE) + 
    scale_fill_viridis(na.value='gray',limits=qvals, option = 'cividis') + 
    theme_bw() +
    coord_equal() + 
    labs(fill=varname_nice) + 
    theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
    xlab("") + 
    ylab("")
  
  if (add_annotation==TRUE)
  {
    g = g + annotation_scale(location = "bl", height = unit(0.1, "cm")) +
      annotation_north_arrow(location = "br",
                             style = north_arrow_minimal,
                             height = unit(0.5, "cm"))
  }
  
  return(g)
}

r_t1 = plot_topo(raster_cos_aspect,"cos_aspect","Cosine aspect")
r_t2 = plot_topo(raster_slope,"slope","Slope (°)")
r_t3 = plot_topo(raster_elevation,"elevation","Elevation (m)")
r_t4 = plot_topo(raster_height,"height","Canopy height (m)")

g_topo = ggarrange(r_t1, r_t2, r_t3, r_t4,align='hv',nrow=2,ncol=2, legend='bottom',labels='auto')
ggsave(g_topo, file='g_topo.png',width=6,height=8)



plot_genetics <- function(r, varname, varname_nice, add_annotation=FALSE, col_low, col_high)
{
  g = ggR(raster(r),geom_raster = TRUE) + 
    scale_fill_gradient2(low=col_low,high=col_high,mid='gray',midpoint=0.5) + 
    theme_bw() +
    coord_equal() + 
    labs(fill=varname_nice) + 
    theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank()) +
    xlab("") + 
    ylab("")
  
  if (add_annotation==TRUE)
  {
    g = g + annotation_scale(location = "bl", height = unit(0.1, "cm")) +
      annotation_north_arrow(location = "br",
                             style = north_arrow_minimal,
                             height = unit(0.5, "cm"))
  }
  
  return(g)
}

raster_cytotype = rast('rasters_for_plotting/r_pred_is_diploid.tif')
raster_sex = rast('rasters_for_plotting/r_pred_is_diploid.tif'); warning('need to replace with real map when phil makes it')



r_cytotype = plot_genetics(raster_cytotype, 'cytotype', "Cytotype (fraction diploid)", col_low = 'red',col_high = 'blue')
r_sex = plot_genetics(raster_sex, 'sex', "Sex (fraction female)", col_low = 'purple',col_high = 'orange')

g_genetics = ggarrange(r_cytotype, r_sex,
          align='hv',nrow=1,ncol=2, legend='bottom',labels='auto')
          
ggsave(g_genetics, file='g_genetics.png',width=6,height=4)





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
                   label.padding=0.05,
                   alpha=0.5) +
  geom_polygon(data = s_ws_ll_df, aes(x=long, y=lat, group=group),fill='purple',col='purple',size=0.01) +
  geom_point(data=data_cities,aes(x=long,y=lat),size=0.2) +
  theme_bw() +
  coord_equal(ylim=c(34.5,42),xlim = c(-109,-102)) +
  xlab("Longitude (°)") +
  ylab("Latitude (°)")

g_map = ggarrange(g_map_big, g_map_zoom,align='hv',labels = 'auto')
ggsave(g_map, file='g_map.png',width=8,height=5)





# get some distributions
df_all = read.csv('df_all.csv')

df_all_for_plotting = df_all %>%
  #group_by(year, cytotype) %>%
  dplyr::select(year, cytotype, OGI, OGMn, GSL.50,EVImax) %>%
  melt(id.vars=c("year","cytotype")) %>%
  group_split(variable)

nice_names_pheno = c(OGI='Greenup date (doy)', OGMn='Greendown date (doy)', GSL.50='Growing season length (days)', EVImax='Maximum greenness')

g_by_year = ggarrange(plotlist=lapply(df_all_for_plotting, function(df_ss) {
  ggplot(df_ss, aes(x=value,col=cytotype)) +
    geom_density() +
    facet_wrap(~year,scales='free_y',nrow=4,ncol=1) +
    xlab(nice_names_pheno[df_ss$variable[1]]) +
    scale_color_manual(values=c("blue","red"),name='Cytotype') +
    theme_bw() +
    ylab("Density") +
    theme(axis.text.y = element_blank())
}),nrow=1,ncol=4, common.legend = TRUE,legend='bottom')
ggsave(g_by_year,file='g_by_year.png',width=8,height=6)





# draw the PDPs
rm(list=ls())
load('workspace all.Rdata')

nice_names_pheno = c(OGI='Greenup date (doy)', OGMn='Greendown date (doy)', GSL.50='Growing season length (days)', EVImax='Maximum greenness')
nice_names_preds = c(`elevation`='Elevation (m)', 
                     `cos_aspect`='Cosine aspect',  
                     `slope`='Slope (°)', 
                     `height_canopy`='Canopy height (m)',
                     `tmax.q99.t.minus.2`='Maximum temperature (t-2) (°C)', 
                     `tmax.q99.t.minus.1`='Maximum temperature (t-1) (°C)', 
                     `tmax.q99.t.minus.0`='Maximum temperature (t-0) (°C)', 
                     `snowmelt.date.t.minus.2`='Snowmelt date (t-2) (doy)', 
                     `snowmelt.date.t.minus.1`='Snowmelt date (t-1) (doy)', 
                     `snowmelt.date.t.minus.0`='Snowmelt date (t-0) (doy)', 
                     `soilmoisture0.1.q01.t.minus.2`='Minimum soil moisture (t-2) (g/g)', 
                     `soilmoisture0.1.q01.t.minus.1`='Minimum soil moisture (t-1) (g/g)', 
                     `soilmoisture0.1.q01.t.minus.0`='Minimum soil moisture (t-0) (g/g)', 
                     `soilmoisture0.1.run.med.t.minus.2`='DELETE2', 
                     `soilmoisture0.1.run.med.t.minus.1`='DELETE1', 
                     `soilmoisture0.1.run.med.t.minus.0`='DELETE0')

# summarize each phenology response variable


make_all_pdps_1d <- function(model_list)
{
  xvars = setdiff(model_list$xvar,c("year","cytotype"))
  
  pdps = lapply(xvars, function(xvar) {
    print(xvar)
    
    result = do_pdp(
      model_list = model_list, 
      pred.vars = c(xvar,"year","cytotype"),
      grid.resolution=5,
      categorical = FALSE,
      df_train = df_all_for_rf) 
    
    result$yvar = model_list$result$yvar[1]
    
    return(result)
  })
  names(pdps) = xvars
  
  return(pdps)
}

plot_pdp_1d <- function(df_pdp_summary, ylim)
{
  xvar = names(df_pdp_summary)[1]
  names(df_pdp_summary)[1] = "xvar"
  
  yvar = df_pdp_summary$yvar[1]
  
  df_pdp_summary$year = factor(df_pdp_summary$year)
  
  g = ggplot(df_pdp_summary, aes(x=xvar,
                                 y=yhat.q50,
                                 ymin=yhat.q05,
                                 ymax=yhat.q95,
                                 col=year,
                                 linetype=cytotype)) +
    geom_point() +
    geom_line() +
    geom_errorbar() +
    theme_bw() + 
    xlab(nice_names_preds[xvar]) +
    ylab(nice_names_pheno[yvar]) +
    ylim(ylim) +
    facet_wrap(~cytotype) +
    scale_color_manual(values=piratepal(palette = 'xmen')[c(2,4,3,1)] %>% as.character,name='Year') +
    scale_linetype(name='Cytotype')
  
  return(g)
}


plot_pdp_1d_all <- function(pdps, name)
{
  pdps_summaries = lapply(pdps, summarize_pdp)
  
  ylim_all = range(as.numeric(sapply(pdps_summaries, function(x) {x$yhat.q50})))
  
  plots_pdp_1d_all = ggarrange(plotlist = lapply(pdps_summaries, plot_pdp_1d, ylim=ylim_all),
                               common.legend = TRUE,
                               legend='bottom')
  
  ggsave(plots_pdp_1d_all, file=sprintf('g_pdp_1d_all_%s.png',name),width=12,height=12)
}

# write out the plots
lapply(1:length(pdps_all), function(i) { 
  print(i)
  plot_pdp_1d_all(pdps_all[[i]], names(pdps_all)[i])
})



# make r2 plots
rf_r2 = t(do.call("rbind",lapply(models_all, function(x) {x$result$r2})))

apply(rf_r2, 2, function(x) {c(mean=mean(x), sd=sd(x))})

rf_r2_melted = melt(rf_r2)
rf_r2_melted$nicename = sapply(strsplit(nice_names_pheno[rf_r2_melted$X2], " \\("), head, 1)
g_r2 = ggplot(rf_r2_melted, aes(x=nicename,y=value)) + 
  geom_boxplot() +
  ylim(0,1) +
  theme_bw() +
  ylab(expression(paste("Cross-validated ", R^2))) +
  xlab("Response variable") 

ggsave(g_r2,file='g_r2.png',width=6,height=5)
