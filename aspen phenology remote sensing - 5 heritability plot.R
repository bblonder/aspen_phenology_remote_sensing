library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# load data
df_heritability <- read_tsv("paper_data_jan_2023.tsv") %>%
  mutate(ploidy=factor(ploidy))

nice_names_pheno = c(OGI='Greenup date (doy)', 
                     OGMn='Greendown date (doy)', 
                     GSL='Growing season length (days)', 
                     EVImax='Maximum greenness (fraction)')

nice_names_pheno_short = sapply(nice_names_pheno, function(x) {  strsplit(x,split=" \\(")[[1]][1]  })
nice_names_pheno_short = nice_names_pheno_short[order(names(nice_names_pheno_short))]
nice_names_pheno_short_with_labels = paste("(",letters[1:4], ") ", nice_names_pheno_short,sep="")
names(nice_names_pheno_short_with_labels) = names(nice_names_pheno_short)

#make plot
make_plot <- function(df)
{
  plot_vgvp <- ggplot(df, aes(x=ploidy, y=Variance,color=ploidy)) +
    geom_point(size=3, position=position_dodge(width=0.5)) +
    geom_errorbar(aes(ymin=lower_bound, ymax=upper_bound),
                  position=position_dodge(width=0.25)) +
    facet_wrap(~phenology, labeller = as_labeller(nice_names_pheno_short_with_labels),nrow=1) +
    theme(axis.text.x = element_text(angle = 75, hjust=1)) + 
    labs(y="Heritability", x='Cytotype') +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() +
    theme(legend.position='none') +
    scale_color_manual(values=c(`all`='gray',`diploid`=brewer.pal(n=3,name='BrBG')[3],`triploid`=brewer.pal(n=3,name='BrBG')[1])) +
    scale_x_discrete(drop=FALSE) +
    geom_text(aes(x=as.numeric(ploidy)+0.25,y=Variance,label=n),color='black',size=2)
}

#save plot
lapply(paste("cover",c(25,30,50,70),sep=""), function(cover_this)
{
  cat('.')
  ggsave(make_plot(df_heritability %>% filter(cover==cover_this)), file=sprintf("figures/g_vgvp_%s.png",cover_this), 
         width=7, height=4)
})
