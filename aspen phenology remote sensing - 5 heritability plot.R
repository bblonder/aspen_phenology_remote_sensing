library(tidyverse)
library(ggplot2)
library(viridis)

# load data
df_heritability <- read_tsv("heritability_paper_data.tsv") %>%
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
    geom_errorbar(aes(ymin=(Variance-SE), ymax=(Variance+SE)),
                  position=position_dodge(width=0.4)) +
    facet_wrap(~phenology, labeller = as_labeller(nice_names_pheno_short_with_labels),nrow=1) +
    theme(axis.text.x = element_text(angle = 75, hjust=1)) + 
    labs(y="Vg/Vp (Field heritability)", x='Cytotype') +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() +
    theme(legend.position='none') +
    scale_color_manual(values=c(gray(0.3),viridis(n=3)[3],viridis(n=3)[1])) +
    scale_x_discrete(drop=FALSE)
}

#save plot
ggsave(make_plot(df_heritability %>% filter(cover=="cover25")), file="figures/g_vgvp_25.png", 
       width=7, height=4)
ggsave(make_plot(df_heritability %>% filter(cover=="cover50")), file="figures/g_vgvp_50.png", 
       width=7, height=4)
