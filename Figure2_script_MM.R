###LOAD LIBRARIES
{
  library('vegan')
  library('ggplot2')
  library('dplyr')
  library('reshape2')
  library('RVAideMemoire')
  library('ggpubr')
  library('Hmisc')
  library('tidyverse')
  library('scales')
  library(psych)
  library(biotools)
  library(ez)
  library(afex) 
  library(emmeans)
  library(rstatix)
  library("devtools")
  library("grid")
  library("ComplexHeatmap")
  library("circlize")
  library("cluster")
  library("dendextend")
  library("matrixStats")
  library(rstatix)
  library(Maaslin2)
  library('ggrepel')
  library(egg)
  library(magick)
  library('patchwork')
  library(EBImage)
}



#Figure 2 (g)


shares <- read.csv("sharing_events_input.csv")

library('ggplot2')

# Bar chart side by side
shares <- within(shares, 
                   Species <- factor(Species, 
                                      levels=names(sort(table(Species), 
                                                        increasing=TRUE))))

shares$Method = factor(shares$Method, levels=c('Culture Isolation','Mixed (inStrain and Culture isolation)', 'inStrain',  'StrainPhlAn3 (Refined Threshold)','Combined (All Methods)'))

shares$Species = factor(shares$Species, levels=c('Veillonella parvula',
                                                 'Sutterella wadsworthensis',
                                                 'Streptococcus salivarius',
                                                 'Staphylococcus epidermidis',
                                                 'Ruthenibacterium lactatiformans',
                                                 'Ruminococcus gnavus',
                                                 'Prevotella copri',
                                                 'Phascolarctobacterium faecium',
                                                 'Parabacteroides distasonis',
                                                 'Megasphaera sp DISK 18',
                                                 'Megasphaera massiliensis',
                                                 'Lactobacillus paragasseri',
                                                 'Lactobacillus gasseri',
                                                 'Lactobacillus crispatus',
                                                 'GGB32884 SGB58188',
                                                 'Fusicatenibacter saccharivorans',
                                                 'Escherichia virus P1',
                                                 'Escherichia coli',
                                                 'Enterobacter asburiae',
                                                 'Collinsella aerofaciens',
                                                 'Clostridium clostridioforme',
                                                 'Barnesiella intestinihominis',
                                                 'Alistipes onderdonkii',
                                                 'Akkermansia muciniphila',
                                                 'Bifidobacterium pseudocatenulatum',
                                                 'Bifidobacterium longum',
                                                 'Bifidobacterium infantis',
                                                 'Bifidobacterium catenulatum',
                                                 'Bifidobacterium breve',
                                                 'Bifidobacterium bifidum',
                                                 'Bifidobacterium adolescentis',
                                                 'Bacteroides vulgatus',
                                                 'Bacteroides uniformis',
                                                 'Bacteroides stercoris',
                                                 'Bacteroides plebeius',
                                                 'Bacteroides massiliensis',
                                                 'Bacteroides intestinalis',
                                                 'Bacteroides fragilis',
                                                 'Bacteroides faecis',
                                                 'Bacteroides eggerthii',
                                                 'Bacteroides dorei',
                                                 'Bacteroides coprophilus',
                                                 'Bacteroides coprocola',
                                                 'Bacteroides clarus',
                                                 'Bacteroides cellulosilyticus',
                                                 'Bacteroides caccae'))




tight<- ggplot(shares, aes(x = Species, fill = Method, line = "black", border = "black")) +
  geom_bar(colour="black") +
  facet_grid(. ~ Method, scales = "fixed") +
  theme(axis.text.y = element_text(face = "italic", vjust = 0.2),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.3)) +
  scale_fill_manual(values=c('#b35806','#e08214','#fdb863','#fee0b6','#d8daeb','#b2abd2')) +
  coord_flip() + scale_y_continuous(expand = c(0.02, 0)) +
  theme(panel.spacing = unit(0.8, "lines")) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "white")) +
  theme(legend.position = "none")



gp <- ggplotGrob(tight)

facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]

y.var <- sapply(ggplot_build(tight)$layout$panel_scales_y,
                function(l) length(l$range$range))
gp$widths[facet.columns] <- gp$widths[facet.columns] * y.var

grid::grid.draw(gp)

ggsave('Fig_2G.pdf', width = 16.5, height = 1.94, dpi = 600) #use h=10 for full, h=1.94 for Bif





######reduced fig 3 panel
shares <- read.csv("Fig3A_input.csv")
shares <- within(shares, 
                 Species <- factor(Species, 
                                   levels=names(sort(table(Species), 
                                                     increasing=TRUE))))

shares$Method = factor(shares$Method, levels=c('Transferred strains'))

shares$Species = factor(shares$Species, levels=c("Akkermansia muciniphila",
                                                 "Alistipes onderdonkii",
                                                 "Bacteroides caccae",
                                                 "Bacteroides cellulosilyticus",
                                                 "Bacteroides coprocola",
                                                 "Bacteroides coprophilus",
                                                 "Bacteroides dorei",
                                                 "Bacteroides eggerthii",
                                                 "Bacteroides fragilis",
                                                 "Bacteroides intestinalis",
                                                 "Bacteroides massiliensis",
                                                 "Bacteroides plebeius",
                                                 "Bacteroides stercoris",
                                                 "Bacteroides uniformis",
                                                 "Bacteroides vulgatus",
                                                 "Barnesiella intestinihominis",
                                                 "Bifidobacterium adolescentis",
                                                 "Bifidobacterium bifidum",
                                                 "Bifidobacterium breve",
                                                 "Bifidobacterium catenulatum",
                                                 "Bifidobacterium infantis",
                                                 "Bifidobacterium longum",
                                                 "Bifidobacterium pseudocatenulatum",
                                                 "Clostridium clostridioforme",
                                                 "Collinsella aerofaciens",
                                                 "Enterobacter asburiae",
                                                 "Escherichia coli",
                                                 "GGB32884 SGB58188",
                                                 "Lactobacillus gasseri",
                                                 "Lactobacillus paragasseri",
                                                 "Megasphaera massiliensis",
                                                 "Megasphaera sp DISK 18",
                                                 "Parabacteroides distasonis",
                                                 "Phascolarctobacterium faecium",
                                                 "Ruminococcus gnavus",
                                                 "Ruthenibacterium lactatiformans",
                                                 "Streptococcus salivarius",
                                                 "Sutterella wadsworthensis",
                                                 "Veillonella parvula"))




tight<- ggplot(shares, aes(x = Species, fill = Method, line = "black", border = "black")) +
  geom_bar(colour="black") +
  facet_grid(. ~ Method, scales = "fixed") +
  theme(axis.text.x = element_text(face = "italic", angle = 90, hjust = 0.95),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.3)) +
  scale_fill_manual(values=c('#b2abd2')) +
  scale_y_continuous(expand = c(0, 0.1)) +
  theme(panel.spacing = unit(0.8, "lines")) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "white")) +
  theme(legend.position = "none")



gp <- ggplotGrob(tight)

facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]

y.var <- sapply(ggplot_build(tight)$layout$panel_scales_y,
                function(l) length(l$range$range))
gp$widths[facet.columns] <- gp$widths[facet.columns] * y.var

grid::grid.draw(gp)


ggsave('Fig_3G.pdf', width = 11, height = 7.5, dpi = 600) #use h=10 for full, h=1.94 for Bif
