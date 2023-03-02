#Libraries
{
library(tidyverse)
library("ggpubr")
library(ggplot2)
library(forcats)
library(reshape2)
library(gghighlight)
library(qdapTools)
library(ape)
library(readxl)
library(writexl)
library(RColorBrewer)
library(scales)
library(pheatmap)
library(patchwork)
library(genoPlotR)
#hrbrthemes::import_roboto_condensed()
library(hrbrthemes)
library(Matrix)
library(UpSetR)
library(ggplotify)
}


#######
#Fig 4 a
#######
{
data4a<-read.csv("Figure_4A_input_MM.csv")  

data4a$Sample.source <- factor(data4a$Sample.source , levels=c("Maternal", "Infant", "Breast Milk"))


ggplot(data4a, aes(`mapped.reads.coverage.of.reference`,`Dyad_Spec`))+
  geom_point(alpha=0.5,aes(size=Average_coverage,color=`Sample.source`,fill=`Sample.source`),na.rm = T,shape=21)+
  theme_classic(base_size = 12)+
  scale_fill_manual(breaks = c("Maternal","Infant","Breast Milk"),values = c("#70AD47","#AFABAB","#2F5597"))+
  scale_color_manual(breaks = c("Maternal","Infant","Breast Milk"),values = c("#70AD47","#AFABAB","#2F5597"))+
  scale_size_continuous(breaks = c(0,2,20,200,400),name="Average read \nmapped per base")+
  geom_path(aes(group=Dyad_Spec), color="gray",linetype="dotted")+
  guides(colour = guide_legend(override.aes = list(shape = 15,size=5, alpha=1)))+
  #labs(title="Metagenome reads mapped to isolated genomes - Transmitted strains")+
  xlab("Genome Coverage (%)")+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1),axis.title.x = element_blank()) +
  theme(panel.border = element_rect(colour="black", fill=NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))

ggsave("Figure_4A_MM.png",w=9,h=5.5,dpi=600)
}

#######
#Fig 4 b
#######
{
data4b<-read.csv("Figure_4B_input_MM.csv")

data4b$isolation.type <- factor(data4b$isolation.type , levels=c("Untargeted isolation", "Targeted isolation"))

ggplot(data4b,aes(Final.snpdist,isolation.type,))+geom_boxplot()+geom_point()+coord_flip()+theme_classic2()+ylab("")+xlab("SNP distance") +
  theme(panel.border = element_rect(colour="black", fill=NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave("Figure_4B_MM.png",w=3,h=6,dpi=600)
}


#######
#Fig 4 c
#######
{
transmissions.upset.data<-c(
  Isolations=16,
  inStrain=1,
  StrainPhlAn=10,
  MixedMethods=5,
  "Isolations&StrainPhlAn&inStrain"=4,
  "inStrain&StrainPhlAn"=6,
  "Isolations&StrainPhlAn"=7,
  "StrainPhlAn&MixedMethods"=4,
  "Isolations&MixedMethods"=0,
  "inStrain&StrainPhlAn&MixedMethods"=0,
  "Isolations&StrainPhlAn&MixedMethods"=0,
  "Isolations&inStrain&MixedMethods"=0,
  "Isolations&StrainPhlAn&inStrain&MixedMethods"=0
)




upset(fromExpression(transmissions.upset.data),
      nintersects = 10, 
      nsets = 4, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8,
      line.size = 1,
      sets.x.label = "Total transmissions detected by method",
      mainbar.y.label = "Transmissions detected",
      main.bar.color = c('#615e5e', '#615e5e', '#615e5e','#615e5e','#615e5e','#615e5e','#615e5e','#615e5e'),
      sets.bar.color =c('#6b6464','#6b6464','#6b6464', '#6b6464') )

#Exported figure as Figure_4C_MM.svg from plot preview w = 1000, h = 500

}

#######
#Fig 4d
#######
{
shares <- read.csv("Figure_4D_input_MM.csv")

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
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.3)) +
  scale_fill_manual(values=c('#636363','#636363','#636363','#636363','#636363','#636363')) +
  coord_flip() + scale_y_continuous(expand = c(0.02, 0)) +
  theme(panel.spacing = unit(0.8, "lines")) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    linewidth = 1, linetype = "solid"),
    panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                    colour = "white")) +
  theme(legend.position = "none")



gp <- ggplotGrob(tight)

facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]

y.var <- sapply(ggplot_build(tight)$layout$panel_scales_y,
                function(l) length(l$range$range))
gp$widths[facet.columns] <- gp$widths[facet.columns] * y.var

grid::grid.draw(gp)

ggsave('Figure_4D_MM.png', width = 16.5, height = 1.94, dpi = 600) 
}
