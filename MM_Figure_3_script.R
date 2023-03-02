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

#######
#Figure 3 a
#######
{
data3<-read.csv("Figure_3A_input_MM.csv") 

data3a<- data3 %>% filter(Sample.Strain==Sample.Matched,Sample.time.Strain==Sample.time.Matched,Strain!=Matched,Species.Strain==Species.Matched)

ggplot(data3a, aes(`ANI.`,fill=`Sample.Type.Strain`))+
  geom_histogram(binwidth = 0.1)+
  scale_x_continuous(name="ANI (%)",breaks = seq(97.7,100,by=0.1))+
  scale_y_continuous(limits = c(0,400), expand = c(0, 0)) +
  theme_classic(base_size = 12)+scale_fill_manual(breaks = c("Breast Milk","Infant","Maternal"),values = c("#2F5597","#AFABAB","#70AD47"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ylab("No. of pairwise comparisons")+
  labs(fill="Sample type")+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size=10)) +
  theme(panel.border = element_rect(colour="black", fill=NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(title = "Sample Type", title.position = "top", title.hjust = 0.5))

ggsave("Figure_3A_MM.png",w=4,h=4,dpi=600)
}

#######
#Figure 3 b
#######
{
  data3b<- data3 %>% filter(Sample.Strain==Sample.Matched,Strain!=Matched,Species.Strain==Species.Matched,`ANI.`>=99.9, str_detect(Strain,"^M"),str_detect(Matched,"^M"),str_detect(Sample.Strain,"^P")) 
  #write.csv(data3b,"data3b.csv")
  
  ggplot(data3b,aes(`snpdist`,`Sample.Type.Strain`,colour=`Sample.Type.Strain`))+
    geom_boxplot()+
    geom_jitter()+
    scale_x_continuous(name="SNP distance",n.breaks=10)+theme_classic(base_size = 12)+
    scale_color_manual(breaks = c("Breast Milk","Infant","Maternal"),values = c("#2F5597","#AFABAB","#70AD47"))+
    ylab("Sample")+
    theme(legend.position = "bottom")+
    labs(colour="Sample type")+
    ylab("Individual")+coord_flip()+
    theme(panel.border = element_rect(colour="black", fill=NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(title = "Sample Type", title.position = "top", title.hjust = 0.5))
  
  ggsave("Figure_3B_MM.png",w=5,h=4,dpi=600)
}




