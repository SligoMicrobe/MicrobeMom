library(ggstatsplot)
library("ggpval")


#######
#Figure 2 a
#######
{
relabun<-read.csv("Figure_2A_input_MM.csv")

relabun$Sample_Type <- factor(relabun$Sample_Type , levels=c("Breast Milk", "Infant stool", "Maternal stool"))


vbx<-ggplot(relabun, aes(x=Sample_Type, y=RA_species_source_microbiome, fill=seq_type, colour = seq_type)) + 
  geom_violin(position=position_dodge(1),lwd=1)+
  geom_boxplot(width=0.1,position=position_dodge(1),lwd = 1.5) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values=c("white", "white")) +
  scale_color_manual(values=(c("#d9d9d9", "#525252"))) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    linewidth = 1, linetype = "solid"),
    panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                    colour = "white")) +
  ylab("Relative abundance of species in source microbiome (%)") +
  xlab("") +
  labs(colour = "Genome type") +
  theme(legend.position="bottom")+ guides(fill = "none") +coord_flip() 


vbx
stat.test<-compare_means(RA_species_source_microbiome~seq_type, data = subset(relabun,ave(RA_species_source_microbiome,Sample_Type,seq_type,FUN=length)>1), method="t.test", group.by = "Sample_Type")

stat.test

vbx + stat_pvalue_manual(stat.test,label = "p.format", y.position=2)


ggsave('Figure_2A_MM.png', height = 6, width = 6, dpi = 600)
}


#######
#Figure 2 b
#######
{
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
  library(hrbrthemes)
  library(ggplot2)
  library(Quandl)
  Quandl.api_key("XXX")
  library(ggpmisc)
}

MAGScheck<-read.csv("Figure_2B_input_MM.csv")
dim(MAGScheck)
MAGScheck<- MAGScheck %>% filter(Reads=='Good')
Mstool<- MAGScheck %>% filter(Sample_type=='Maternal_stool')
Istool<- MAGScheck %>% filter(Sample_type=='Infant_stool')
Vaginal<- MAGScheck %>% filter(Sample_type=='Vaginal')
Oral<- MAGScheck %>% filter(Sample_type=='Maternal_oral')
Breast<- MAGScheck %>% filter(Sample_type=='Breast_milk')


lmMstool <- lm(log_reads~MAG_count, data = Mstool)
lmIstool <- lm(log_reads~MAG_count, data = Istool)
lmVaginal <- lm(log_reads~MAG_count, data = Vaginal)
lmOral <- lm(log_reads~MAG_count, data = Oral)
lmBreast <- lm(log_reads~MAG_count, data = Breast)




lm_eqn <- function(df){
  m <- lm(MAG_count ~ log_reads, Mstool);
  eq <- substitute(italic(y) == a + b~~italic(.x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


PB<- ggplot(Mstool, aes(x=log_reads, y=MAG_count)) + geom_point(shape=1) + geom_smooth(method=lm, se=TRUE) +
  facet_wrap(~Sample_type) +
  ggtitle("") +
  scale_x_continuous(name = "High quality metagenomic reads (log10)") +
  scale_y_continuous(name = "Number of MAGS") +
  geom_vline(xintercept=c(6.45591), linetype="dotted", colour="blue") +
  annotate("rect", xmin = 5.5, xmax = 7.5, ymin = 23.5, ymax = 26.5, fill="white", colour="red") +
  geom_text(x = 6.5, y = 25, label = lm_eqn(df), parse = TRUE) +
  theme_bw() 

PB 

lm_eqn <- function(df){
  m <- lm(MAG_count ~ log_reads, Istool);
  eq <- substitute(italic(y) == a + b~~italic(.x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}
PC<- ggplot(Istool, aes(x=log_reads, y=MAG_count)) + geom_point(shape=1) + geom_smooth(method=lm, se=TRUE) +
  facet_wrap(~Sample_type) +
  ggtitle("") +
  scale_x_continuous(name = "High quality metagenomic reads (log10)") +
  scale_y_continuous(name = "Number of MAGS") +
  geom_vline(xintercept=c(6.0716), linetype="dotted", colour="blue") +
  annotate("rect", xmin = 5, xmax = 7.45, ymin = 9.2, ymax = 10.8, fill="white", colour="red") +
  geom_text(x = 6.2, y = 10, label = lm_eqn(df), parse = TRUE) +
  theme_bw()
PC

lm_eqn <- function(df){
  m <- lm(MAG_count ~ log_reads, Vaginal);
  eq <- substitute(italic(y) == a + b~~italic(.x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}
PV<- ggplot(Vaginal, aes(x=log_reads, y=MAG_count)) + geom_point(shape=1) + geom_smooth(method=lm, se=TRUE) +
  facet_wrap(~Sample_type) +
  ggtitle("") +
  scale_x_continuous(name = "High quality metagenomic reads (log10)") +
  scale_y_continuous(name = "Number of MAGS") +
  geom_vline(xintercept=c(5.64540), linetype="dotted", colour="blue") +
  annotate("rect", xmin = 5, xmax = 6.95, ymin = 4.5, ymax = 5.5, fill="white", colour="red") +
  geom_text(x = 6, y = 5, label = lm_eqn(df), parse = TRUE) +
  theme_bw()
PV


lm_eqn <- function(df){
  m <- lm(MAG_count ~ log_reads, Oral);
  eq <- substitute(italic(y) == a + b~~italic(.x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}
PO<- ggplot(Oral, aes(x=log_reads, y=MAG_count)) + geom_point(shape=1) + geom_smooth(method=lm, se=TRUE) +
  facet_wrap(~Sample_type) +
  ggtitle("") +
  scale_x_continuous(name = "High quality metagenomic reads (log10)") +
  scale_y_continuous(name = "Number of MAGS") +
  geom_vline(xintercept=c(5.7921), linetype="dotted", colour="blue") +
  annotate("rect", xmin = 5.15, xmax = 6.6, ymin = 3.3, ymax = 3.9, fill="white", colour="red") +
  geom_text(x = 5.9, y = 3.6, label = lm_eqn(df), parse = TRUE) +
  theme_bw()
PO

lm_eqn <- function(df){
  m <- lm(MAG_count ~ log_reads, Breast);
  eq <- substitute(italic(y) == a + b~~italic(.x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}
PM<- ggplot(Breast, aes(x=log_reads, y=MAG_count)) + geom_point(shape=1) + geom_smooth(method=lm, se=TRUE) +
  facet_wrap(~Sample_type) +
  ggtitle("") +
  scale_x_continuous(name = "High quality metagenomic reads (log10)") +
  scale_y_continuous(name = "Number of MAGS") +
  geom_vline(xintercept=c(5.8783), linetype="dotted", colour="blue") +
  annotate("rect", xmin = 5, xmax = 7.9, ymin = 4, ymax = 5, fill="white", colour="red") +
  geom_text(x = 6.5, y = 4.5, label = lm_eqn(df), parse = TRUE) +
  theme_bw() 
PM
summary(lmMstool) 
summary(lmIstool) 
summary(lmVaginal) 
summary(lmOral) 
summary(lmBreast) 


MAGSvReads<- ggarrange(PB, PC, PO, PV,PM, nrow = 1, ncol = 5,
                       font.label = list(size = 20),
                       common.legend = TRUE, legend = "bottom") 




ggsave("Figure_2B_MM.png", w=17, h=5, dpi=600)
}


#######
#Figure 2 c
#######
{
library("ggvenn")


#Make the plot
vennMI<-read.csv("Figure_2C_input_MM.csv",header=T)
head(vennMI)
VMI<-apply(vennMI,2,as.list)


ggvenn(VMI, columns = c("Isolate", "MAG"),fill_color = c("#d9d9d9", "#525252"),
       stroke_size = 0.5, auto_scale = TRUE)

ggsave("Figure_2C_MM.png", w=10, h=10, dpi=600)

}

#######
#Supplementary Figure 2 b
#######
{

dataS2B<-read.csv("Supplementary_Figure_2B_input_MM.csv")

dataS2B$isolation.type <- factor(dataS2B$isolation_type , levels=c("Untargeted isolation", "Targeted isolation"))

ggplot(dataS2B,aes(NormDistEL,iso_type,))+geom_boxplot()+geom_point()+theme_classic2()+ylab("Identification method")+xlab("Normalised phylogenetic distance") +
  theme(panel.border = element_rect(colour="black", fill=NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  scale_x_continuous(trans='log10') +
  geom_vline(xintercept=c(1E-02), linetype="dotted", colour="blue")

ggsave("Supplementary_Figure_2B_MM.png",w=6,h=6,dpi=600)

}



#######
#Supplementary Figure 2 c
#######
#Figure S2c:
{
dataS2C<-read.csv("Supplementary_Figure_2C_input_MM.csv")


plot<-ggplot(dataS2C,aes(`snpdist`,`comparison.type`,colour=`comparison.type`))+
  geom_boxplot()+
  geom_jitter()+
  scale_x_continuous(trans="log10",name="Log10(SNP distance)",n.breaks=10)+theme_classic(base_size = 12)+
  scale_color_manual(breaks = c("Within dyad","Between dyad"),values = c("#2F5597","#AFABAB"))+
  ylab("Log10(SNP distance)")+
  labs(colour="")+
  ylab("")+coord_flip()+
  theme(panel.border = element_rect(colour="black", fill=NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())


stat.test2<-compare_means(snpdist~comparison.type, data = dataS2C, method="t.test")
stat.test2

plot + annotate("", x = 1, y = 30005, label = "p = 1.1e-12")


ggsave("Supplementary_Figure_2C_MM.png",w=7,h=6,dpi=600)

}




