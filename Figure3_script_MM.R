{
  library(tidyverse)
  library(ggplot2)
  library(lsr)
  library(rcompanion)
  library(MASS)
  library(corrr)
  library("corrplot")
  library(RColorBrewer)
  library(RVAideMemoire)
  library(ggpubr)
  library(rstatix)
  library(gridExtra)
  library(latticeExtra)
  library(pdp)
  library(patchwork)
  library(svglite)
}



##########################################################################
#¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬#
# Figure 3 (b)                                                           #
#________________________________________________________________________#
##########################################################################
shares<- read.csv("csv_where_each_row_is_a_shared_species_event_and_each_column_is_metadata_associated_with_sample", header=T)


get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


mode<- shares %>%
  group_by(Species, Delivery_mode) %>%
  summarise(Number = length(Species)) %>%
  filter(!any(is.na(Delivery_mode))) 

A1 <- ggplot(mode, aes(x=Species, y=Number, fill=Delivery_mode)) +
  geom_bar(position="dodge", colour="black", stat = "identity") + coord_flip() +
  scale_fill_manual(values = c('#b35806','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788')) +
  theme (axis.text.y = element_text(size = 15, colour = "black", vjust = 0.5, hjust = 1, face = "italic")) +
  ylab("Total") + guides(fill=guide_legend("Delivery mode"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))+
  scale_y_continuous(expand = c(0, 0))
A1

A1_legend <- get_legend(A1)

ggsave('Delivery_mode_strains_shared.png', height = 10, width = 12, dpi = 300)


sex<- shares %>%
  group_by(Species, Child_sex) %>%
  summarise(Number = length(Species)) %>%
  filter(!any(is.na(Child_sex))) 

A2 <- ggplot(sex, aes(x=Species, y=Number, fill=Child_sex)) +
  geom_bar(position="dodge", colour="black", stat = "identity") + coord_flip() +
  scale_fill_manual(values = c('#b35806','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788')) +
  theme (axis.text.y = element_blank()) +
  ylab("Total") + guides(fill=guide_legend("Infant sex")) +
  scale_y_continuous(expand = c(0, 0))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
A2

A2_legend <- get_legend(A2)

ggsave('Child_sex_strains_shared.png', height = 10, width = 12, dpi = 300)

rupture<- shares %>%
  group_by(Species, Rupture_of_membranes) %>%
  summarise(Number = length(Species)) %>%
  filter(!any(is.na(Rupture_of_membranes))) 

A3 <- ggplot(rupture, aes(x=Species, y=Number, fill=Rupture_of_membranes)) +
  geom_bar(position="dodge", colour="black", stat = "identity") + coord_flip() +
  scale_fill_manual(values = c('#b35806','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788')) +
  theme (axis.text.y = element_blank()) +
  ylab("Total") + guides(fill=guide_legend("Rupture of membranes"))+
  scale_y_continuous(expand = c(0, 0))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
A3

A3_legend <- get_legend(A3)

ggsave('ROM_strains_shared.png', height = 10, width = 12, dpi = 300)

abxlab<- shares %>%
  group_by(Species, ABX_in_labour) %>%
  summarise(Number = length(Species)) %>%
  filter(!any(is.na(ABX_in_labour))) 

A4 <- ggplot(abxlab, aes(x=Species, y=Number, fill=ABX_in_labour)) +
  geom_bar(position="dodge", colour="black", stat = "identity") + coord_flip() +
  scale_fill_manual(values = c('#b35806','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788')) +
  theme (axis.text.y = element_blank()) +
  ylab("Total") + guides(fill=guide_legend("ABX in labour"))+
  scale_y_continuous(expand = c(0, 0))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
A4

A4_legend <- get_legend(A4)

ggsave('ABX_strains_shared.png', height = 10, width = 12, dpi = 300)


laborsetlab<- shares %>%
  group_by(Species, Labour_onset) %>%
  summarise(Number = length(Species)) %>%
  filter(!any(is.na(Labour_onset))) 

A5 <- ggplot(laborsetlab, aes(x=Species, y=Number, fill=Labour_onset)) +
  geom_bar(position="dodge", colour="black", stat = "identity") + coord_flip() +
  scale_fill_manual(values = c('#b35806','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788')) +
  theme (axis.text.y = element_blank()) +
  ylab("Total") + guides(fill=guide_legend("Labour onset"))+
  scale_y_continuous(expand = c(0, 0))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
A5

A5_legend <- get_legend(A5)

ggsave('ABX_strains_shared.png', height = 10, width = 12, dpi = 300)



Species_shared <- ggarrange(A1,A2+rremove("ylab"),A3+ rremove("ylab"),A4+ rremove("ylab"), A5+rremove("ylab"),
                            ncol = 5, nrow = 1, common.legend = FALSE, legend = "none",
                            widths = c(2.2,1,1,1,1))
Species_shared




ggsave('Fig3_B.png', height = 10, width = 20, dpi = 300)


ALegends <- ggarrange(A1_legend,A2_legend,A3_legend,A4_legend,A5_legend, ncol=5)
ALegends
ggsave('legend_for_Fig3_B.png', h = 5, w = 20)

##########################################################################
#¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬#
#Figure 3 (c)                                                            #
#________________________________________________________________________#
##########################################################################

df <- read.csv('csv_file_where_each_row_represents_each_dyad_and_each_column_has_corresponding_metadata_including_whether_sharing_event_occured_or_not')


##
#did individual chi-square-tests on variables (only sig shown here)
chisq_test(df$Sharing_event_Y_N,df$ABX_in_labour)

chisq_test(df$Sharing_event_Y_N,df$Child_sex)

chisq_test(df$Sharing_event_Y_N,df$Delivery_mode)

chisq_test(df$Sharing_event_Y_N,df$Labour_onset)

chisq_test(df$Sharing_event_Y_N,df$Rupture_of_membranes)
##


df$Sharing_event_Y_N <- factor(df$Sharing_event_Y_N , levels=c("Yes", "No"))
df$ABX_in_labour <- factor(df$ABX_in_labour , levels=c("Yes", "No"))

stat.test <-chisq_test(df$Sharing_event_Y_N,df$Delivery_mode)

Delmode <- table(df$Sharing_event_Y_N, df$Delivery_mode) 
Delmode <-as.data.frame((Delmode))
names(Delmode) = c('Sharing_event_Y_N', 'Delivery_mode', 'Freq')

B1<-ggplot(Delmode, aes(fill=Sharing_event_Y_N, y=Freq, x=Delivery_mode)) + 
  geom_bar(position="dodge", stat="identity", colour = "black")+
  xlab("Delivery mode") +
  scale_fill_manual(values= c('#b35806','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788'))+
  scale_y_continuous(name ="Dyads",breaks=seq(0,60,10), expand = c(0, 0))+
  guides(fill=guide_legend("Sharing event")) + labs(subtitle = get_test_label(stat.test, detailed = TRUE))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
B1 

ggsave("deliverymode_sharing.png")



stat.test <-chisq_test(df$Sharing_event_Y_N,df$Rupture_of_membranes)
ROM <- table(df$Sharing_event_Y_N, df$Rupture_of_membranes) 
ROM <-as.data.frame((ROM))
names(ROM) = c('Sharing_event_Y_N', 'Rupture_of_membranes', 'Freq')

B2<-ggplot(ROM, aes(fill=Sharing_event_Y_N, y=Freq, x=Rupture_of_membranes)) + 
  geom_bar(position="dodge", stat="identity", colour = "black")+
  xlab("Rupture of membranes") +
  scale_fill_manual(values= c('#b35806','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788'))+
  scale_y_continuous(name ="Dyads",breaks=seq(0,60,10), expand = c(0, 0))+
  guides(fill=guide_legend("Sharing event")) + labs(subtitle = get_test_label(stat.test, detailed = TRUE))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
B2

ggsave("ROM_sharing.png")


stat.test <-chisq_test(df$Sharing_event_Y_N,df$Child_sex)
sex <- table(df$Sharing_event_Y_N, df$Child_sex) 
sex <-as.data.frame((sex))
names(sex) = c('Sharing_event_Y_N', 'Child_sex', 'Freq')

B3<-ggplot(sex, aes(fill=Sharing_event_Y_N, y=Freq, x=Child_sex)) + 
  geom_bar(position="dodge", stat="identity", colour = "black")+
  xlab("Infant sex") +
  scale_fill_manual(values= c('#b35806','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788'))+
  scale_y_continuous(name ="Dyads",breaks=seq(0,60,10), expand = c(0, 0))+
  guides(fill=guide_legend("Sharing event")) + labs(subtitle = get_test_label(stat.test, detailed = TRUE))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
B3

ggsave("sex_sharing.png")


stat.test <-chisq_test(df$Sharing_event_Y_N,df$ABX_in_labour)
abx <- table(df$Sharing_event_Y_N, df$ABX_in_labour) 
abx <-as.data.frame((abx))
names(abx) = c('Sharing_event_Y_N', 'ABX_in_labour', 'Freq')

B4<-ggplot(abx, aes(fill=Sharing_event_Y_N, y=Freq, x=ABX_in_labour)) + 
  geom_bar(position="dodge", stat="identity", colour = "black")+
  xlab("Antibiotics during labour") +
  scale_fill_manual(values= c('#b35806','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788'))+
  scale_y_continuous(name ="Dyads",breaks=seq(0,60,10), expand = c(0, 0))+
  guides(fill=guide_legend("Sharing event")) + labs(subtitle = get_test_label(stat.test, detailed = TRUE))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
B4
ggsave("ABX_sharing.png")



stat.test <-chisq_test(df$Sharing_event_Y_N,df$Labour_onset)
lab <- table(df$Sharing_event_Y_N, df$Labour_onset) 
lab <-as.data.frame((lab))
names(lab) = c('Sharing_event_Y_N', 'Labour_onset', 'Freq')

B5<-ggplot(lab, aes(fill=Sharing_event_Y_N, y=Freq, x=Labour_onset)) + 
  geom_bar(position="dodge", stat="identity", colour = "black")+
  xlab("Labour onset") +
  scale_fill_manual(values= c('#b35806','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788'))+
  scale_y_continuous(name ="Dyads",breaks=seq(0,60,10), expand = c(0, 0))+
  guides(fill=guide_legend("Sharing event")) + labs(subtitle = get_test_label(stat.test, detailed = TRUE))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
B5
ggsave("Labout_onset_sharing.png")

Sig_sharing_factors <- ggarrange(B1,B2,B3,B4,B5, ncol = 2, nrow = 3, common.legend = TRUE, legend = "bottom")

Sig_sharing_factors
ggsave("Fig3_C.png", w = 10, h = 10)

plot_list<- list(Species_shared,Sig_sharing_factors)


Figure3BC <- ggarrange(Species_shared,Sig_sharing_factors, ALegends,align = "v",
                            ncol = 2, nrow = 2, common.legend = TRUE, legend = "right",
                     heights = c(1.5,0.4,2.0), labels = c("(a)","(b)"))

Figure3BC <- cowplot::ggdraw(Figure3BC) + 
  theme(plot.background = element_rect(fill="white", color = NA))

Figure3BC

ggsave("Fig3_BC.png", w = 22, h = 12)


#Supplementary Figure 3 (f)
alldata <- read.csv('Fig_S3f_input.csv', header = T)



###>>>Infant 1 month<<<###

df <- alldata 


empty_rows <- df[rowSums(df[,-c(1:101)]) == 0, ]
df <- df[rowSums(df[,-c(1:101)]) > 0, ]


meta <- df[,c(7:42,47)]
species <- df[,-c(1:101)]

infant_1M.mds <- metaMDS(species, distance = "bray", autotransform = FALSE)

envfit_meta <- envfit(infant_1M.mds, cbind(df[,c(7:42,47)]), na.rm = T)
envfit_meta 

str(envfit_meta)

envfit_meta_dffactors<- data.frame((envfit_meta$factors)$r, (envfit_meta$factors)$pvals )
(names(envfit_meta_dffactors)[1] <- "R2")
(names(envfit_meta_dffactors)[2] <- "pval")


envfit_meta_factors <- as.data.frame(scores(envfit_meta, display = "factors"))
envfit_meta_factors$Variable <- gsub('Sharing_event_Y_N', 'Sharing_event_Y_N',rownames(envfit_meta_factors))  

data.scores <- as.data.frame(scores(infant_1M.mds))


Cred<-subset(envfit_meta_dffactors,pval<0.05)

tops <- envfit_meta_factors[-c(1:34,35,37,39,40,41,42,43:48,49,51:76),]




p <- ggplot(data.scores) +
  geom_point(mapping = aes(x = NMDS1, y = NMDS2),  size = 4, colour = "#808080") +
  
  coord_fixed() + ## need aspect ratio of 1!
  
  geom_segment(data = tops,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "red") +
  geom_text(data = tops, aes(x = NMDS1, y = NMDS2, label = Variable),
            size = 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black",
                                        size = 1, linetype = "solid"))
p  

ggsave('FigS3_f.png', p, height = 8, width = 10, dpi = 300)  


#Supplementary Figure 3 (g)

data<-read.csv("abundance_of_transfered_strains.csv")
data_long <- gather(data, Species, Relative_abundance, 'Veillonella_parvula':'Akkermansia_muciniphila', factor_key=TRUE)
write.csv(data_long,"abund_transfer_long.csv")
#manually add a genus column to left of specie column, remove underscore from species names
data_long<-read.csv("abund_transfer_long.csv")
View(data_long)



stat.test_alpha <- data_long %>%
  group_by(Species) %>%
  t_test(Relative_abundance ~ Sharing_event_Y_N) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
View(stat.test_alpha)




p <- ggplot(data_long, aes(Species, Relative_abundance))
p + geom_boxplot(notch = TRUE)



p<-ggplot(data_long, aes(x=Species, y=Relative_abundance, fill=Sharing_event_Y_N)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Species, scales = "free") +
  theme(axis.text.x = element_text(angle = 0, vjust = 01, hjust=1)) +
  theme(axis.text.y = element_text(face = "italic")) +
  scale_fill_manual(values = c("#7F7F7F", "#5AB4AC")) +
  theme(panel.background = element_rect(fill = "white", colour = "black",
                                        size = 1, linetype = "solid")) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  coord_flip()


p


p_tag <- c("", "","","","","","","","","","","","","*","","","",
           "","","","","","","","","","","","","","","")

sig_plots<-tag_facet(p,
                     x = 0.9, y = 20,
                     vjust = 0, hjust = 0,
                     open = "", close = "",
                     fontface = 4,
                     size = 5,
                     color = "black",
                     family = "serif",
                     tag_pool = p_tag)


sig_plots<-sig_plots + theme(strip.text = element_text(face="italic", size=8),
                             strip.background = element_rect(size=2))
sig_plots


ggsave('FigS3_g.png', height = 8, width = 12, dpi = 300)

#Supplementary Figure 3 (h)
###all 1M infant species above 1% realtvie abundance in at least 1 sample
df <- read.csv("Transfer_effect_heatmaps.csv", header = T)
dim(df)
matrix <- as.matrix(df[,2:119], "numeric")
rownames(matrix)<-df[,1]


col_fun <- colorRamp2(c(0, 10, 40, 50,60,70, 80,90, 100), c('#F2F2F2','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506','#54278f','#3f007d'))

pheatmap_meta_all <-read.csv("Transfer_effect_meta.csv")
dim(pheatmap_meta_all)


pheatmap_meta_all$Bifidobacterium_breve_transfer <- factor(pheatmap_meta_all$Bifidobacterium_breve_transfer , levels=c("Yes", "No"))
pheatmap_meta_all$Bifidobacterium_bifidum_transfer <- factor(pheatmap_meta_all$Bifidobacterium_bifidum_transfer , levels=c("Yes", "No"))
pheatmap_meta_all$Escherichia_coli_transfer <- factor(pheatmap_meta_all$Escherichia_coli_transfer , levels=c("Yes", "No"))






ann <- data.frame(pheatmap_meta_all$Bifidobacterium_breve_transfer,pheatmap_meta_all$Bifidobacterium_bifidum_transfer,pheatmap_meta_all$Escherichia_coli_transfer)
colnames(ann) <- c("Bifidobacterium_breve_transfer", "Bifidobacterium_bifidum_transfer","Escherichia_coli_transfer")
colours <- list("Bifidobacterium_breve_transfer"=c("Yes"="#ef8a62",
                                                   "No"="#252525"),
                "Bifidobacterium_bifidum_transfer"=c("Yes"="#ef8a62",
                                                     "No"="#252525"),
                "Escherichia_coli_transfer"=c("Yes"="#ef8a62",
                                              "No"="#252525"))


colAnn <- HeatmapAnnotation(df=ann, which="col", col=colours,
                            annotation_width=unit(c(1, 2), "cm"),
                            gap=unit(1, "mm"),
                            simple_anno_size = unit(0.3, "cm"))

pdf("Test_heatmap4.pdf", width = 20, height = 15)
png("Test_heatmap4.png",width=20,height=15,units="in",res=600)

ht <- Heatmap(matrix, border = T, name = "Relative abund", col = col_fun,
              heatmap_legend_param = list(
                at = c(0, 25, 50, 75, 100),
                title = "Relative abundance",
                legend_height = unit(5, "cm"),
                direction = "horizontal",
                title_position = "topcenter"),
              column_title = "title", column_title_side = "bottom",
              column_title_gp = gpar(col = "white", lwd = 0.1),
              column_dend_height = unit(2, "cm"),
              row_title = "Species",
              cluster_rows = TRUE, show_row_dend = TRUE,
              clustering_distance_columns = "euclidean", column_dend_reorder = TRUE,
              show_column_names = FALSE,
              row_names_gp = gpar(fontsize = 8, fontface = "italic"),
              column_gap = unit(2, "mm"),
              cluster_column_slices = FALSE,
              heatmap_height = unit(30, "cm"),
              heatmap_width = unit(25, "cm"),
              top_annotation=colAnn, use_raster = TRUE, raster_quality = 5)

draw(ht, heatmap_legend_side = "bottom")
dev.off()
}
  


