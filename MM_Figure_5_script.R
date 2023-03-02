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



#######
#Fig 5 a
#######
{
df <- read.csv("Figure_5A_input_MM.csv", header = T)
dim(df)
matrix <- as.matrix(df[,2:14], "numeric")
rownames(matrix)<-df[,1]


col_fun <- circlize::colorRamp2(c(0, 2, 5, 10,15,20), c('white','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6'))


##top annotation
heatmap_meta_all <-read.csv("Figure_5A_metadata_top_MM.csv")
dim(heatmap_meta_all)

heatmap_meta_all$Category <- factor(heatmap_meta_all$Category , levels=c("Delivery mode", "Infant sex", "Rupture of membranes", "ABX in labour", "Labour onset"))


ann <- data.frame(heatmap_meta_all$Category)

colnames(ann) <- c("Category")
colours <- list("Category"=c("Delivery mode"="#8dd3c7",
                             "Infant sex"="#ffffb3",
                             "Rupture of membranes"="#bebada",
                             "ABX in labour"="#fb8072",
                             "Labour onset"="#80b1d3"))


colAnn <- HeatmapAnnotation(df=ann, which="col", col=colours,
                            annotation_width=unit(c(8, 8), "cm"),
                            gap=unit(1, "mm"),
                            simple_anno_size = unit(0.8, "cm"))

##side annotation
sideann<-read.csv("Figure_5A_metadata_side_MM.csv")
rowann <- data.frame(sideann$Total_shared_count)
colnames(rowann) <- c("Total_shared_count")
row_ha = rowAnnotation(Total_shared_count=anno_barplot(rowann),width = unit(4, "cm") )

pdf("Figure_5A_MM.pdf", width = 10, height = 12)
png("Figure_5A_MM.png",width=10,height=18,units="in",res=600)

ht <- Heatmap(matrix, border = T, name = "Sharing occurences", col = col_fun,
              rect_gp = gpar(col = "black", lwd = 1),
              width = ncol(matrix)*unit(10, "mm"), 
              height = nrow(matrix)*unit(10, "mm"),
              heatmap_legend_param = list(
                at = c(0, 5, 10, 15, 20),
                title = "Sharing occurences",
                legend_height = unit(8, "cm"),
                direction = "horizontal",
                title_position = "topcenter"),
              column_title = "title", column_title_side = "bottom",
              column_title_gp = gpar(col = "white", lwd = 0.1),
              row_title = "Species",
              row_title_gp = gpar(fontsize = 20),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = TRUE,
              column_names_rot = 45,
              row_names_gp = gpar(fontsize = 12, fontface = "italic"),
              row_names_side = "left",
              column_split = factor(heatmap_meta_all$Category , levels=c("Delivery_mode", "Infant_sex", "Rupture_of_membranes", "ABX_in_labour", "Labour_onset"),
                                    labels=c("Delivery_mode", "Infant_sex", "Rupture_of_membranes", "ABX_in_labour", "Labour_onset")),
              column_gap = unit(2, "mm"),
              cluster_column_slices = FALSE,
              top_annotation=colAnn,
              use_raster = TRUE, raster_quality = 5,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(matrix[i, j] > 0.1)
                  grid.text(sprintf("%.1f", matrix[i, j]), x, y, gp = gpar(fontsize = 10, col = "black"))
              },
              right_annotation = row_ha)

ht<-draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")


dev.off()
}


#######
#Fig 5 b
#######
{
df <- read.csv('Figure_5B_input_MM.csv')


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
  scale_fill_manual(values= c('#252525','#bdbdbd'))+
  scale_y_continuous(name ="Dyads",breaks=seq(0,60,10), expand = c(0, 0))+
  guides(fill=guide_legend("Sharing event")) + labs(subtitle = get_test_label(stat.test, detailed = TRUE))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
B1 

#ggsave("deliverymode_sharing.png")



stat.test <-chisq_test(df$Sharing_event_Y_N,df$Rupture_of_membranes)
ROM <- table(df$Sharing_event_Y_N, df$Rupture_of_membranes) 
ROM <-as.data.frame((ROM))
names(ROM) = c('Sharing_event_Y_N', 'Rupture_of_membranes', 'Freq')

B2<-ggplot(ROM, aes(fill=Sharing_event_Y_N, y=Freq, x=Rupture_of_membranes)) + 
  geom_bar(position="dodge", stat="identity", colour = "black")+
  xlab("Rupture of membranes") +
  scale_fill_manual(values= c('#252525','#bdbdbd'))+
  scale_y_continuous(name ="Dyads",breaks=seq(0,60,10), expand = c(0, 0))+
  guides(fill=guide_legend("Sharing event")) + labs(subtitle = get_test_label(stat.test, detailed = TRUE))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
B2

#ggsave("ROM_sharing.png")


stat.test <-chisq_test(df$Sharing_event_Y_N,df$Child_sex)
sex <- table(df$Sharing_event_Y_N, df$Child_sex) 
sex <-as.data.frame((sex))
names(sex) = c('Sharing_event_Y_N', 'Child_sex', 'Freq')

B3<-ggplot(sex, aes(fill=Sharing_event_Y_N, y=Freq, x=Child_sex)) + 
  geom_bar(position="dodge", stat="identity", colour = "black")+
  xlab("Infant sex") +
  scale_fill_manual(values= c('#252525','#bdbdbd'))+
  scale_y_continuous(name ="Dyads",breaks=seq(0,60,10), expand = c(0, 0))+
  guides(fill=guide_legend("Sharing event")) + labs(subtitle = get_test_label(stat.test, detailed = TRUE))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
B3

#ggsave("sex_sharing.png")


stat.test <-chisq_test(df$Sharing_event_Y_N,df$ABX_in_labour)
abx <- table(df$Sharing_event_Y_N, df$ABX_in_labour) 
abx <-as.data.frame((abx))
names(abx) = c('Sharing_event_Y_N', 'ABX_in_labour', 'Freq')

B4<-ggplot(abx, aes(fill=Sharing_event_Y_N, y=Freq, x=ABX_in_labour)) + 
  geom_bar(position="dodge", stat="identity", colour = "black")+
  xlab("Antibiotics during labour") +
  scale_fill_manual(values= c('#252525','#bdbdbd'))+
  scale_y_continuous(name ="Dyads",breaks=seq(0,60,10), expand = c(0, 0))+
  guides(fill=guide_legend("Sharing event")) + labs(subtitle = get_test_label(stat.test, detailed = TRUE))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
B4
#ggsave("ABX_sharing.png")



stat.test <-chisq_test(df$Sharing_event_Y_N,df$Labour_onset)
lab <- table(df$Sharing_event_Y_N, df$Labour_onset) 
lab <-as.data.frame((lab))
names(lab) = c('Sharing_event_Y_N', 'Labour_onset', 'Freq')

B5<-ggplot(lab, aes(fill=Sharing_event_Y_N, y=Freq, x=Labour_onset)) + 
  geom_bar(position="dodge", stat="identity", colour = "black")+
  xlab("Labour onset") +
  scale_fill_manual(values= c('#252525','#bdbdbd'))+
  scale_y_continuous(name ="Dyads",breaks=seq(0,60,10), expand = c(0, 0))+
  guides(fill=guide_legend("Sharing event")) + labs(subtitle = get_test_label(stat.test, detailed = TRUE))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
B5
#ggsave("Labout_onset_sharing.png")

Sig_sharing_factors <- ggarrange(B1,B2,B3,B4,B5, ncol = 2, nrow = 3, common.legend = TRUE, legend = "bottom")

Sig_sharing_factors

ggsave("Figure_5B_MM.png", w = 10, h = 10)
}

#######
#Fig 5 c
#######
{
df <- read.csv('Figure_5C_input_MM.csv')


##
#did individual chi-square-tests on variables (only sig shown here)
chisq_test(df$Labour_onset,df$Bacteroides_dorei)

chisq_test(df$Labour_onset,df$Bifidobacterium_bifidum)

chisq_test(df$Labour_onset,df$Bifidobacterium_breve)

chisq_test(df$Labour_onset,df$Bifidobacterium_longum)

chisq_test(df$Labour_onset,df$Bacteroides_uniformis)

chisq_test(df$Labour_onset,df$Bacteroides_vulgatus)

chisq_test(df$Labour_onset,df$Bacteroides_coprophilus)

chisq_test(df$Labour_onset,df$Barnesiella_intestinihominis)

chisq_test(df$Labour_onset,df$Collinsella_aerofaciens)

chisq_test(df$Labour_onset,df$Megasphaera_massiliensis)

chisq_test(df$Labour_onset,df$Megasphaera_sp_DISK_18)

chisq_test(df$Labour_onset,df$Veillonella_parvula)

chisq_test(df$Labour_onset,df$Akkermansia_muciniphila)

chisq_test(df$Labour_onset,df$Bacteroides_eggerthii)

chisq_test(df$Labour_onset,df$Bifidobacterium_adolescentis)

chisq_test(df$Labour_onset,df$Sutterella_wadsworthensis)

chisq_test(df$Labour_onset,df$Ruminococcus_gnavus)

chisq_test(df$Labour_onset,df$Bacteroides_intestinalis)

chisq_test(df$Labour_onset,df$Bacteroides_stercoris)

chisq_test(df$Labour_onset,df$Bacteroides_massiliensis)

chisq_test(df$Labour_onset,df$Alistipes_onderdonkii)

chisq_test(df$Labour_onset,df$Bacteroides_fragilis)

chisq_test(df$Labour_onset,df$Bacteroides_caccae)

chisq_test(df$Labour_onset,df$Escherichia_coli)

chisq_test(df$Labour_onset,df$Phascolarctobacterium_faecium)

chisq_test(df$Labour_onset,df$Bacteroides_cellulosilyticus)

chisq_test(df$Labour_onset,df$GGB32884_SGB58188)

chisq_test(df$Labour_onset,df$Bacteroides_plebeius)

chisq_test(df$Labour_onset,df$Clostridium_clostridioforme)

chisq_test(df$Labour_onset,df$Ruthenibacterium_lactatiformans)

chisq_test(df$Labour_onset,df$Lactobacillus_gasseri)

chisq_test(df$Labour_onset,df$Lactobacillus_paragasseri)

chisq_test(df$Labour_onset,df$Streptococcus_salivarius)

chisq_test(df$Labour_onset,df$Ruthenibacterium_lactatiformans)

chisq_test(df$Labour_onset,df$Bacteroides_coprocola)

chisq_test(df$Labour_onset,df$Bifidobacterium_pseudocatenulatum)

chisq_test(df$Labour_onset,df$Bifidobacterium_catenulatum)
##

df$Delivery_mode <- factor(df$Delivery_mode , levels=c("Vaginal", "Caesarian"))
df$Bacteroides_vulgatus <- factor(df$Bacteroides_vulgatus , levels=c("Yes", "No"))
df$Bifidobacterium_bifidum <- factor(df$Bifidobacterium_bifidum , levels=c("Yes", "No"))
df$ABX_in_labour <- factor(df$ABX_in_labour , levels=c("Yes", "No"))


#delivery mode sig species transferred#
stat.test <-chisq_test(df$Delivery_mode,df$Bacteroides_vulgatus)

vulgatusDelmode <- table(df$Delivery_mode, df$Bacteroides_vulgatus) 
vulgatusDelmode <-as.data.frame((vulgatusDelmode))
names(vulgatusDelmode) = c('Delivery_mode', 'Bacteroides_vulgatus', 'Freq')

D1<-ggplot(vulgatusDelmode, aes(fill=Delivery_mode, y=Freq, x=Bacteroides_vulgatus), size = 40) + 
  geom_bar(position="dodge", stat="identity", colour = "black")+
  xlab("Delivery mode") +
  scale_fill_manual(values= c('#252525','#bdbdbd'))+
  scale_y_continuous(name ="Dyads",breaks=seq(0,100,20), expand = c(0, 0))+
  guides(fill=guide_legend("B. vulgatus\n sharing event")) + labs(subtitle = get_test_label(stat.test, detailed = TRUE))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
D1 

ggsave("Figure_5C1_i_MM.png", w = 5, h = 3)

#abx in labour sig species transferred#

stat.test <-chisq_test(df$ABX_in_labour,df$Bacteroides_vulgatus)

vulgatusABX <- table(df$ABX_in_labour, df$Bacteroides_vulgatus) 
vulgatusABX <-as.data.frame((vulgatusABX))
names(vulgatusABX) = c('ABX_in_labour', 'Bacteroides_vulgatus', 'Freq')

D2<-ggplot(vulgatusABX, aes(fill=ABX_in_labour, y=Freq, x=Bacteroides_vulgatus), size = 40) + 
  geom_bar(position="dodge", stat="identity", colour = "black")+
  xlab("Antibiotics during labour") +
  scale_fill_manual(values= c('#252525','#bdbdbd'))+
  scale_y_continuous(name ="Dyads",breaks=seq(0,100,20), expand = c(0, 0))+
  guides(fill=guide_legend("B. vulgatus\n sharing event")) + labs(subtitle = get_test_label(stat.test, detailed = TRUE))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
D2 

ggsave("Figire_5C2_ii_MM.png", w = 5, h = 3)


stat.test <-chisq_test(df$ABX_in_labour,df$Bifidobacterium_bifidum)

bifidumABX <- table(df$ABX_in_labour, df$Bifidobacterium_bifidum) 
bifidumABX <-as.data.frame((bifidumABX))
names(bifidumABX) = c('ABX_in_labour', 'Bifidobacterium_bifidum', 'Freq')

D3<-ggplot(bifidumABX, aes(fill=ABX_in_labour, y=Freq, x=Bifidobacterium_bifidum), size = 40) + 
  geom_bar(position="dodge", stat="identity", colour = "black")+
  xlab("Antibiotics during labour") +
  scale_fill_manual(values= c('#252525','#bdbdbd'))+
  scale_y_continuous(name ="Dyads",breaks=seq(0,100,20), expand = c(0, 0))+
  guides(fill=guide_legend("B. bifidum\n sharing event")) + labs(subtitle = get_test_label(stat.test, detailed = TRUE))+
  theme(
    panel.background = element_rect(fill = "white", colour = "black",
                                    size = 1, linetype = "solid"))
D3 

ggsave("Figure_5C_iii_MM.png", w = 5, h = 3)
}

