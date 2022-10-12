
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
  library(purrr)
}


#Figure 1 (a) was just made in powerpoint
#Figure 1 (b) was made in excel from the calculated relative abundances of Bifiodbacterium

###### Figure 1 (c) Heatmap ####
{
  df <- read.csv("csv_file_with_samples_as_columns_species_as_rows.csv", header = T)
  dim(df)
  #there were 1012 columns
  matrix <- as.matrix(df[,2:1012], "numeric")
  rownames(matrix)<-df[,1]
  
  
  col_fun <- colorRamp2(c(0, 10, 40, 50,60,70, 80,90, 100), c('#F2F2F2','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506','#54278f','#3f007d'))
  
  heatmap_meta_all <-read.csv("csv_file_where_first_column_are_same_sample_names_as_df_and_each_next_column_corresponding_metadata.csv")
  dim(heatmap_meta_all)
  
#only looked at 5 metadata classes
  heatmap_meta_all$Timepoint <- factor(heatmap_meta_all$Timepoint , levels=c("16_weeks", "34_weeks", "Delivery", "1_week", "1_month"))
  heatmap_meta_all$Sample_Type <- factor(heatmap_meta_all$Sample_Type , levels=c("Maternal_stool", "Infant_stool", "Maternal_vaginal", "Maternal_oral", "Maternal_EBM"))
  heatmap_meta_all$ABX_in_labour <- factor(heatmap_meta_all$ABX_in_labour , levels=c("Yes", "No"))
  heatmap_meta_all$Sharing_event <- factor(heatmap_meta_all$Sharing_event , levels=c("Yes", "No"))
  heatmap_meta_all$Delivery_mode <- factor(heatmap_meta_all$Delivery_mode , levels=c("Vaginal", "Caesarian"))
  
  ann <- data.frame(heatmap_meta_all$Sample_Type,
                    heatmap_meta_all$Timepoint,
                    heatmap_meta_all$Sharing_event,
                    heatmap_meta_all$Rupture_of_membranes,
                    heatmap_meta_all$Child_sex,
                    heatmap_meta_all$ABX_in_labour,
                    heatmap_meta_all$Delivery_mode)
  colnames(ann) <- c("Sample_Type","Timepoint",
                     "Sharing_event", "Rupture_of_membranes",
                     "Child_sex","ABX_in_labour",
                     "Delivery_mode")
  colours <- list("ABX_in_labour"=c("Yes"="#ef8a62",
                                    "No"="#252525"),
                  "Sharing_event"=c("Yes"="#fec44f",
                                    "No"="#252525"),
                  "Rupture_of_membranes"=c("POM"="#5ab4ac",
                                           "SROM"="#d8b365",
                                           "At LSCS"="#9ecae1",
                                           "ARM"="#6baed6",
                                           "SROM and ARM"="#252525"),
                  "Delivery_mode"=c("Vaginal"="#d9d9d9",
                                    "Caesarian"="#252525"),
                  "Child_sex"=c("Male"="#252525",
                                 "Female"="#542788"),
                  "Timepoint"=c("16_weeks" = "#fcc5c0",
                                "34_weeks" = "#fa9fb5",
                                "Delivery"="#f768a1",
                                "1_week"="#c51b8a",
                                "1_month"="#7a0177"),
                  "Sample_Type"=c("Infant_stool" = "#ce1256",
                                  "Maternal_EBM" = "#2b83ba",
                                  "Maternal_vaginal" = "#addd8e",
                                  "Maternal_oral" = "#abd9e9",
                                  "Maternal_stool" = "#238443"))

  
  
  alpha_meta <-heatmap_meta_all
  alpha_data <-t(df)
  colnames(alpha_data)<-alpha_data[1,]
  alpha_data <-alpha_data[-1,]
  write.csv(alpha_data, "alpha_data.csv")
  alpha_data<-read.csv("alpha_data.csv", header=T)
  alpha_data <- alpha_data[,-c(1:1)]
  shannon <- setNames(data.frame(alpha_meta,
                               diversity(alpha_data, index = 'shannon')),
                    c(colnames(alpha_meta), "Diversity"))
  shannon <-shannon[,-c(2:8)]
  shannon<- data.frame(shannon$Diversity)
  
  colAnn <- HeatmapAnnotation(df=ann, 'Shannon diversity' = anno_points(shannon, size = unit(0.15,"cm"), gp = gpar(col = "#228b22")), which="col", col=colours,
                              annotation_width=unit(c(1, 2), "cm"),
                              gap=unit(1, "mm"),
                              simple_anno_size = unit(0.3, "cm"))
  
  pdf("Fig_1c.pdf", width = 20, height = 12)
  png("Fig_1c.png",width=20,height=12,units="in",res=600)
  
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
                column_split = factor(heatmap_meta_all$Timepoint , levels=c("16_weeks", "34_weeks", "Delivery", "1_week", "1_month")),
                column_gap = unit(2, "mm"),
                cluster_column_slices = FALSE,
                heatmap_height = unit(25, "cm"),
                heatmap_width = unit(30, "cm"),
                top_annotation=colAnn, use_raster = TRUE, raster_quality = 5)
  
  draw(ht, heatmap_legend_side = "bottom")
  dev.off()
}

#Figure 1 (d) was made in excel
#Figure 1 (e) was plotted using iTOL


#Figure 1 (f)




######################################
#   Supplementary Figure 1 plots     #
######################################
  
# Supplementary Figure 1 (a) Heatmap #
  {
    df <- read.csv("csv_file_with_samples_as_columns_species_as_rows.csv", header = T)
    matrix <- as.matrix(df[,2:1012], "numeric")
    rownames(matrix)<-df[,1]
    
    
    col_fun <- colorRamp2(c(0, 10, 40, 50,60,70, 80,90, 100), c('#F2F2F2','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506','#54278f','#3f007d'))
    
    #hmp_strict_meta <- strict_meta[,c(1,2,6,7)]
    #write.csv(hmp_strict_meta,"top53_hmp_meta.csv")
    heatmap_meta_all <-read.csv("csv_file_where_first_column_are_same_sample_names_as_df_and_each_next_column_corresponding_metadata.csv")
    dim(heatmap_meta_all)
    
    heatmap_meta_all$Timepoint <- factor(heatmap_meta_all$Timepoint , levels=c("16_weeks", "34_weeks", "Delivery", "1_week", "1_month"))
    heatmap_meta_all$Sample_Type <- factor(heatmap_meta_all$Sample_Type , levels=c("Maternal_stool", "Infant_stool", "Maternal_vaginal", "Maternal_oral", "Maternal_EBM"))
    heatmap_meta_all$ABX_in_labour <- factor(heatmap_meta_all$ABX_in_labour , levels=c("Yes", "No"))
    heatmap_meta_all$Sharing_event <- factor(heatmap_meta_all$Sharing_event , levels=c("Yes", "No"))
    heatmap_meta_all$Delivery_mode <- factor(heatmap_meta_all$Delivery_mode , levels=c("Vaginal", "Caesarian"))
    
    ann <- data.frame(heatmap_meta_all$Sample_Type,
                      heatmap_meta_all$Timepoint,
                      heatmap_meta_all$Sharing_event,
                      heatmap_meta_all$Rupture_of_membranes,
                      heatmap_meta_all$Child_sex,
                      heatmap_meta_all$ABX_in_labour,
                      heatmap_meta_all$Delivery_mode)
    colnames(ann) <- c("Sample_Type","Timepoint",
                       "Sharing_event", "Rupture_of_membranes",
                       "Child_sex","ABX_in_labour",
                       "Delivery_mode")
    colours <- list("ABX_in_labour"=c("Yes"="#ef8a62",
                                      "No"="#252525"),
                    "Sharing_event"=c("Yes"="#fec44f",
                                      "No"="#252525"),
                    "Rupture_of_membranes"=c("POM"="#5ab4ac",
                                             "SROM"="#d8b365",
                                             "At LSCS"="#9ecae1",
                                             "ARM"="#6baed6",
                                             "SROM and ARM"="#252525"),
                    "Delivery_mode"=c("Vaginal"="#d9d9d9",
                                      "Caesarian"="#252525"),
                    "Child_sex"=c("Male"="#252525",
                                   "Female"="#542788"),
                    "Timepoint"=c("16_weeks" = "#fcc5c0",
                                  "34_weeks" = "#fa9fb5",
                                  "Delivery"="#f768a1",
                                  "1_week"="#c51b8a",
                                  "1_month"="#7a0177"),
                    "Sample_Type"=c("Infant_stool" = "#ce1256",
                                    "Maternal_EBM" = "#2b83ba",
                                    "Maternal_vaginal" = "#addd8e",
                                    "Maternal_oral" = "#abd9e9",
                                    "Maternal_stool" = "#238443"))
    
    alpha_meta <-heatmap_meta_all
    alpha_data <-t(df)
    colnames(alpha_data)<-alpha_data[1,]
    alpha_data <-alpha_data[-1,]
    write.csv(alpha_data, "top_alpha_data.csv")
    alpha_data<-read.csv("top_alpha_data.csv", header=T)
    alpha_data <- alpha_data[,-c(1:1)]
    shannon <- setNames(data.frame(alpha_meta,
                                   diversity(alpha_data, index = 'shannon')),
                        c(colnames(alpha_meta), "Diversity"))
    shannon <-shannon[,-c(2:8)]
    shannon<- data.frame(shannon$Diversity)
    
    colAnn <- HeatmapAnnotation(df=ann, 'Shannon diversity' = anno_points(shannon, size = unit(0.15,"cm"), gp = gpar(col = "#228b22")), which="col", col=colours,
                                annotation_width=unit(c(1, 2), "cm"),
                                gap=unit(1, "mm"),
                                simple_anno_size = unit(0.3, "cm"))
    

    
    pdf("FigS1_A.pdf", width = 15, height = 12)
    png("FigS1_A.png",width=15,height=12,units="in",res=1200)
    
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
                  column_split = factor(heatmap_meta_all$Timepoint , levels=c("16_weeks", "34_weeks", "Delivery", "1_week", "1_month")),
                  column_gap = unit(2, "mm"),
                  cluster_column_slices = FALSE,
                  heatmap_height = unit(25, "cm"),
                  heatmap_width = unit(30, "cm"),
                  top_annotation=colAnn, use_raster = TRUE, raster_quality = 5)
    
    draw(ht, heatmap_legend_side = "bottom")
    dev.off()
  }


#Supplementary Figure 1 (b) #

  
  alldata <- read.csv('Filtered_cleaned_HQ_species.csv', header = T)
  alldata_sample <- alldata[,-c(1:35)]
  alldata_meta <- alldata[,1:35]
  
  sampledata <- alldata %>% filter(Sample_Call == "Sample")
  sampledata_sample <- sampledata[,-c(1:35)]
  sampledata_meta <- sampledata[,1:35]
  write.csv(sampledata, "sampledata.csv")
  
  stooldata <- sampledata %>% filter(Sample_association == "Stool")
  stooldata_sample <- stooldata[,-c(1:35)]
  stooldata_meta <- stooldata[,1:35]
  write.csv(stooldata, "stooldata.csv")  

  
  
  alpha <- setNames(data.frame(stooldata_meta,
                               diversity(stooldata_sample, index = 'shannon')),
                    c(colnames(stooldata_meta), "Shannon"))
  
  alpha_l <- melt(alpha, id.vars = colnames(stooldata_meta), variable.name = 'Measure', value.name = 'Shannon Diversity')
  write.csv(alpha_l, 'stool_samples_alphas.csv')
  stooldata<- alpha_l %>% filter(Sample_Type == "Infant_stool")
  is.numeric(stooldata$Diversity)
  
  alpha_l$Timepoint_Source <- factor(alpha_l$Timepoint_Source , levels=c("Maternal_16_weeks", "Maternal_34_weeks","Maternal_1_month","Infant_Delivery", "Infant_1_week", "Infant_1_month"))
  alpha_l$Sample_Type <- factor(alpha_l$Sample_Type , levels=c("Maternal_stool", "Infant_stool"))
  
  my_comparisons <- list( c("Maternal_16_weeks", "Infant_Delivery"), c("Maternal_16_weeks", "Infant_1_week"), c("Maternal_16_weeks","Infant_1_month"),
                          c("Maternal_34_weeks","Infant_Delivery"), c("Maternal_34_weeks", "Infant_1_week"), c("Maternal_34_weeks","Infant_1_month"),
                          c("Maternal_1_month","Infant_Delivery"), c("Maternal_1_month", "Infant_1_week"), c("Maternal_1_month","Infant_1_month"))
  
  bxp <- ggboxplot( alpha_l, x = "Timepoint_Source", y = "Shannon Diversity", fill="Sample_Type") +
    scale_fill_manual(values = c("#636363", "#cccccc")) +
    theme(panel.background = element_rect(fill = "#e0e0e0",
                                          colour = "#e0e0e0",
                                          size = 0.25, linetype = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white")) +
    stat_compare_means(method = "anova", label.y = 1)+      
    stat_compare_means(label = "p.signif", method = "t.test",
                       ref.group = "Maternal_16_weeks") + theme(legend.position="bottom")+
    scale_x_discrete(labels=c("Maternal_16_weeks" = "16 Weeks", "Maternal_34_weeks" = "34 Weeks", "Maternal_1_month" = "1 Month",
                              "Infant_Delivery" = "Delivery", "Infant_1_week" = "1 Week", "Infant_1_month" = "1 Month"))+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    theme(
      panel.background = element_rect(fill = "white", colour = "black",
                                      size = 1, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      colour = "white"))
  bxp
  
  ggsave('FigS1_B.png', height = 6, width = 6, dpi = 300)
  

#Supplementary Figure1 (e)
  
  
  ##inputting data and arranging
  alldata <- read.csv('Filtered_cleaned_HQ_species.csv', header = T)
  
  
  ###################################################
  #Species samples
  ###################################################
  
  #Infant Delivery
  {
    infantstooldata <- alldata %>% filter(Family != "Control" & Source != "Maternal" & Timepoint_Source =="Infant_Delivery")
    #View(infantstooldata)
    #infantstooldata<-infantstooldata[-c(7,51,52,53,55,56,59),]
    all_infant_sample <- infantstooldata[,-c(1:35)]
    all_infant_sample<-all_infant_sample %>% 
      select_if(negate(function(col) is.numeric(col) && sum(col) < 0.000001))
    all_infant_sample <- model.matrix( ~ ., all_infant_sample)
    dim(all_infant_sample)
    
    
    all_infant_meta <- infantstooldata[,c(14:35)]
    
    infant_Del.mds <- metaMDS(all_infant_sample, distance = "bray", autotransform = FALSE)
    plot(infant_Del.mds)
    
    ordiplot(infant_Del.mds,type="n")
    orditorp(infant_Del.mds,display="species",col="red",air=0.01)
    orditorp(infant_Del.mds,display="sites",cex=1.25,air=0.01)
    
  
    envfit_infant_Del<- envfit(infant_Del.mds, all_infant_meta, perm=1000, na.rm=TRUE)
    envfit_infant_Del
    str(envfit_infant_Del)
    
    envInfDel_dfvars <- data.frame((envfit_infant_Del$vectors)$r, (envfit_infant_Del$vectors)$pvals)
    
    names(envInfDel_dfvars)[1] <- "R2"
    names(envInfDel_dfvars)[2] <- "pval"
    
    
    envInfDel_dffactors<- data.frame((envfit_infant_Del$factors)$r, (envfit_infant_Del$factors)$pvals)
    (names(envInfDel_dffactors)[1] <- "R2")
    (names(envInfDel_dffactors)[2] <- "pval")
    
    write.csv(rbind(envInfDel_dfvars, envInfDel_dffactors), "envInfdel_df.csv")
    infdelRs<-(rbind(envInfDel_dfvars, envInfDel_dffactors))
    

    
  }
  ##Infant 1 week
  {
    infantstooldata <- alldata %>% filter(Family != "Control" & Source != "Maternal" & Timepoint_Source =="Infant_1_week")
    all_infant_sample <- infantstooldata[,-c(1:35)]
    all_infant_meta <- infantstooldata[,c(14:35)]
    
    infant_1W.mds <- metaMDS(all_infant_sample, distance = "bray", autotransform = FALSE)
    plot(infant_1W.mds)
    
    ordiplot(infant_1W.mds,type="n")
    orditorp(infant_1W.mds,display="species",col="red",air=0.01)
    orditorp(infant_1W.mds,display="sites",cex=1.25,air=0.01)
    
    
    
    envfit_infant_1W<- envfit(infant_1W.mds, all_infant_meta, perm=1000, na.rm=TRUE)
    envfit_infant_1W
    str(envfit_infant_1W)
    
    envInf1W_dfvars <- data.frame((envfit_infant_1W$vectors)$r, (envfit_infant_1W$vectors)$pvals)
    
    names(envInf1W_dfvars)[1] <- "R2"
    names(envInf1W_dfvars)[2] <- "pval"
    
    
    envInf1W_dffactors<- data.frame((envfit_infant_1W$factors)$r, (envfit_infant_1W$factors)$pvals)
    (names(envInf1W_dffactors)[1] <- "R2")
    (names(envInf1W_dffactors)[2] <- "pval")
    
    write.csv(rbind(envInf1W_dfvars, envInf1W_dffactors), "envInf1W_df.csv")
    inf1wRs<-(rbind(envInf1W_dfvars, envInf1W_dffactors))
    

    
  }
  
  #Infant 1 month
  {
    infantstooldata <- alldata %>% filter(Family != "Control" & Source != "Maternal" & Timepoint_Source =="Infant_1_month")
    all_infant_sample <- infantstooldata[,-c(1:35)]
    all_infant_meta <- infantstooldata[,c(14:35)]
    
    infant_1M.mds <- metaMDS(all_infant_sample, distance = "bray", autotransform = FALSE)
    
    plot(infant_1M.mds)
    
    ordiplot(infant_1M.mds,type="n")
    orditorp(infant_1M.mds,display="species",col="red",air=0.01)
    orditorp(infant_1M.mds,display="sites",cex=1.25,air=0.01)
    
    
    
    envfit_infant_1M<- envfit(infant_1M.mds, all_infant_meta, perm=1000, na.rm=TRUE)
    envfit_infant_1M
    str(envfit_infant_1M)
    
    envInf1M_dfvars <- data.frame((envfit_infant_1M$vectors)$r, (envfit_infant_1M$vectors)$pvals)
    
    names(envInf1M_dfvars)[1] <- "R2"
    names(envInf1M_dfvars)[2] <- "pval"
    
    
    envInf1M_dffactors<- data.frame((envfit_infant_1M$factors)$r, (envfit_infant_1M$factors)$pvals)
    (names(envInf1M_dffactors)[1] <- "R2")
    (names(envInf1M_dffactors)[2] <- "pval")
    
    write.csv(rbind(envInf1M_dfvars, envInf1M_dffactors), "envInf1M_df.csv")
    inf1mRs<-(rbind(envInf1M_dfvars, envInf1M_dffactors))
    
    
    
  }
  
  #Maternal stool 16W
  {
    maternalstooldata <- alldata %>% filter(Family != "Control" & Source != "Infant" & Timepoint_Source =="Maternal_16_weeks" & Sample_association =="Stool")
    all_maternal_sample <- maternalstooldata[,-c(1:35)]
    all_maternal_meta <- maternalstooldata[,14:35]
    
    
    mum_16W.mds <- metaMDS(all_maternal_sample, distance = "bray", autotransform = FALSE)
    
    
    plot(mum_16W.mds)
    
    ordiplot(mum_16W.mds,type="n")
    orditorp(mum_16W.mds,display="species",col="red",air=0.01)
    orditorp(mum_16W.mds,display="sites",cex=1.25,air=0.01)
    
 
    
    envfit_mum_16W<- envfit(mum_16W.mds, all_maternal_meta, perm=1000, na.rm=TRUE)
    envfit_mum_16W
    str(envfit_mum_16W)
    
    envmum16W_dfvars <- data.frame((envfit_mum_16W$vectors)$r, (envfit_mum_16W$vectors)$pvals)
    
    names(envmum16W_dfvars)[1] <- "R2"
    names(envmum16W_dfvars)[2] <- "pval"
    
    
    envmum16W_dffactors<- data.frame((envfit_mum_16W$factors)$r, (envfit_mum_16W$factors)$pvals)
    (names(envmum16W_dffactors)[1] <- "R2")
    (names(envmum16W_dffactors)[2] <- "pval")
    
    write.csv(rbind(envmum16W_dfvars, envmum16W_dffactors), "envMum16W_df.csv")
    mum16Rs<-(rbind(envmum16W_dfvars, envmum16W_dffactors))
    
    
  }
  
  #Maternal stool 34w
  {
    maternalstooldata <- alldata %>% filter(Family != "Control" & Source != "Infant" & Timepoint_Source =="Maternal_34_weeks" & Sample_association =="Stool")
    all_maternal_sample <- maternalstooldata[,-c(1:35)]
    all_maternal_meta <- maternalstooldata[,14:35]
    
    mum_34W.mds <- metaMDS(all_maternal_sample, distance = "bray", autotransform = FALSE)
    
    
    plot(mum_34W.mds)
    
    ordiplot(mum_34W.mds,type="n")
    orditorp(mum_34W.mds,display="species",col="red",air=0.01)
    orditorp(mum_34W.mds,display="sites",cex=1.25,air=0.01)
    
    envfit_mum_34W<- envfit(mum_34W.mds, all_maternal_meta, perm=1000, na.rm=TRUE)
    envfit_mum_34W
    str(envfit_mum_34W)
    
    envmum34W_dfvars <- data.frame((envfit_mum_34W$vectors)$r, (envfit_mum_34W$vectors)$pvals)
    
    names(envmum34W_dfvars)[1] <- "R2"
    names(envmum34W_dfvars)[2] <- "pval"
    
    
    envmum34W_dffactors<- data.frame((envfit_mum_34W$factors)$r, (envfit_mum_34W$factors)$pvals)
    (names(envmum34W_dffactors)[1] <- "R2")
    (names(envmum34W_dffactors)[2] <- "pval")
    
    write.csv(rbind(envmum34W_dfvars, envmum34W_dffactors), "envMum34W_df.csv")
    mum34Rs<-(rbind(envmum34W_dfvars, envmum34W_dffactors))
    
    
  }
  
  #Maternal stool 1M
  {
    maternalstooldata <- alldata %>% filter(Family != "Control" & Source != "Infant" & Timepoint_Source =="Maternal_1_month" & Sample_association =="Stool")
    all_maternal_sample <- maternalstooldata[,-c(1:35)]
    all_maternal_meta <- maternalstooldata[,14:35]
    
    mum_1M.mds <- metaMDS(all_maternal_sample, distance = "bray", autotransform = FALSE)
    
    plot(mum_1M.mds)
    
    ordiplot(mum_1M.mds,type="n")
    orditorp(mum_1M.mds,display="species",col="red",air=0.01)
    orditorp(mum_1M.mds,display="sites",cex=1.25,air=0.01)
    
    envfit_mum_1M<- envfit(mum_1M.mds, all_maternal_meta, perm=1000, na.rm=TRUE)
    envfit_mum_1M
    str(envfit_mum_1M)
    
    envmum1M_dfvars <- data.frame((envfit_mum_1M$vectors)$r, (envfit_mum_1M$vectors)$pvals)
    
    names(envmum1M_dfvars)[1] <- "R2"
    names(envmum1M_dfvars)[2] <- "pval"
    
    
    envmum1M_dffactors<- data.frame((envfit_mum_1M$factors)$r, (envfit_mum_1M$factors)$pvals)
    (names(envmum1M_dffactors)[1] <- "R2")
    (names(envmum1M_dffactors)[2] <- "pval")
    
    write.csv(rbind(envmum1M_dfvars, envmum1M_dffactors), "envMum1M_df.csv")
    mum1Rs<-(rbind(envmum1M_dfvars, envmum1M_dffactors))
    
    
    
  }
  
  #Maternal EBM 1M
  {
    maternalEBMdata <- alldata %>% filter(Family != "Control" & Source != "Infant" & Timepoint_Source =="Maternal_1_month" & Sample_association =="EBM")
    all_maternal_sample <- maternalEBMdata[,-c(1:35)]
    all_maternal_meta <- maternalEBMdata[,c(14:35)]
    
    
    mum_EBM.mds <- metaMDS(all_maternal_sample, distance = "bray", autotransform = FALSE)
    
    plot(mum_EBM.mds)
    
    ordiplot(mum_EBM.mds,type="n")
    orditorp(mum_EBM.mds,display="species",col="red",air=0.01)
    orditorp(mum_EBM.mds,display="sites",cex=1.25,air=0.01)
    
    envfit_mum_EBM<- envfit(mum_EBM.mds, all_maternal_meta, perm=1000, na.rm=TRUE)
    envfit_mum_EBM
    str(envfit_mum_EBM)
    
    envmumEBM_dfvars <- data.frame((envfit_mum_EBM$vectors)$r, (envfit_mum_EBM$vectors)$pvals)
    
    names(envmumEBM_dfvars)[1] <- "R2"
    names(envmumEBM_dfvars)[2] <- "pval"
    
    
    envmumEBM_dffactors<- data.frame((envfit_mum_EBM$factors)$r, (envfit_mum_EBM$factors)$pvals)
    (names(envmumEBM_dffactors)[1] <- "R2")
    (names(envmumEBM_dffactors)[2] <- "pval")
    
    write.csv(rbind(envmumEBM_dfvars, envmumEBM_dffactors), "envMumEBM_df.csv")
    ebmRs<-(rbind(envmumEBM_dfvars, envmumEBM_dffactors))
   
    
    
  }
  
  #Maternal Vaginal 16W
  {
    maternalVaginal_16Wdata <- alldata %>% filter(Family != "Control" & Source != "Infant" & Timepoint_Source =="Maternal_16_weeks" & Sample_association =="Vaginal")
    all_maternal_sample <- maternalVaginal_16Wdata[,-c(1:35)]
    all_maternal_sample <- model.matrix( ~ ., all_maternal_sample)
    
    
    
    all_maternal_meta <- maternalVaginal_16Wdata[,14:35]
    
    
    
    mum_Vaginal_16W.mds <- metaMDS(all_maternal_sample, distance = "bray", autotransform = FALSE)
    envfit_mum_Vaginal_16W<- envfit(mum_Vaginal_16W.mds, all_maternal_meta, perm=1000, na.rm=TRUE)
    envfit_mum_Vaginal_16W
    str(envfit_mum_Vaginal_16W)
    
    envmumVaginal_16W_dfvars <- data.frame((envfit_mum_Vaginal_16W$vectors)$r, (envfit_mum_Vaginal_16W$vectors)$pvals)
    
    names(envmumVaginal_16W_dfvars)[1] <- "R2"
    names(envmumVaginal_16W_dfvars)[2] <- "pval"
    
    
    envmumVaginal_16W_dffactors<- data.frame((envfit_mum_Vaginal_16W$factors)$r, (envfit_mum_Vaginal_16W$factors)$pvals)
    (names(envmumVaginal_16W_dffactors)[1] <- "R2")
    (names(envmumVaginal_16W_dffactors)[2] <- "pval")
    
    write.csv(rbind(envmumVaginal_16W_dfvars, envmumVaginal_16W_dffactors), "envMumVaginal_16W_df.csv")
    vag16Rs<-(rbind(envmumVaginal_16W_dfvars, envmumVaginal_16W_dffactors))
   
    
  }
  
  #Maternal Vaginal 34W
  {
    maternalVaginal_34Wdata <- alldata %>% filter(Family != "Control" & Source != "Infant" & Timepoint_Source =="Maternal_34_weeks" & Sample_association =="Vaginal")
    all_maternal_sample <- maternalVaginal_34Wdata[,-c(1:35)]
    all_maternal_meta <- maternalVaginal_34Wdata[,14:35]
    
    
    mum_Vaginal_34W.mds <- metaMDS(all_maternal_sample, distance = "bray", autotransform = FALSE)
    envfit_mum_Vaginal_34W<- envfit(mum_Vaginal_34W.mds, all_maternal_meta, perm=1000, na.rm=TRUE)
    envfit_mum_Vaginal_34W
    str(envfit_mum_Vaginal_34W)
    
    envmumVaginal_34W_dfvars <- data.frame((envfit_mum_Vaginal_34W$vectors)$r, (envfit_mum_Vaginal_34W$vectors)$pvals)
    
    names(envmumVaginal_34W_dfvars)[1] <- "R2"
    names(envmumVaginal_34W_dfvars)[2] <- "pval"
    
    
    envmumVaginal_34W_dffactors<- data.frame((envfit_mum_Vaginal_34W$factors)$r, (envfit_mum_Vaginal_34W$factors)$pvals)
    (names(envmumVaginal_34W_dffactors)[1] <- "R2")
    (names(envmumVaginal_34W_dffactors)[2] <- "pval")
    
    vag34Rs<-(rbind(envmumVaginal_34W_dfvars, envmumVaginal_34W_dffactors))
  
  }
  
  
  #Maternal Oral 16W
  {
    maternalOral_16Wdata <- alldata %>% filter(Family != "Control" & Source != "Infant" & Timepoint_Source =="Maternal_16_weeks" & Sample_association =="Oral")
    all_maternal_sample <- maternalOral_16Wdata[,-c(1:35)]
    all_maternal_meta <- maternalOral_16Wdata[,14:35]
    
    mum_Oral_16W.mds <- metaMDS(all_maternal_sample, distance = "bray", autotransform = FALSE)
    envfit_mum_Oral_16W<- envfit(mum_Oral_16W.mds, all_maternal_meta, perm=1000, na.rm=TRUE)
    envfit_mum_Oral_16W
    str(envfit_mum_Oral_16W)
    
    envmumOral_16W_dfvars <- data.frame((envfit_mum_Oral_16W$vectors)$r, (envfit_mum_Oral_16W$vectors)$pvals)
    
    names(envmumOral_16W_dfvars)[1] <- "R2"
    names(envmumOral_16W_dfvars)[2] <- "pval"
    
    
    envmumOral_16W_dffactors<- data.frame((envfit_mum_Oral_16W$factors)$r, (envfit_mum_Oral_16W$factors)$pvals)
    (names(envmumOral_16W_dffactors)[1] <- "R2")
    (names(envmumOral_16W_dffactors)[2] <- "pval")
    
    write.csv(rbind(envmumOral_16W_dfvars, envmumOral_16W_dffactors), "envMumOral_16W_df.csv")
    oralRs<-(rbind(envmumOral_16W_dfvars, envmumOral_16W_dffactors))
    
    
    
  }
  
  
  ###########################################
  ##plotting ENVFIT
  #########################################
  csvs <- c('infdelRs','inf1wRs','inf1mRs',
            'mum16Rs','mum34Rs','mum1Rs',
            'ebmRs',
            'vag16Rs','vag34Rs',
            'oralRs')
  
  
  my_fun <- function(csvs){
    tmp <- paste0(csvs,'$period <- \'', csvs, '\'',sep = "")
    eval(parse(text = tmp))
    df <- eval(parse(text=csvs))
    return(df)
  }
  
  dfs <- lapply(csvs, FUN = my_fun)
  names(dfs) <- csvs
  all_envfit_results<-as.data.frame(lapply(names(dfs), function(x) assign(x, dfs[[x]], envir = .GlobalEnv)))
  #View(all_envfit_results)
  
  
  covariate_categories <- data.frame("Covariate" =c('Birthweight',
                                                    'Gestation days',
                                                    'Labour length min',
                                                    'Maternal age at recruitment',
                                                    'Placental weight',
                                                    'ABX in labour',
                                                    'ABX in NICU',
                                                    'Blood_group',
                                                    'BMI class',
                                                    'Breast feeding class short',
                                                    'Child sex',
                                                    'Delivery mode',
                                                    'HP index label',
                                                    'Infant ABX in 1 month',
                                                    'Labour onset',
                                                    'Lewis status',
                                                    'Maternal ABX in 1 month',
                                                    'NICU admission',
                                                    'Parity category',
                                                    'Rupture of membranes',
                                                    'Secretor status',
                                                    'Trial group'),
                                     "Covariate type" = c('Infant',
                                                          'Maternal Perinatal',
                                                          'Maternal Perinatal',
                                                          'Maternal Health & Lifestyle',
                                                          'Infant',
                                                          'Maternal Perinatal',
                                                          'Infant',
                                                          'Maternal Health & Lifestyle',
                                                          'Maternal Health & Lifestyle',
                                                          'Infant',
                                                          'Infant',
                                                          'Maternal Perinatal',
                                                          'Maternal Health & Lifestyle',
                                                          'Infant',
                                                          'Maternal Perinatal',
                                                          'Maternal Health & Lifestyle',
                                                          'Maternal Health & Lifestyle',
                                                          'Infant',
                                                          'Maternal Perinatal',
                                                          'Maternal Perinatal',
                                                          'Maternal Health & Lifestyle',
                                                          'Maternal Dietary'))
  
  all_envfit_results<-cbind(covariate_categories,all_envfit_results)
  write.csv(all_envfit_results,"envfitPvalues.csv")
  
  names(all_envfit_results)
  all_envfit_R2 <- all_envfit_results[ -c(4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32) ]
  colnames(all_envfit_R2)<- c('Covariate',
                              'Covariate_type',
                              'Infant_Delivery',
                              'Infant_1W',
                              'Infant_1_month',
                              'Maternal_stool_16W',
                              'Maternal_stool_34W',
                              'Maternal_stool_1_month',
                              'Maternal_EBM_1_month',
                              'Maternal_HVS_16W',
                              'Maternal_HVS_34W',
                              'Maternal_OralRinse_16W')
  names(all_envfit_R2)
  
  envfitR2_long <- gather(all_envfit_R2, Sample_type, R2, 'Infant_Delivery':'Maternal_OralRinse_16W', factor_key=TRUE)
  envfitR2_long$Covariate_type = factor(envfitR2_long$Covariate_type, levels=c('Infant','Maternal Health & Lifestyle', 'Maternal Perinatal', 'Maternal Dietary'))
  envfitR2_long$Sample_type = factor(envfitR2_long$Sample_type, levels=c('Maternal_stool_16W','Maternal_stool_34W','Maternal_stool_1_month',
                                                                         'Maternal_HVS_16W','Maternal_HVS_34W', 'Maternal_OralRinse_16W','Maternal_EBM_1_month',
                                                                         'Infant_Delivery', 'Infant_1W', 'Infant_1_month'))
  head(envfitR2_long)
  
  colPalette<-c("#FC8D59", "#998EC3", "#91BFDB", "#810F7C")
  
  
  
  element_textbox <- function(...) {
    el <- element_text(...)
    class(el) <- c("element_textbox", class(el))
    el
  }
  
  element_grob.element_textbox <- function(element, ...) {
    text_grob <- NextMethod()
    rect_grob <- element_grob(calc_element("strip.background", theme_grey()))
    
    ggplot2:::absoluteGrob(
      grid::gList(
        element_grob(calc_element("strip.background", theme_grey())),
        text_grob
      ),
      height = grid::grobHeight(text_grob), 
      width = grid::unit(1, "npc")
    )
  }
  
  
  
  envplots<-ggplot(envfitR2_long, aes(x = Covariate, y = R2, fill = Covariate_type))+
    facet_wrap(~Sample_type, nrow = 1) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values=colPalette) +
    theme (axis.text.x = element_text(size = 15, colour = "black", angle = 45)) +
    theme (axis.text.y = element_text(size = 20, colour = "black"))+
    theme(axis.title = element_text(size = 20)) +
    theme(legend.position = "bottom",
          legend.text = element_text(color = "black", size = 20),
          legend.title = element_text(color = "black", size = 22)) +
    theme(strip.text = element_text(face="bold", size=20),
          strip.background = element_rect(size=2)) +
    theme(panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.25, linetype = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "white")) +
    xlab('Covariate')+
    ylab(bquote(''~R^2)) +
    coord_flip()+ labs( title= "Relationships with taxonomy")+
    theme(plot.title = element_textbox(size =30, hjust = 0.5, margin = margin(t = 5, b = 5)))
  
  
  
  library(egg)
  #birthweight
  BW_tag <- c("", "","","","","","","*","","")
  
  sig_envplots<-tag_facet(envplots,
                          x = 0.9, y = 0.18,
                          vjust = -1, hjust = -0.25,
                          open = "", close = "",
                          fontface = 4,
                          size = 15,
                          color = "red",
                          family = "serif",
                          tag_pool = BW_tag)
  
  #abx_in_labour
  ABXLAB_tag <- c("", "","","","","","","","*","")
  
  sig_envplots<-tag_facet(sig_envplots, 
                          x = 0, y = 0.3, 
                          vjust = -1, hjust = -0.25,
                          open = "", close = "",
                          fontface = 4,
                          size = 15,
                          color = "red",
                          family = "serif",
                          tag_pool = ABXLAB_tag)
  
  #lewis factor
  lewis_tag <- c("", "","","","","","","","","*")
  
  sig_envplots<-tag_facet(sig_envplots, 
                          x = 11.8, y = 0.06, 
                          vjust = -1, hjust = -0.25,
                          open = "", close = "",
                          fontface = 4,
                          size = 15,
                          color = "red",
                          family = "serif",
                          tag_pool = lewis_tag)
  
  
  #placental weight
  PLAC_tag <- c("", "","*","","","","","","","")
  
  sig_envplots<-tag_facet(sig_envplots, 
                          x = 16.6, y = 0.13, 
                          vjust = -1, hjust = -0.25,
                          open = "", close = "",
                          fontface = 4,
                          size = 15,
                          color = "red",
                          family = "serif",
                          tag_pool = PLAC_tag)
  
  
  #parity
  parity_tag <- c("", "","","*","*","","","*","","")
  
  sig_envplots<-tag_facet(sig_envplots, 
                          x = 15.8, y = 0.12, 
                          vjust = -1, hjust = -0.25,
                          open = "", close = "",
                          fontface = 4,
                          size = 15,
                          color = "red",
                          family = "serif",
                          tag_pool = parity_tag)
  
  #BMI class
  BMI_tag <- c("", "","","","","","*","","","")
  
  sig_envplots<-tag_facet(sig_envplots, 
                          x = 2.9, y = 0.11, 
                          vjust = -1, hjust = -0.25,
                          open = "", close = "",
                          fontface = 4,
                          size = 15,
                          color = "red",
                          family = "serif",
                          tag_pool = BMI_tag)
  
  #Breast feeding
  BREAST_tag <- c("", "*","","","","","","","","")
  
  sig_envplots<-tag_facet(sig_envplots, 
                          x = 3.6, y = 0.07, 
                          vjust = -1, hjust = -0.25,
                          open = "", close = "",
                          fontface = 4,
                          size = 15,
                          color = "red",
                          family = "serif",
                          tag_pool = BREAST_tag)
  
  
  #Delivery mode
  DELIVERY_tag <- c("", "","*","","","","","","","")
  
  sig_envplots<-tag_facet(sig_envplots, 
                          x = 5.9, y = 0.1, 
                          vjust = -1, hjust = -0.25,
                          open = "", close = "",
                          fontface = 4,
                          size = 15,
                          color = "red",
                          family = "serif",
                          tag_pool = DELIVERY_tag)
  
  #INFANT ABX
  ABX_tag <- c("", "","","","*","","","","","")
  
  sig_envplots<-tag_facet(sig_envplots, 
                          x = 8.9, y = 0.1, 
                          vjust = -1, hjust = -0.25,
                          open = "", close = "",
                          fontface = 4,
                          size = 15,
                          color = "red",
                          family = "serif",
                          tag_pool = ABX_tag)
  
  
  sig_envplots<-sig_envplots + theme(strip.text = element_text(face="bold", size=20),
                                     strip.background = element_rect(size=2))
  sig_envplots
  ggsave('FigS1_E.png',width=60,height=9,units="in", limitsize=FALSE)
  
  
  
  