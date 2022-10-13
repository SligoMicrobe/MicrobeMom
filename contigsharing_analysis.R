library('tidyverse')

##### Counting the number of contigs screened for transmission/sharing
infantsamplelist <- list.files(path = 'kraken_outputs/', pattern = 'PC.*_output_named_formatted.txt', full.names = F) %>%
  str_replace(., '_output_named_formatted.txt', '')
infantcontigs <- data.frame()
for (i in infantsamplelist){
  df <- read.delim(paste0('kraken_outputs/', i, '_output_named_formatted.txt'), header = F, sep = '\t') %>%
    mutate(V2 = str_replace(V2, '_cov.*', '')) %>%
    mutate(V2 = str_replace(V2, '_length', '_')) %>%
    mutate(V2 = paste0(i, '__', V2)) %>%
    pull(V2)
  
  infantcontigs <- append(infantcontigs, df)
}

infantcontigs <- infantcontigs %>%
  unlist() %>%
  data.frame() %>%
  separate('.', into = c('Sample', 'Contig', 'Length'), sep = '__')

infantcontigs %>%
  filter(Length > 1500) %>%
  dim()

# Counting the number of contigs deemed to be shared between mother and infant
trans_events <- data.frame()
dyadlist <- list.files(path = 'Blast_Outputs/', pattern = '.*_infant_vsmaternal_blastout.txt', full.names = F) %>%
  str_replace(., '_infant_vsmaternal_blastout.txt', '')
samplelist <- read.delim('samplelist.txt', header = T) %>%
  data.frame() 
for (i in dyadlist){
  if(file.info(paste0('Blast_Outputs/', i, '_infant_vsmaternal_blastout.txt'))$size == 0){next}
  trans <- read.delim(paste0('Blast_Outputs/', i, '_infant_vsmaternal_blastout.txt'), header = F) %>%
    select(V1, V2, V4, V5) %>%
    rename('Infant_Contig' = 'V1', 'Maternal_Contig' = 'V2', 'PCov' = 'V4', 'PID' = 'V5') %>%
    filter(PID >= 99.82928) %>% # BASED THIS CUTOFF ON THE RELATIONSHIP BETWEEN CONANI AND POPANI FROM INSTRAIN OUTPUTS
    select(Infant_Contig, Maternal_Contig) %>%
    distinct()
  
  if (dim(trans)[1] > 0){
    trans_events <- rbind(trans_events, trans)
  }
}

trans_events %>%
  mutate(MaternalSource = str_sub(Maternal_Contig, 1, 2)) %>%
  group_by(MaternalSource) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = 100*(Count/sum(Count)))

contigsharing_families <- trans_events %>%
  mutate(Family = str_sub(Infant_Contig, 3, 5)) %>%
  pull(Family) %>%
  unique()
100*(length(contigsharing_families)/length(dyadlist))

faecalcontigsharing_families <- trans_events %>%
  filter(str_sub(Maternal_Contig, 1, 2) == 'PB') %>%
  mutate(Family = str_sub(Infant_Contig, 3, 5)) %>%
  pull(Family) %>%
  unique()
length(faecalcontigsharing_families)
100*(length(faecalcontigsharing_families)/length(dyadlist))

vaginalcontigsharing_families <- trans_events %>%
  filter(str_sub(Maternal_Contig, 1, 2) == 'PV') %>%
  mutate(Family = str_sub(Infant_Contig, 3, 5)) %>%
  pull(Family) %>%
  unique()
length(vaginalcontigsharing_families)
100*(length(vaginalcontigsharing_families)/length(dyadlist))

milkcontigsharing_families <- trans_events %>%
  filter(str_sub(Maternal_Contig, 1, 2) == 'PM') %>%
  mutate(Family = str_sub(Infant_Contig, 3, 5)) %>%
  pull(Family) %>%
  unique()
length(milkcontigsharing_families)
100*(length(milkcontigsharing_families)/length(dyadlist))

oralcontigsharing_families <- trans_events %>%
  filter(str_sub(Maternal_Contig, 1, 2) == 'PO') %>%
  mutate(Family = str_sub(Infant_Contig, 3, 5)) %>%
  pull(Family) %>%
  unique()
length(oralcontigsharing_families)
100*(length(oralcontigsharing_families)/length(dyadlist))

VennDiagram::venn.diagram(x = list(faecalcontigsharing_families,
                                   vaginalcontigsharing_families,
                                   oralcontigsharing_families,
                                   milkcontigsharing_families),
                          category.names = c('Faecal', 'Vaginal', 'Oral', 'Breast Milk'),
                          filename = 'venn_contigsharing.png',
                          width = 1800,
                          height = 1500)

library('ComplexHeatmap')
library('UpSetR')

input <- list(Faecal = faecalcontigsharing_families,
              Vaginal = vaginalcontigsharing_families,
              Oral = oralcontigsharing_families,
              `Breast Milk` = milkcontigsharing_families)

upset_sharing <- upset(fromList(input),
           order.by = 'freq',
           mainbar.y.label = "Number of Dyads",
           sets.x.label = "Dyads",
           point.size = 3.5, 
           text.scale = 1.5,
           sets.bar.color = c('#006a4e',
                              '#d4af37',
                              '#1167b1',
                              '#522d80'))

pdf('upset_contigsharing.pdf', width = 6, height = 3.2)
upset_sharing
grid::grid.text('Dyads Exhibiting Contig Sharing', x = 0.65, y = 0.95, gp = gpar(fontsize = 10))
dev.off()

# Redoing the analysis from above but only including contigs where the maternal samples were pre-delivery (definite transmission)
vert_trans_events <- trans_events %>%
  mutate(MaternalSource = str_sub(Maternal_Contig, 1, 2)) %>%
  mutate(MaternalTimepoint = str_replace(Maternal_Contig, '_NODE.*', '')) %>%
  mutate(MaternalTimepoint = str_replace(MaternalTimepoint, '.*_', '')) %>%
  filter(MaternalTimepoint %in% c('E', 'L'))

vert_trans_events %>%
  group_by(MaternalSource) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = 100*(Count/sum(Count)))

contigsharing_families <- vert_trans_events %>%
  mutate(Family = str_sub(Infant_Contig, 3, 5)) %>%
  pull(Family) %>%
  unique()
100*(length(contigsharing_families)/length(dyadlist))

faecalcontigsharing_families <- vert_trans_events %>%
  filter(str_sub(Maternal_Contig, 1, 2) == 'PB') %>%
  mutate(Family = str_sub(Infant_Contig, 3, 5)) %>%
  pull(Family) %>%
  unique()
length(faecalcontigsharing_families)
100*(length(faecalcontigsharing_families)/length(dyadlist))

vaginalcontigsharing_families <- vert_trans_events %>%
  filter(str_sub(Maternal_Contig, 1, 2) == 'PV') %>%
  mutate(Family = str_sub(Infant_Contig, 3, 5)) %>%
  pull(Family) %>%
  unique()
length(vaginalcontigsharing_families)
100*(length(vaginalcontigsharing_families)/length(dyadlist))

milkcontigsharing_families <- vert_trans_events %>%
  filter(str_sub(Maternal_Contig, 1, 2) == 'PM') %>%
  mutate(Family = str_sub(Infant_Contig, 3, 5)) %>%
  pull(Family) %>%
  unique()
length(milkcontigsharing_families)
100*(length(milkcontigsharing_families)/length(dyadlist))

oralcontigsharing_families <- vert_trans_events %>%
  filter(str_sub(Maternal_Contig, 1, 2) == 'PO') %>%
  mutate(Family = str_sub(Infant_Contig, 3, 5)) %>%
  pull(Family) %>%
  unique()
length(oralcontigsharing_families)
100*(length(oralcontigsharing_families)/length(dyadlist))

VennDiagram::venn.diagram(x = list(faecalcontigsharing_families,
                                   vaginalcontigsharing_families,
                                   oralcontigsharing_families),
                          category.names = c('Faecal', 'Vaginal', 'Oral'),
                          filename = 'venn_contigtranmission.png',
                          width = 1800,
                          height = 1500)

library('ComplexHeatmap')
library('UpSetR')

input <- list(Faecal = faecalcontigsharing_families,
              Vaginal = vaginalcontigsharing_families,
              Oral = oralcontigsharing_families,
              `Breast Milk` = milkcontigsharing_families)

upset_transmission <- upset(fromList(input),
                            order.by = 'freq',
                            mainbar.y.label = "Number of Dyads",
                            sets.x.label = "Dyads",
                            point.size = 3.5, 
                            text.scale = 1.5,
                            sets.bar.color = c('#006a4e',
                                               '#d4af37',
                                               '#522d80'))

pdf('upset_contigtransmission.pdf', width = 6, height = 3.2)
upset_transmission
grid::grid.text('Dyads Exhibiting Contig Transfer', x = 0.65, y = 0.95, gp = gpar(fontsize = 10))
dev.off()


##### COMBINING BLAST AND ABRICATE OUTPUTS TO ID SHARED CONTIGS WITH FEATURES OF INTEREST
abricate_trans <- data.frame()
trans_events <- data.frame()
dyadlist <- list.files(path = 'Blast_Outputs/', pattern = '.*_infant_vsmaternal_blastout.txt', full.names = F) %>%
  str_replace(., '_infant_vsmaternal_blastout.txt', '')
samplelist <- read.delim('samplelist.txt', header = T) %>%
  data.frame() 

for (i in dyadlist){
  if(file.info(paste0('Blast_Outputs/', i, '_infant_vsmaternal_blastout.txt'))$size == 0){next}
  trans <- read.delim(paste0('Blast_Outputs/', i, '_infant_vsmaternal_blastout.txt'), header = F) %>%
    select(V1, V2, V4, V5) %>%
    rename('Infant_Contig' = 'V1', 'Maternal_Contig' = 'V2', 'PCov' = 'V4', 'PID' = 'V5') %>%
    filter(PID >= 99.82928) %>% # BASED THIS CUTOFF ON THE RELATIONSHIP BETWEEN CONANI AND POPANI FROM INSTRAIN OUTPUTS
    select(Infant_Contig, Maternal_Contig) %>%
    distinct()
  
  if (dim(trans)[1] > 0){
    trans_events <- rbind(trans_events, trans)
  }
    
  amr <- read.delim(paste0('Blast_Outputs/', i, '_card.txt'), header = T, check.names = F) %>%
    filter(`%COVERAGE` >= 80)
  amr_trans <- amr %>%
    filter(SEQUENCE %in% trans$Infant_Contig)
    
  vf <- read.delim(paste0('Blast_Outputs/', i, '_vfdb.txt'), header = T, check.names = F) %>%
    filter(`%COVERAGE` >= 80)
  vf_trans <- vf %>%
    filter(SEQUENCE %in% trans$Infant_Contig)
    
  plasmid <- read.delim(paste0('Blast_Outputs/', i, '_plasmid.txt'), header = T, check.names = F) %>%
    filter(`%COVERAGE` >= 80)
  plasmid_trans <- plasmid %>%
    filter(SEQUENCE %in% trans$Infant_Contig)
    
  if (dim(rbind(amr_trans, vf_trans, plasmid_trans))[1] > 0){
    abricate_trans <- rbind(abricate_trans,
                            rbind(amr_trans, vf_trans, plasmid_trans))
  }
}

trans_events <- trans_events %>%
  mutate(Dyad = substr(Infant_Contig, 3, 5)) %>%
  mutate(TimeSource_Pair = paste0(gsub('_NODE.*', '', Infant_Contig), '__', gsub('_NODE.*', '', Maternal_Contig))) %>%
  mutate(TimeSource_Pair = gsub('[0-9][0-9][0-9]_', '', TimeSource_Pair)) %>%
  mutate(DyadTimeSource_Pair = paste0(Dyad, '__', TimeSource_Pair)) %>%
  select(-Dyad, -TimeSource_Pair) %>%
  relocate(DyadTimeSource_Pair, .before = Infant_Contig)

abricate_trans <- abricate_trans %>%
  mutate(`#FILE` = gsub('_allcontigs.fasta' ,'', `#FILE`)) %>%
  rename('Dyad' = `#FILE`)

samplepairs <- samplelist %>%
  rename('Sample1' = Sample) %>%
  mutate(Sample2 = Sample1) %>%
  expand(Sample1, Sample2) %>%
  filter(Sample1 != Sample2) %>%
  filter((str_sub(Sample1, 1, 2) == 'PC') & (str_sub(Sample2, 1, 2) != 'PC')) %>%
  filter(str_sub(Sample1, 3, 5) == str_sub(Sample2, 3, 5)) %>%
  mutate(DyadTimeSource_Pair = paste0(str_sub(Sample1, 3, 5),
                                      '__',
                                      str_replace(Sample1, '[0-9][0-9][0-9]_', ''),
                                      '__',
                                      str_replace(Sample2, '[0-9][0-9][0-9]_', ''))) %>%
  relocate(DyadTimeSource_Pair, .before = Sample1) %>%
  select(DyadTimeSource_Pair)

trans_events %>%
  group_by(DyadTimeSource_Pair) %>%
  summarise(Events = n()) %>%
  arrange(DyadTimeSource_Pair) %>%
  right_join(samplepairs) %>%
  write.table('sharing_events.tsv', quote = F, row.names = F, sep = '\t')

trans_events %>%
  group_by(DyadTimeSource_Pair) %>%
  summarise(Events = n()) %>%
  arrange(DyadTimeSource_Pair) %>%
  right_join(samplepairs) %>%
  mutate(Events = replace_na(Events, 0)) %>%
  mutate(TimeSource_Pair = str_replace(DyadTimeSource_Pair, '^[0-9][0-9][0-9]__', '')) %>%
  group_by(TimeSource_Pair) %>%
  summarise(Mean = mean(Events),
            Median = median(Events),
            SD = sd(Events)) %>%
  arrange(-Mean) %>%
  write.table('sharing_events_summary.tsv', quote = F, row.names = F, sep = '\t')

abricate_trans %>%
  write.table('sharing_abricate.tsv', quote = F, row.names = F, sep = '\t')

##### PLOTTING CONTIG SHARING RESULTS BETWEEN SAMPLE ORIGINS (SOURCE AND TIMEPOINT)
sample_colours <- c('Faeces Early' =	'#b8d4cd',
                    'Faeces Late' =	'#5ca08e',
                    'Faeces 1 Month' =	'#006a4e',
                    'Infant Faeces Delivery' =	'#faa19b',
                    'Infant Faeces 1 Week' =	'#f66257',
                    'Infant Faeces 1 Month' =	'#d63b2f',
                    'Vaginal Swab Early' =	'#ae87d0',
                    'Vaginal Swab Late' =	'#522d80',
                    'Breast Milk 1 Month' =	'#1167b1',
                    'Oral Rinse Early' =	'#d4af37')

events <- read.delim('sharing_events_summary.tsv', header = T, check.names = F, sep = '\t') %>%
  separate(TimeSource_Pair, into = c('Infant_Sample', 'Maternal_Sample'), sep = '__') %>%
  mutate(Maternal_Sample = gsub('^PB', 'Faeces ', Maternal_Sample)) %>%
  mutate(Maternal_Sample = gsub('^PM', 'Breast Milk ', Maternal_Sample)) %>%
  mutate(Maternal_Sample = gsub('^PO', 'Oral Rinse ', Maternal_Sample)) %>%
  mutate(Maternal_Sample = gsub('^PV', 'Vaginal Swab ', Maternal_Sample)) %>%
  mutate(Infant_Sample = gsub('^PC', 'Infant Faeces ', Infant_Sample)) %>%
  mutate(Maternal_Sample = gsub(' E$', ' Early', Maternal_Sample)) %>%
  mutate(Maternal_Sample = gsub(' L$', ' Late', Maternal_Sample)) %>%
  mutate(Maternal_Sample = gsub(' 1M$', ' 1 Month', Maternal_Sample)) %>%
  mutate(Infant_Sample = gsub(' DEL$', ' Delivery', Infant_Sample)) %>%
  mutate(Infant_Sample = gsub(' 1W$', ' 1 Week', Infant_Sample)) %>%
  mutate(Infant_Sample = gsub(' 1M$', ' 1 Month', Infant_Sample))

p_cont <- events %>%
  mutate(Infant_Sample = factor(Infant_Sample, levels = c('Infant Faeces Delivery',
                                                          'Infant Faeces 1 Week',
                                                          'Infant Faeces 1 Month'))) %>%
  mutate(Maternal_Sample = factor(Maternal_Sample, levels = c('Faeces Early',
                                                             'Faeces Late',
                                                             'Faeces 1 Month',
                                                             'Breast Milk 1 Month',
                                                             'Vaginal Swab Early',
                                                             'Vaginal Swab Late',
                                                             'Oral Rinse Early'))) %>%
  ggplot(aes(x = forcats::fct_rev(Maternal_Sample), y = Mean)) +
  geom_bar(aes(fill = Maternal_Sample), stat = 'identity', colour = 'black', show.legend = F) +
  coord_flip() +
  ggforce::facet_col(~Infant_Sample, space = 'free') +
  scale_fill_manual(values = sample_colours) +
  labs(x = '', y = 'Mean Contig-Sharing Events') +
  theme(panel.background = element_rect(fill = NA, colour = 'black'),
        strip.text = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold'),
        panel.grid = element_blank())

pdf('sharing_events.pdf', width = 6, height = 5)
p_cont
dev.off()


##### SUMMARISING SHARED GENES OF INTEREST
events <- read.delim('sharing_events.tsv', header = T, check.names = F, sep = '\t') %>%
  mutate(Dyad = str_sub(DyadTimeSource_Pair, 1, 3))
annot <- read.delim('sharing_abricate.tsv', header = T, check.names = F, sep = '\t') %>%
  select(Dyad, SEQUENCE, GENE, `%COVERAGE`, `%IDENTITY`, DATABASE) %>%
  mutate(Dyad = str_sub(Dyad, 1, 3))

annot %>%
  filter(DATABASE == 'card') %>%
  group_by(GENE) %>%
  summarise(Count = length(GENE)) %>%
  arrange(-Count) %>%
  write.table('transmitted_genes_card.tsv', quote = F, row.names = F, sep = '\t')

annot %>%
  filter(DATABASE == 'vfdb') %>%
  group_by(GENE) %>%
  summarise(Count = length(GENE)) %>%
  arrange(-Count) %>%
  write.table('transmitted_genes_vf.tsv', quote = F, row.names = F, sep = '\t')

annot %>%
  filter(DATABASE == 'plasmidfinder') %>%
  group_by(GENE) %>%
  summarise(Count = length(GENE)) %>%
  arrange(-Count) %>%
  write.table('transmitted_genes_plasmid.tsv', quote = F, row.names = F, sep = '\t')

##### PLOTTING TOP SHARED GENES OF INTEREST
events <- read.delim('sharing_events.tsv', header = T, check.names = F, sep = '\t') %>%
  mutate(Dyad = str_sub(DyadTimeSource_Pair, 1, 3))
annot <- read.delim('sharing_abricate.tsv', header = T, check.names = F, sep = '\t') %>%
  select(Dyad, SEQUENCE, GENE, `%COVERAGE`, `%IDENTITY`, DATABASE) %>%
  mutate(Dyad = str_sub(Dyad, 1, 3))
p <- annot %>%
  filter(DATABASE == 'card') %>%
  group_by(GENE) %>%
  summarise(Count = length(GENE)) %>%
  arrange(-Count) %>%
  slice_max(Count, n = 10) %>%
  mutate(GENE = gsub('_conferring_resistance_.*', '', GENE)) %>%
  ggplot(aes(x = reorder(GENE, Count), y = Count)) +
  geom_bar(stat = 'identity', colour = 'black', fill = 'steelblue') +
  theme(panel.background = element_rect(fill = NA, colour = 'black'),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = 'bold', size = 9),
        axis.text = element_text(face = 'bold'),
        axis.title = element_text(face = 'bold')) +
  labs(x = 'Gene', y = 'Instances of Sharing') +
  ggtitle('Sharing of AMR-Associated Genes Between Mother and Infant')
pdf('sharing_card.pdf', width = 5, height = 3.5)
p
dev.off()

p <- annot %>%
  filter(DATABASE == 'vfdb') %>%
  group_by(GENE) %>%
  summarise(Count = length(GENE)) %>%
  arrange(-Count) %>%
  slice_max(Count, n = 10) %>%
  ggplot(aes(x = reorder(GENE, Count), y = Count)) +
  geom_bar(stat = 'identity', colour = 'black', fill = 'steelblue') +
  theme(panel.background = element_rect(fill = NA, colour = 'black'),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = 'bold', size = 9),
        axis.text = element_text(face = 'bold'),
        axis.title = element_text(face = 'bold')) +
  labs(x = 'Gene', y = 'Instances of Sharing') +
  ggtitle('Sharing of Virulence-Associated Genes Between Mother and Infant')
pdf('sharing_vfdb.pdf', width = 10, height = 3.5)
p
dev.off()

p <- annot %>%
  filter(DATABASE == 'plasmidfinder') %>%
  group_by(GENE) %>%
  summarise(Count = length(GENE)) %>%
  arrange(-Count) %>%
  ggplot(aes(x = reorder(GENE, Count), y = Count)) +
  geom_bar(stat = 'identity', colour = 'black', fill = 'steelblue') +
  theme(panel.background = element_rect(fill = NA, colour = 'black'),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = 'bold', size = 9),
        axis.text = element_text(face = 'bold'),
        axis.title = element_text(face = 'bold')) +
  labs(x = 'Gene', y = 'Instances of Sharing') +
  ggtitle('Sharing of Plasmid-Associated Genes Between Mother and Infant')
pdf('sharing_plasmidfinder.pdf', width = 10, height = 3.5)
p
dev.off()

p_abr <- annot %>%
  group_by(GENE, DATABASE) %>%
  summarise(Count = length(GENE)) %>%
  arrange(-Count) %>%
  group_by(DATABASE) %>%
  slice_max(Count, n = 10) %>%
  mutate(GENE = gsub('_conferring_resistance_.*', '', GENE)) %>%
  mutate(GENE = gsub('Escherichia_coli', 'E.coli', GENE)) %>%
  mutate(GENE = gsub('Bifidobacteria_intrinsic', 'Bif_intrinsic', GENE)) %>%
  ggplot(aes(x = reorder(GENE, Count), y = Count)) +
  geom_bar(stat = 'identity', colour = 'black', fill = 'steelblue') +
  ggforce::facet_col(~DATABASE, space = 'free', scales = 'free') +
  theme(panel.background = element_rect(fill = NA, colour = 'black'),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = 'bold', size = 9),
        axis.title = element_text(face = 'bold'),
        strip.text = element_text(face = 'bold')) +
  labs(x = 'Gene', y = '') 
  
pdf('sharing_abricate.pdf', width = 7, height = 8)
p_abr
dev.off()

##### MAKING A LIST OF ALL INFANT CONTIGS THAT ARE INHERITED FROM THE MOTHER
inherited <- events %>%
  select(Infant_Contig, TimeSource_Pair) %>%
  mutate(Infant_Contig = str_replace(Infant_Contig, '_length_.*', '')) %>%
  mutate(Infant_Contig = str_replace(Infant_Contig, '_NODE', '__NODE')) %>%
  rename(SampleContig = Infant_Contig) %>%
  mutate(Sample = str_replace(SampleContig, '__NODE.*', '')) %>%
  relocate(Sample, .before = SampleContig) 

inherited %>%
  write.table('inheritedcontigs.tsv', quote = F, row.names = F, sep = '\t')

##### ASSIGNING TAXONOMY TO INHERITED CONTIGS
inherited <- read.delim('inheritedcontigs.tsv', header = T, check.names = F, sep = '\t')

infantsamplelist <- list.files(path = 'kraken_outputs/', pattern = 'PC.*_output_named_formatted.txt', full.names = F) %>%
  str_replace(., '_output_named_formatted.txt', '')

infantcontig_taxonomy <- data.frame()
for (i in infantsamplelist){
  df <- read.delim(paste0('kraken_outputs/', i, '_output_named_formatted.txt'), header = F, sep = '\t') %>%
    select(V2, V6) %>%
    mutate(V2 = str_replace(V2, '_length.*', '')) %>%
    mutate(V2 = paste0(i, '__', V2)) %>%
    rename(Taxonomy = V6)
  
  infantcontig_taxonomy <- rbind(infantcontig_taxonomy, df)
}

taxonomicRanks <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

inherited_taxonomy <- inherited %>%
  left_join(infantcontig_taxonomy, by = c('SampleContig' = 'V2')) %>%
  separate(Taxonomy, sep = ';',
           into = taxonomicRanks)

inherited_taxonomy %>%
  write.table('taxonomy_inherited.tsv', quote = F, row.names = F, sep = '\t')

# MAKING A TABLE DESCRIBING THE NUMBER OF TRANSMISSION EVENTS PER SAMPLE
inherited_taxonomy <- read.delim('taxonomy_inherited.tsv', header = T, check.names = F)

inherited_taxonomy %>%
  group_by(Sample, TimeSource_Pair) %>%
  summarise(Count = n()) %>%
  write.table('sharingevents_persample.tsv', quote = F, row.names = F, sep = '\t')

inherited_taxonomy %>%
  group_by(Sample, TimeSource_Pair, Species) %>%
  summarise(Count = n()) %>%
  write.table('sharingevents_speciespersample.tsv', quote = F, row.names = F, sep = '\t')
  
inherited_taxonomy %>%
  group_by(Sample, TimeSource_Pair, Genus) %>%
  summarise(Count = n()) %>%
  write.table('sharingevents_genuspersample.tsv', quote = F, row.names = F, sep = '\t')

inherited_taxonomy %>%
  group_by(Sample, TimeSource_Pair, Family) %>%
  summarise(Count = n()) %>%
  write.table('sharingevents_familypersample.tsv', quote = F, row.names = F, sep = '\t')

# showing the distribution of the number of sample inherited by sample source and timepoint
inherited_taxonomy <- read.delim('taxonomy_inherited.tsv', header = T, check.names = F)
summarised <- inherited_taxonomy %>%
  select(Sample, TimeSource_Pair, Species) %>%
  distinct() %>%
  filter(Species != '') %>%
  group_by(Sample, TimeSource_Pair) %>%
  summarise(SpeciesShared = n()) %>%
  separate(TimeSource_Pair, into = c('Infant_Sample', 'Maternal_Sample'), sep = '__') %>%
  mutate(Maternal_Sample = gsub('^PB', 'Faeces ', Maternal_Sample)) %>%
  mutate(Maternal_Sample = gsub('^PM', 'Breast Milk ', Maternal_Sample)) %>%
  mutate(Maternal_Sample = gsub('^PO', 'Oral Rinse ', Maternal_Sample)) %>%
  mutate(Maternal_Sample = gsub('^PV', 'Vaginal Swab ', Maternal_Sample)) %>%
  mutate(Infant_Sample = gsub('^PC', 'Infant Faeces ', Infant_Sample)) %>%
  mutate(Maternal_Sample = gsub(' E$', ' Early', Maternal_Sample)) %>%
  mutate(Maternal_Sample = gsub(' L$', ' Late', Maternal_Sample)) %>%
  mutate(Maternal_Sample = gsub(' 1M$', ' 1 Month', Maternal_Sample)) %>%
  mutate(Infant_Sample = gsub(' DEL$', ' Delivery', Infant_Sample)) %>%
  mutate(Infant_Sample = gsub(' 1W$', ' 1 Week', Infant_Sample)) %>%
  mutate(Infant_Sample = gsub(' 1M$', ' 1 Month', Infant_Sample)) %>%
  mutate(Infant_Sample = factor(Infant_Sample, levels = c('Infant Faeces Delivery',
                                                          'Infant Faeces 1 Week',
                                                          'Infant Faeces 1 Month'))) %>%
  mutate(Maternal_Sample = factor(Maternal_Sample, levels = c('Faeces Early',
                                                              'Faeces Late',
                                                              'Faeces 1 Month',
                                                              'Breast Milk 1 Month',
                                                              'Vaginal Swab Early',
                                                              'Vaginal Swab Late',
                                                              'Oral Rinse Early')))

p <- summarised %>%
  group_by(Infant_Sample, Maternal_Sample) %>%
  summarise(Mean = mean(SpeciesShared)) %>%
  ggplot(aes(x = forcats::fct_rev(Maternal_Sample), y = Mean)) +
  geom_bar(aes(fill = Maternal_Sample), stat = 'identity', colour = 'black', show.legend = F) +
  coord_flip() +
  ggforce::facet_col(~Infant_Sample, space = 'free') +
  scale_fill_manual(values = sample_colours) +
  labs(x = 'Maternal Sample',
       y = 'Mean Number of Species') +
  theme(panel.background = element_rect(fill = NA, colour = 'black'),
        strip.text = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold'),
        panel.grid = element_blank()) +
  ggtitle('Species Exhibiting Strain-Sharing Events Per Dyad',
          subtitle = 'As determined by contig-sharing')

pdf('sharing_species.pdf', width = 6, height = 5)
p
dev.off()
svg('sharing_species.svg', width = 6, height = 5)
p
dev.off()

# where are the shared contigs coming from
inherited_taxonomy <- read.delim('taxonomy_inherited.tsv', header = T, check.names = F)

summarised <- inherited_taxonomy %>%
  group_by(TimeSource_Pair) %>%
  summarise(Count = n()) %>%
  separate(TimeSource_Pair, into = c('InfantTimeSource', 'MaternalTimeSource'), sep = '__') %>%
  mutate(Perc= 100*(Count/sum(Count)))

summarised %>%
  group_by(str_sub(MaternalTimeSource, 1, 2)) %>%
  summarise(Count = sum(Count)) %>%
  ungroup() %>%
  mutate(Perc = 100*(Count/sum(Count)))

summarised %>%
  group_by(InfantTimeSource) %>%
  summarise(Count = sum(Count)) %>%
  ungroup() %>%
  mutate(Perc = 100*(Count/sum(Count)))

summarised %>%
  group_by(InfantTimeSource, str_sub(MaternalTimeSource, 1, 2)) %>%
  summarise(Count = sum(Count)) %>%
  group_by(InfantTimeSource) %>%
  mutate(Perc = 100*(Count/sum(Count)))

summarised %>%
  group_by(InfantTimeSource) %>%
  pivot_wider(id_cols = MaternalTimeSource, names_from = InfantTimeSource, values_from = Count)
  
# SPECIES TRANSMISSION PREVALENCE
inherited_taxonomy %>%
  group_by(Sample, TimeSource_Pair, Species) %>%
  summarise(Count = n()) %>%
  mutate(Species = if_else(Species == '', 'Unclassified', Species)) %>%
  pivot_wider(id_cols = c('Sample', 'TimeSource_Pair'), names_from = Species, values_from = Count, values_fill = 0) %>%
  write.table('counts_inherited_species.tsv', quote = F, row.names = F, sep = '\t')

inherited_taxonomy %>%
  group_by(Sample, TimeSource_Pair, Genus) %>%
  summarise(Count = n()) %>%
  mutate(Genus = if_else(Genus == '', 'Unclassified', Genus)) %>%
  pivot_wider(id_cols = c('Sample', 'TimeSource_Pair'), names_from = Genus, values_from = Count, values_fill = 0) %>%
  write.table('counts_inherited_genus.tsv', quote = F, row.names = F, sep = '\t')

inherited_taxonomy %>%
  group_by(Sample, TimeSource_Pair, Family) %>%
  summarise(Count = n()) %>%
  mutate(Family = if_else(Family == '', 'Unclassified', Family)) %>%
  pivot_wider(id_cols = c('Sample', 'TimeSource_Pair'), names_from = Family, values_from = Count, values_fill = 0) %>%
  write.table('counts_inherited_family.tsv', quote = F, row.names = F, sep = '\t')

# SUMMARISING THE AMOUNT OF CONTIGS INHERITED AT EACH TAXONOMIC RANK BY TIMESOURCE_PAIR
inherited_taxonomy %>%
  pivot_longer(cols = taxonomicRanks, names_to = 'Rank', values_to = 'Taxon') %>% 
  group_by(TimeSource_Pair, Rank, Taxon) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = 100*(Count/sum(Count))) %>%
  arrange(TimeSource_Pair, Rank, -Percentage) %>%
  mutate(Taxon = if_else(Taxon == '', 'Unclassified', Taxon)) %>%
  write.table('taxonomy_inherited_timesourcepair_allranks.tsv', row.names = F, quote = F, sep = '\t')

inherited_taxonomy %>%
  group_by(TimeSource_Pair, Species) %>%
  summarise(Count = n()) %>% 
  mutate(Species = if_else(Species == '', 'Unclassified', Species)) %>%
  group_by(TimeSource_Pair) %>%
  mutate(Percentage = 100*(Count/sum(Count))) %>%
  arrange(TimeSource_Pair, -Percentage) %>%
  write.table('taxonomy_inherited_timesourcepair_species.tsv', row.names = F, quote = F, sep = '\t')

inherited_taxonomy %>%
  group_by(TimeSource_Pair, Genus) %>%
  summarise(Count = n()) %>% 
  mutate(Genus = if_else(Genus == '', 'Unclassified', Genus)) %>%
  group_by(TimeSource_Pair) %>%
  mutate(Percentage = 100*(Count/sum(Count))) %>%
  arrange(TimeSource_Pair, -Percentage) %>%
  write.table('taxonomy_inherited_timesourcepair_genus.tsv', row.names = F, quote = F, sep = '\t')

inherited_taxonomy %>%
  group_by(TimeSource_Pair, Family) %>%
  summarise(Count = n()) %>% 
  mutate(Family = if_else(Family == '', 'Unclassified', Family)) %>%
  group_by(TimeSource_Pair) %>%
  mutate(Percentage = 100*(Count/sum(Count))) %>%
  arrange(TimeSource_Pair, -Percentage) %>%
  write.table('taxonomy_inherited_timesourcepair_family.tsv', row.names = F, quote = F, sep = '\t')

# summarising the 10 most frequenctly inherited taxa (number of dyads) at each TimeSource_Pair
# at family, genus, and species level
inherited_taxonomy <- read.delim('taxonomy_inherited.tsv', header = T, check.names = F)

topSpecies <- data.frame()
for (i in unique(inherited_taxonomy$TimeSource_Pair)){
  df <- inherited_taxonomy %>%
    select(Sample, TimeSource_Pair, Species) %>% 
    filter(TimeSource_Pair == i & Species != '') %>%
    distinct() %>%
    group_by(Species) %>%
    summarise(Dyads = n()) %>%
    slice_max(Dyads, n = if_else(str_sub(i, 1, 5) == 'PCDEL',
                                 5, 10)) %>%
    mutate(TimeSource_Pair = i) %>%
    separate(TimeSource_Pair,
             into = c('Infant_Sample', 'Maternal_Sample'), 
             remove = F)
  
  topSpecies <- rbind(topSpecies, df)
}

topGenera <- data.frame()
for (i in unique(inherited_taxonomy$TimeSource_Pair)){
  df <- inherited_taxonomy %>%
    select(Sample, TimeSource_Pair, Genus) %>% 
    filter(TimeSource_Pair == i & Genus != '') %>%
    distinct() %>%
    group_by(Genus) %>%
    summarise(Dyads = n()) %>%
    slice_max(Dyads, n = if_else(str_sub(i, 1, 5) == 'PCDEL',
                                 5, 10)) %>%
    mutate(TimeSource_Pair = i) %>%
    separate(TimeSource_Pair,
             into = c('Infant_Sample', 'Maternal_Sample'), 
             remove = F)
  
  topGenera <- rbind(topGenera, df)
}

topFamilies <- data.frame()
for (i in unique(inherited_taxonomy$TimeSource_Pair)){
  df <- inherited_taxonomy %>%
    select(Sample, TimeSource_Pair, Family) %>% 
    filter(TimeSource_Pair == i & Family != '') %>%
    distinct() %>%
    group_by(Family) %>%
    summarise(Dyads = n()) %>%
    slice_max(Dyads, n = if_else(str_sub(i, 1, 5) == 'PCDEL',
                                 5, 10)) %>%
    mutate(TimeSource_Pair = i) %>%
    separate(TimeSource_Pair,
             into = c('Infant_Sample', 'Maternal_Sample'), 
             remove = F)
  
  topFamilies <- rbind(topFamilies, df)
}

topSpecies %>%
  write.table('taxonomy_inherited_topspecies.tsv', quote = F, row.names = F, sep = '\t')
topGenera %>%
  write.table('taxonomy_inherited_topgenera.tsv', quote = F, row.names = F, sep = '\t')
topFamilies %>%
  write.table('taxonomy_inherited_topfamilies.tsv', quote = F, row.names = F, sep = '\t')

# PLOTTING THE MOST INHERITED GENUS AT EACH INFANT TIMEPOINT
inherited_genus <- read.delim('taxonomy_inherited_topgenera.tsv', header = T, sep = '\t') %>%
  mutate(Maternal_Sample = gsub('^PB', 'Faeces ', Maternal_Sample)) %>%
  mutate(Maternal_Sample = gsub('^PM', 'Breast Milk ', Maternal_Sample)) %>%
  mutate(Maternal_Sample = gsub('^PO', 'Oral Rinse ', Maternal_Sample)) %>%
  mutate(Maternal_Sample = gsub('^PV', 'Vaginal Swab ', Maternal_Sample)) %>%
  mutate(Infant_Sample = gsub('^PC', 'Infant Faeces ', Infant_Sample)) %>%
  mutate(Maternal_Sample = gsub(' E$', ' Early', Maternal_Sample)) %>%
  mutate(Maternal_Sample = gsub(' L$', ' Late', Maternal_Sample)) %>%
  mutate(Maternal_Sample = gsub(' 1M$', ' 1 Month', Maternal_Sample)) %>%
  mutate(Infant_Sample = gsub(' DEL$', ' Delivery', Infant_Sample)) %>%
  mutate(Infant_Sample = gsub(' 1W$', ' 1 Week', Infant_Sample)) %>%
  mutate(Infant_Sample = gsub(' 1M$', ' 1 Month', Infant_Sample)) %>%
  mutate(Infant_Sample = factor(Infant_Sample, levels = c('Infant Faeces Delivery',
                                                          'Infant Faeces 1 Week',
                                                          'Infant Faeces 1 Month'))) %>%
  mutate(Maternal_Sample = factor(Maternal_Sample, levels = c('Faeces Early',
                                                              'Faeces Late',
                                                              'Faeces 1 Month',
                                                              'Breast Milk 1 Month',
                                                              'Vaginal Swab Early',
                                                              'Vaginal Swab Late',
                                                              'Oral Rinse Early')))

p_del <- inherited_genus %>%
  filter(Infant_Sample == 'Infant Faeces Delivery') %>%
  mutate(Genus = factor(Genus, levels = sort(unique(Genus), decreasing = T))) %>%
  ggplot(aes(x = Genus, y = Dyads)) +
  geom_bar(aes(fill = Maternal_Sample), colour = 'black', stat = 'identity', show.legend = F) +
  facet_wrap(~ Maternal_Sample, nrow = 1) +
  coord_flip() +
  ylim(0, 60) +
  scale_fill_manual(values = sample_colours) +
  theme_bw() +
  labs(x = '')
p_1w <- inherited_genus %>%
  filter(Infant_Sample == 'Infant Faeces 1 Week') %>%
  mutate(Genus = factor(Genus, levels = sort(unique(Genus), decreasing = T))) %>%
  ggplot(aes(x = Genus, y = Dyads)) +
  geom_bar(aes(fill = Maternal_Sample), colour = 'black', stat = 'identity', show.legend = F) +
  facet_wrap(~ Maternal_Sample, nrow = 1) +
  coord_flip() +
  ylim(0, 60) +
  scale_fill_manual(values = sample_colours) +
  theme_bw() +
  labs(x = '')
p_1m <- inherited_genus %>%
  filter(Infant_Sample == 'Infant Faeces 1 Month') %>%
  mutate(Genus = factor(Genus, levels = sort(unique(Genus), decreasing = T))) %>%
  ggplot(aes(x = Genus, y = Dyads)) +
  geom_bar(aes(fill = Maternal_Sample), colour = 'black', stat = 'identity', show.legend = F) +
  facet_wrap(~ Maternal_Sample, nrow = 1) +
  coord_flip() +
  ylim(0, 60) +
  scale_fill_manual(values = sample_colours) +
  theme_bw() +
  labs(x = '')

pdf('inherited_taxonomy_genus.pdf', width = 14, height = 15)
cowplot::plot_grid(plotlist = list(p_del, p_1w, p_1m),
                   ncol = 1,
                   rel_heights = c(0.5, 0.8, 1),
                   labels = c('Delivery', '1 Week', '1 Month'))
dev.off()

svg('inherited_taxonomy_genus.svg', width = 14, height = 15)
cowplot::plot_grid(plotlist = list(p_del, p_1w, p_1m),
                   ncol = 1,
                   rel_heights = c(0.5, 0.8, 1),
                   labels = c('Delivery', '1 Week', '1 Month'))
dev.off()

# adding taxonomy to the shared abricate output
shared_taxonomy <- read.delim('taxonomy_inherited.tsv', header = T, sep = '\t')
shared_abricate <- read.delim('sharing_abricate.tsv', header = T, sep = '\t')

shared_abrtax <- shared_abricate %>%
  select(SEQUENCE, GENE, DATABASE) %>%
  mutate(SEQUENCE = str_replace(SEQUENCE, '_length_.*', '')) %>%
  mutate(SEQUENCE = str_replace(SEQUENCE, '_NODE', '__NODE')) %>%
  left_join(shared_taxonomy %>%
              select(-Sample),
            by = c('SEQUENCE' = 'SampleContig')) %>%
  filter(str_sub(SEQUENCE, 7, 8) != '2W')

shared_abrtax %>%
  write.table('shared_abricate_taxonomy.tsv', row.names = F, quote = F, sep = '\t')

p <- shared_abrtax %>%
  group_by(TimeSource_Pair, GENE, DATABASE, Genus) %>%
  summarise(Count = n()) %>%
  arrange(-Count) %>%
  filter(DATABASE == 'card') %>%
  mutate(GENE = str_replace(GENE, '_conferring_.*', '')) %>%
  mutate(Genus = if_else(Genus == '', 'Unclassified', Genus)) %>%
  ggplot(aes(x = GENE, y = Count)) +
  geom_bar(aes(fill = Genus), colour = 'black', stat = 'identity') +
  facet_wrap(~TimeSource_Pair, nrow = 1) +
  theme_bw() +
  coord_flip() +
  scale_fill_manual(values = c("#b5aa44",
                               "#7065cc",
                               "#65b347",
                               "#bc50bb",
                               "#617b35",
                               "#cd447c",
                               "#4fb28e",
                               "#cf4938",
                               "#6092d1",
                               "#c47f3b",
                               "#bc7aba",
                               "#c96d6f",
                               "#6e6d69"))

pdf('shared_abricate_taxonomy_timesourcepairGenus.pdf', width = 20, height = 9)
p
dev.off()

p <- shared_abrtax %>%
  group_by(TimeSource_Pair, GENE, DATABASE, Genus) %>%
  summarise(Count = n()) %>%
  arrange(-Count) %>%
  filter(DATABASE == 'card') %>%
  mutate(GENE = str_replace(GENE, '_conferring_.*', '')) %>%
  mutate(Genus = if_else(Genus == '', 'Unclassified', Genus)) %>%
  filter(str_sub(TimeSource_Pair, 1, 4) == 'PC1M') %>%
  ggplot(aes(x = GENE, y = Count)) +
  geom_bar(aes(fill = Genus), colour = 'black', stat = 'identity') +
  facet_wrap(~TimeSource_Pair, nrow = 1) +
  theme_bw() +
  coord_flip() +
  scale_fill_manual(values = c("#b5aa44",
                               "#7065cc",
                               "#65b347",
                               "#bc50bb",
                               "#617b35",
                               "#cd447c",
                               "#4fb28e",
                               "#cf4938",
                               "#6092d1",
                               "#c47f3b",
                               "#bc7aba",
                               "#c96d6f",
                               "#6e6d69"))

pdf('shared_abricate_taxonomy_timesourcepairGenus_PC1M.pdf', width = 12, height = 9)
p
dev.off()

p <- shared_abrtax %>%
  group_by(TimeSource_Pair, GENE, DATABASE, Genus) %>%
  summarise(Count = n()) %>%
  arrange(-Count) %>%
  filter(DATABASE == 'card') %>%
  mutate(GENE = str_replace(GENE, '_conferring_.*', '')) %>%
  mutate(Genus = if_else(Genus == '', 'Unclassified', Genus)) %>%
  filter(str_sub(TimeSource_Pair, 1, 4) == 'PC1W') %>%
  ggplot(aes(x = GENE, y = Count)) +
  geom_bar(aes(fill = Genus), colour = 'black', stat = 'identity') +
  facet_wrap(~TimeSource_Pair, nrow = 1) +
  theme_bw() +
  coord_flip() +
  scale_fill_manual(values = c("#b5aa44",
                               "#7065cc",
                               "#65b347",
                               "#bc50bb",
                               "#617b35",
                               "#cd447c",
                               "#4fb28e",
                               "#cf4938",
                               "#6092d1",
                               "#c47f3b",
                               "#bc7aba",
                               "#c96d6f",
                               "#6e6d69"))

pdf('shared_abricate_taxonomy_timesourcepairGenus_PC1W.pdf', width = 12, height = 3)
p
dev.off()

p <- shared_abrtax %>%
  group_by(TimeSource_Pair, GENE, DATABASE, Genus) %>%
  summarise(Count = n()) %>%
  arrange(-Count) %>%
  filter(DATABASE == 'card') %>%
  mutate(GENE = str_replace(GENE, '_conferring_.*', '')) %>%
  mutate(Genus = if_else(Genus == '', 'Unclassified', Genus)) %>%
  filter(str_sub(TimeSource_Pair, 1, 5) == 'PCDEL') %>%
  ggplot(aes(x = GENE, y = Count)) +
  geom_bar(aes(fill = Genus), colour = 'black', stat = 'identity') +
  facet_wrap(~TimeSource_Pair, nrow = 1) +
  theme_bw() +
  coord_flip() +
  scale_fill_manual(values = c("#b5aa44",
                               "#7065cc",
                               "#65b347",
                               "#bc50bb",
                               "#617b35",
                               "#cd447c",
                               "#4fb28e",
                               "#cf4938",
                               "#6092d1",
                               "#c47f3b",
                               "#bc7aba",
                               "#c96d6f",
                               "#6e6d69"))

pdf('shared_abricate_taxonomy_timesourcepairGenus_PCDEL.pdf', width = 12, height = 4)
p
dev.off()

# looking for inherited carbohydrate-activate enzymes and assigning taxonomy
dbcan <- read.delim('dbcan_all_ungrouped.tsv', header = T, check.names = F, sep = '\t')
inherited <- read.delim('taxonomy_inherited.tsv', header = T, sep = '\t')

inherited_tax_dbcan <- inherited %>%
  left_join(dbcan, by = c('SampleContig' = 'Contig'))

inhertied_tax_GH <- inherited_tax_dbcan %>%
  filter(grepl("^GH", CAZY)) %>%
  separate(TimeSource_Pair,
           into = c('Infant_Sample', 'Maternal_Sample'), 
           remove = F)

inherited_speciesGH_summary <- inhertied_tax_GH %>%
  select(Sample, Infant_Sample, Maternal_Sample, Species, CAZY) %>%
  distinct() %>%
  group_by(Infant_Sample, Maternal_Sample, Species, CAZY) %>%
  summarise(Dyads = n()) %>%
  filter(Species != '') %>%
  arrange(-Dyads)
inherited_speciesGH_summary %>%
  write.table('inherited_speciesGH_summary.tsv', row.names = F, quote = F, sep = '\t')

inherited_speciesGH_summary %>%
  group_by(Infant_Sample, Maternal_Sample, Species) %>%
  summarise(NumberOfCazy = n()) %>%
  arrange(-NumberOfCazy)

inherited_speciesGH_summary %>%
  group_by(Infant_Sample, Maternal_Sample, CAZY) %>%
  summarise(NumberOfSpecies = n()) %>%
  arrange(-NumberOfSpecies)

