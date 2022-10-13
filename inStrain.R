library('tidyverse')

instrain <- read.delim('instrain_genomeWide_compare.tsv') %>%
  mutate(genome = str_replace(genome, '.fasta', '')) %>%
  mutate(genome = str_replace(genome, '.fa', '')) %>%
  mutate(name1 = str_replace(name1, '.sorted.bam', '')) %>%
  mutate(name2 = str_replace(name2, '.sorted.bam', ''))

# filtering dataset to only include strain-sharing comparisons based on instrain recommended cutoffs
instrain_samestrain <- instrain %>%
  filter(popANI >= 0.99999 & percent_compared >= 0.5)

instrain_samestrain %>%
  write.table('samestrain.tsv', quote = F, row.names = F, sep = '\t')

# importing taxonomy and sample information of isolates and MAGs
genomeTaxonomy <- read.delim('referencegenometaxonomy.txt', header = T, sep = '\t')

# adding sample and taxonomy information to strain-sharing events
instrain_samestrain <- read.delim('samestrain.tsv', header = T, sep = '\t') %>%
  left_join(genomeTaxonomy %>%
              select(-Sample), by = c('genome' = 'ReferenceGenome')) %>%
  left_join(genomeTaxonomy %>%
              select(-Taxonomy), by = c('name1' = 'ReferenceGenome')) %>%
  left_join(genomeTaxonomy %>%
              select(-Taxonomy), by = c('name2' = 'ReferenceGenome')) %>%
  rename(Sample1 = Sample.x, Sample2 = Sample.y) %>%
  mutate(Sample1 = if_else(str_sub(name1, 1, 1) == 'P', name1, Sample1)) %>%
  mutate(Sample2 = if_else(str_sub(name2, 1, 1) == 'P', name2, Sample2))

instrain_samestrain %>%
  write.table('samestrain_withmeta.tsv', quote = F, row.names = F, sep = '\t')

# filtering some comparisons
instrain_samestrain <- read.delim('samestrain_withmeta.tsv', header = T, sep = '\t') %>%
  filter(Sample1 != Sample2) %>% # removes comparisons within the same sample (eg. isolate vs. metagenome)
  filter(str_sub(Sample1, 3, 5) == str_sub(Sample2, 3, 5)) %>% # removes comparisons from different families
  filter((str_sub(Sample1, 1, 2) == 'PC' & str_sub(Sample2, 1, 2) %in% c('PB', 'PM', 'PO', 'PV')) |
           (str_sub(Sample2, 1, 2) == 'PC' & str_sub(Sample1, 1, 2) %in% c('PB', 'PM', 'PO', 'PV'))) %>% # removes comparisons that aren't between mother and infant
  mutate(SharingType = if_else((str_sub(Sample1, -2, -1) %in% c('_E', '_L') & !str_sub(Sample2, -2, -1) %in% c('_E', '_L')) | 
                                 (str_sub(Sample2, -2, -1) %in% c('_E', '_L') & !str_sub(Sample1, -2, -1) %in% c('_E', '_L')),
                               'Transmission', 'Sharing')) %>% # classifies SharingTypes based on whether the strain was detected in the mother before birth
  mutate(Family = paste0('PX', str_sub(Sample1, 3, 5)))

instrain_samestrain %>%
  write.table('samestrain_sharingevents.tsv', quote = F, row.names = F, sep = '\t')

# summarising sharing events 
instrain_samestrain <- read.delim('samestrain_sharingevents.tsv', header = T, sep = '\t')

instrain_samestrain %>%
  group_by(Family, genome, SharingType) %>%
  summarise(Events = n()) %>%
  pivot_wider(id_cols = c(Family, genome), names_from = SharingType, values_from = Events, values_fill = 0) %>%
  left_join(genomeTaxonomy %>%
              select(-Sample),
            by = c('genome' = 'ReferenceGenome')) %>%
  write.table('summary_sharingTypeGenomeFamily.tsv', row.names = F, quote = F, sep = '\t')

##### BUILDING A LM FOR CONANI VS POPANI
instrain <- read.delim('instrain_genomeWide_compare.tsv') %>%
  mutate(genome = str_replace(genome, '.fasta', '')) %>%
  mutate(genome = str_replace(genome, '.fa', '')) %>%
  mutate(name1 = str_replace(name1, '.sorted.bam', '')) %>%
  mutate(name2 = str_replace(name2, '.sorted.bam', ''))

model <- lm(conANI ~ popANI, data = instrain)
summary(model)
predict(model, data.frame(popANI = c(0.99999)))
coef(model)

ggplot(instrain, aes(x = popANI, y = conANI)) +
  geom_point(shape = 21, fill = 'black', alpha = 0.5) +
  stat_smooth(method = 'lm', col = 'red')

# looking at Strain Clusters
clusters_in <- read.delim('instrain_strain_clusters.tsv', header = T, sep = '\t') %>%
  mutate(sample = str_replace(sample, '\\.sorted\\.bam', '')) %>%
  mutate(genome = str_replace(genome, '\\.fa.*', '')) 

clusters_in %>%
  group_by(cluster, genome) %>%
  summarise(Samples = paste(sample, collapse = ';'),
            ClusterSize = n()) %>%
  filter(ClusterSize > 1) %>%
  arrange(-ClusterSize) %>%
  left_join(read.delim('referencegenometaxonomy.txt', header = T, sep = '\t') %>%
              select(ReferenceGenome, Taxonomy),
            by = c('genome' = 'ReferenceGenome')) %>%
  relocate(c(ClusterSize, Taxonomy), .before = Samples) %>%
  write.table('strain_clusters.tsv', quote = F, row.names = F, sep = '\t')

clusters_in2 <- clusters_in %>%
  left_join(read.delim('mm001-484_strainlist.tsv', header = T, sep = '\t') %>%
              select(Strain, Sample),
            by = c('sample' = 'Strain')) %>%
  mutate(Sample = if_else(is.na(Sample), sample, Sample))

clusters_in2 %>%
  group_by(cluster, genome) %>%
  summarise(Samples = paste(Sample, collapse = ';'),
            ClusterSize = n()) %>%
  filter(ClusterSize > 1) %>%
  arrange(-ClusterSize) %>%
  left_join(read.delim('referencegenometaxonomy.txt', header = T, sep = '\t') %>%
              select(ReferenceGenome, Taxonomy),
            by = c('genome' = 'ReferenceGenome')) %>%
  relocate(c(ClusterSize, Taxonomy), .before = Samples) %>%
  write.table('strain_clusters_genomesrenamed.tsv', quote = F, row.names = F, sep = '\t')

# determining if each cluster represents a single family or multiple families
cluster_familycount <- clusters_in2 %>%
  mutate(Family = str_sub(Sample, 3, 5)) %>%
  select(cluster, Family) %>%
  distinct() %>%
  group_by(cluster) %>%
  summarise(FamilyCount = n()) 
clusters_in %>%
  group_by(cluster, genome) %>%
  summarise(Samples = paste(sample, collapse = ';'),
            ClusterSize = n()) %>%
  filter(ClusterSize > 1) %>%
  arrange(-ClusterSize) %>%
  left_join(read.delim('referencegenometaxonomy.txt', header = T, sep = '\t') %>%
              select(ReferenceGenome, Taxonomy),
            by = c('genome' = 'ReferenceGenome')) %>%
  relocate(c(ClusterSize, Taxonomy), .before = Samples) %>%
  left_join(cluster_familycount) %>%
  relocate(FamilyCount, .before = Taxonomy) %>%
  write.table('strain_clusters.tsv', quote = F, row.names = F, sep = '\t')
clusters_in2 %>%
  group_by(cluster, genome) %>%
  summarise(Samples = paste(Sample, collapse = ';'),
            ClusterSize = n()) %>%
  filter(ClusterSize > 1) %>%
  arrange(-ClusterSize) %>%
  left_join(read.delim('referencegenometaxonomy.txt', header = T, sep = '\t') %>%
              select(ReferenceGenome, Taxonomy),
            by = c('genome' = 'ReferenceGenome')) %>%
  relocate(c(ClusterSize, Taxonomy), .before = Samples) %>%
  left_join(cluster_familycount) %>%
  relocate(FamilyCount, .before = Taxonomy) %>%
  write.table('strain_clusters_genomesrenamed.tsv', quote = F, row.names = F, sep = '\t')

# making a contingency table for each species - number of single/multiple family clusters
singlefam_cont <- clusters_in %>%
  group_by(cluster, genome) %>%
  summarise(Samples = paste(sample, collapse = ';'),
            ClusterSize = n()) %>%
  filter(ClusterSize > 1) %>%
  left_join(read.delim('referencegenometaxonomy.txt', header = T, sep = '\t') %>%
              select(ReferenceGenome, Taxonomy),
            by = c('genome' = 'ReferenceGenome')) %>%
  left_join(cluster_familycount) %>%
  mutate(SingleFamily = if_else(FamilyCount == 1, 'SingleFamily', 'MultipleFamilies')) %>%
  group_by(Taxonomy, SingleFamily) %>%
  summarise(Count = n()) %>%
  pivot_wider(id_cols = Taxonomy, names_from = SingleFamily, values_from = Count, values_fill = 0) %>%
  mutate(Ratio = SingleFamily/MultipleFamilies) %>%
  arrange(-Ratio, -SingleFamily)

singlefam_cont %>%
  write.table('strain_clusters_singlefamily.tsv', row.names = F, quote = F, sep = '\t')

p <- singlefam_cont %>%
  mutate(Sum = SingleFamily + MultipleFamilies) %>%
  mutate(Species = factor(Taxonomy)) %>%
  ungroup() %>%
  select(Species, SingleFamily, MultipleFamilies, Sum) %>%
  pivot_longer(cols = -c(Species, Sum), names_to = 'Families', values_to = 'Count') %>%
  mutate(StrainClusters = if_else(Families == 'SingleFamily', Count, -Count)) %>%
  ggplot(aes(x = reorder(Species, Sum), y = StrainClusters)) +
  geom_bar(aes(fill = Families), stat = 'identity', position = 'identity') +
  coord_flip() +
  theme_bw() +
  theme(legend.position = c(0.75, 0.9)) +
  labs(x = 'Species', y = 'Strain Clusters')

pdf('strain_clusters_singlefamily.pdf', width = 6, height = 25)
p  
dev.off()

p <- singlefam_cont %>%
  mutate(Sum = SingleFamily + MultipleFamilies) %>%
  mutate(Species = factor(Taxonomy)) %>%
  ungroup() %>%
  select(Species, SingleFamily, MultipleFamilies, Sum) %>%
  slice_max(Sum, n = 15) %>%
  pivot_longer(cols = -c(Species, Sum), names_to = 'Families', values_to = 'Count') %>%
  mutate(StrainClusters = if_else(Families == 'SingleFamily', Count, -Count)) %>%
  ggplot(aes(x = reorder(Species, Sum), y = StrainClusters)) +
  geom_bar(aes(fill = Families), stat = 'identity', position = 'identity') +
  coord_flip() +
  theme_bw() +
  theme(legend.position = c(0.7, 0.15)) +
  labs(x = 'Species', y = 'Strain Clusters')

pdf('strain_clusters_singlefamily_top15.pdf', width = 5, height = 4)
p  
dev.off()

p <- singlefam_cont %>%
  filter(MultipleFamilies != 0) %>%
  mutate(Sum = SingleFamily + MultipleFamilies) %>%
  mutate(Species = factor(Taxonomy)) %>%
  ungroup() %>%
  select(Species, SingleFamily, MultipleFamilies, Sum) %>%
  pivot_longer(cols = -c(Species, Sum), names_to = 'Families', values_to = 'Count') %>%
  mutate(StrainClusters = if_else(Families == 'SingleFamily', Count, -Count)) %>%
  ggplot(aes(x = reorder(Species, Sum), y = StrainClusters)) +
  geom_bar(aes(fill = Families), stat = 'identity', position = 'identity') +
  coord_flip() +
  theme_bw() +
  theme(legend.position = c(0.65, 0.2)) +
  labs(x = 'Species', y = 'Strain Clusters')

pdf('strain_clusters_singlefamily_multipleonly.pdf', width = 5, height = 4)
p  
dev.off()

# looking at Strain Clusters identified by metagenomic data only
clusters_in_meta <- read.delim('instrain_strain_clusters.tsv', header = T, sep = '\t') %>%
  mutate(sample = str_replace(sample, '\\.sorted\\.bam', '')) %>%
  mutate(genome = str_replace(genome, '\\.fa.*', '')) %>%
  filter(!str_sub(sample, 1, 2) %in% c('MB', 'MM'))

cluster_familycount_meta <- clusters_in_meta %>%
  mutate(Family = str_sub(sample, 3, 5)) %>%
  select(cluster, Family) %>%
  distinct() %>%
  group_by(cluster) %>%
  summarise(FamilyCount = n()) 

clusters_in_meta %>%
  group_by(cluster, genome) %>%
  summarise(Samples = paste(sample, collapse = ';'),
            ClusterSize = n()) %>%
  filter(ClusterSize > 1) %>%
  arrange(-ClusterSize) %>%
  left_join(read.delim('referencegenometaxonomy.txt', header = T, sep = '\t') %>%
              select(ReferenceGenome, Taxonomy),
            by = c('genome' = 'ReferenceGenome')) %>%
  relocate(c(ClusterSize, Taxonomy), .before = Samples) %>%
  left_join(cluster_familycount_meta) %>%
  relocate(FamilyCount, .before = Taxonomy) %>%
  write.table('strain_clusters_metagenomiconly.tsv', quote = F, row.names = F, sep = '\t')

singlefam_cont_meta <- clusters_in_meta %>%
  group_by(cluster, genome) %>%
  summarise(Samples = paste(sample, collapse = ';'),
            ClusterSize = n()) %>%
  filter(ClusterSize > 1) %>%
  left_join(read.delim('referencegenometaxonomy.txt', header = T, sep = '\t') %>%
              select(ReferenceGenome, Taxonomy),
            by = c('genome' = 'ReferenceGenome')) %>%
  left_join(cluster_familycount_meta) %>%
  mutate(SingleFamily = if_else(FamilyCount == 1, 'SingleFamily', 'MultipleFamilies')) %>%
  group_by(Taxonomy, SingleFamily) %>%
  summarise(Count = n()) %>%
  pivot_wider(id_cols = Taxonomy, names_from = SingleFamily, values_from = Count, values_fill = 0) %>%
  mutate(Ratio = SingleFamily/MultipleFamilies) %>%
  arrange(-Ratio, -SingleFamily)

singlefam_cont_meta %>%
  write.table('strain_clusters_singlefamily_metagenomiconly.tsv', row.names = F, quote = F, sep = '\t')

p <- singlefam_cont_meta %>%
  mutate(Sum = SingleFamily + MultipleFamilies) %>%
  mutate(Species = factor(Taxonomy)) %>%
  ungroup() %>%
  select(Species, SingleFamily, MultipleFamilies, Sum) %>%
  pivot_longer(cols = -c(Species, Sum), names_to = 'Families', values_to = 'Count') %>%
  mutate(StrainClusters = if_else(Families == 'SingleFamily', Count, -Count)) %>%
  ggplot(aes(x = reorder(Species, Sum), y = StrainClusters)) +
  geom_bar(aes(fill = Families), stat = 'identity', position = 'identity') +
  coord_flip() +
  theme_bw() +
  theme(legend.position = c(0.75, 0.9)) +
  labs(x = 'Species', y = 'Strain Clusters')

pdf('strain_clusters_singlefamily_metagenomiconly.pdf', width = 6, height = 25)
p  
dev.off()

p <- singlefam_cont_meta %>%
  mutate(Sum = SingleFamily + MultipleFamilies) %>%
  mutate(Species = factor(Taxonomy)) %>%
  ungroup() %>%
  select(Species, SingleFamily, MultipleFamilies, Sum) %>%
  slice_max(Sum, n = 15) %>%
  pivot_longer(cols = -c(Species, Sum), names_to = 'Families', values_to = 'Count') %>%
  mutate(StrainClusters = if_else(Families == 'SingleFamily', Count, -Count)) %>%
  ggplot(aes(x = reorder(Species, Sum), y = StrainClusters)) +
  geom_bar(aes(fill = Families), stat = 'identity', position = 'identity') +
  coord_flip() +
  theme_bw() +
  theme(legend.position = c(0.7, 0.15)) +
  labs(x = 'Species', y = 'Strain Clusters')

pdf('strain_clusters_singlefamily_top15_metagenomiconly.pdf', width = 5, height = 4)
p  
dev.off()

p <- singlefam_cont_meta %>%
  filter(MultipleFamilies != 0) %>%
  mutate(Sum = SingleFamily + MultipleFamilies) %>%
  mutate(Species = factor(Taxonomy)) %>%
  ungroup() %>%
  select(Species, SingleFamily, MultipleFamilies, Sum) %>%
  pivot_longer(cols = -c(Species, Sum), names_to = 'Families', values_to = 'Count') %>%
  mutate(StrainClusters = if_else(Families == 'SingleFamily', Count, -Count)) %>%
  ggplot(aes(x = reorder(Species, Sum), y = StrainClusters)) +
  geom_bar(aes(fill = Families), stat = 'identity', position = 'identity') +
  coord_flip() +
  theme_bw() +
  theme(legend.position = c(0.6, 0.2)) +
  labs(x = 'Species', y = 'Strain Clusters')

pdf('strain_clusters_singlefamily_multipleonly_metagenomiconly.pdf', width = 5, height = 4)
p  
dev.off()

# determining how much the isolate data added to strain cluster analysis
groupColours <- c('Metagenomic' = 'blue',
                  'Isolate' = 'red',
                  'Combined' = 'purple')
clus_famcount <- clusters_in %>%
  group_by(cluster, genome) %>%
  summarise(Samples = paste(sample, collapse = ';'),
            ClusterSize = n()) %>%
  filter(ClusterSize > 1) %>%
  left_join(read.delim('referencegenometaxonomy.txt', header = T, sep = '\t') %>%
              select(ReferenceGenome, Taxonomy),
            by = c('genome' = 'ReferenceGenome')) %>%
  left_join(cluster_familycount) %>%
  mutate(SingleFamily = if_else(FamilyCount == 1, 'SingleFamily', 'MultipleFamilies')) %>%
  select(cluster, Taxonomy, SingleFamily)
clus_famcount_meta <- clusters_in_meta %>%
  group_by(cluster, genome) %>%
  summarise(Samples = paste(sample, collapse = ';'),
            ClusterSize = n()) %>%
  filter(ClusterSize > 1) %>%
  left_join(read.delim('referencegenometaxonomy.txt', header = T, sep = '\t') %>%
              select(ReferenceGenome, Taxonomy),
            by = c('genome' = 'ReferenceGenome')) %>%
  left_join(cluster_familycount) %>%
  mutate(SingleFamily = if_else(FamilyCount == 1, 'SingleFamily', 'MultipleFamilies')) %>%
  select(cluster, SingleFamily)

df <- clus_famcount %>%
  left_join(clus_famcount_meta, by = 'cluster') %>%
  rename('Combined' = 'SingleFamily.x',
         'Meta' = 'SingleFamily.y') %>%
  mutate(Meta = replace_na(Meta, 'Undetected')) 

clus_famcount_summary <- df %>%
  group_by(Taxonomy, Combined, Meta) %>%
  summarise(StrainClusterCount = n())

# how many different b. breve strains were detected in the dataset?
clusters_in <- read.delim('instrain_strain_clusters.tsv', header = T, sep = '\t') %>%
  mutate(sample = str_replace(sample, '\\.sorted\\.bam', '')) %>%
  mutate(genome = str_replace(genome, '\\.fa.*', '')) %>%
  group_by(cluster, genome) %>%
  summarise(Samples = paste(sample, collapse = ';'),
            ClusterSize = n()) %>%
  arrange(-ClusterSize) %>%
  left_join(read.delim('referencegenometaxonomy.txt', header = T, sep = '\t') %>%
              select(ReferenceGenome, Taxonomy),
            by = c('genome' = 'ReferenceGenome')) %>%
  relocate(c(ClusterSize, Taxonomy), .before = Samples)
  
breve_clusters <- clusters_in %>%
  filter(Taxonomy == 'Bifidobacterium breve')

