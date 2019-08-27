# Load dependencies
deps <- c('vegan', 'shape', 'plotrix', 'reshape2', 'GMD', 'randomForest', 'AUCRF','RColorBrewer', 'gplots','viridis', 'scales','Hmisc','VennDiagram','dtw')
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(deps)

installed.packages('ggplot2')

library(ggplot2)
library(ggdendro)
library(reshape2)
library(grid)
# library(ggplot2)
library(tidyr)
library(AUCRF)
library(randomForest)
library(dplyr)
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(gridExtra)
library(scales)
default_theme <- theme_get()

#### Read in data ####
data_in <- read.table('../Lactobacillus/Data/pFBA_likelihoods_output_3.tsv', sep='\t', header=TRUE, colClasses = "character", row.names=1)

data <- sapply(data_in, as.numeric)
rownames(data) <- rownames(data_in)
rm(data_in)
data <- as.data.frame(data)

## Human Gut, Vaginal, Oral
meta <- read.csv('../Lactobacillus/Data/All_meta_data.csv', sep=",", header=TRUE, colClasses = "character", row.names=1)
meta <- meta[,c("Fermentation","Isolation.Site","species")]

target <- c("enviro","food_alcohol","other","unknown","waste","gut_other", "gut_swine", "gut_rodent","food_wheat", "food_milk", "food_veg","",
            "food_meat","food_other","gut_avian","gut_mammal", "blood_human", "bv_human", "other_human", "stomach_human", "milk_human", "gut_lab")
meta_human <- meta
meta_human$Row.names <- rownames(meta_human)
meta_human <- dplyr::filter(meta_human, !(Isolation.Site %in% target))

target <- c("adolescentis","animalis","bifidum","breve","longum","agilis","parabuchneri","ruminis","mucosae")
meta_human <- dplyr::filter(meta_human, !(species %in% target))
meta_human_global <- meta_human

data <- data[c(meta_human$Row.names), ]

data <- data[complete.cases(data), ]

data_global <- data

#### Metabolomics Ordination Analysis ####
data_dist <- vegdist(data_global, method='bray')
data_nmds <- as.data.frame(metaMDS(data_dist, k=2, trymax=100)$points)

# which(data_nmds$MDS2 > 0.1)

## Add Meta-data, then plot
rownames(meta_human) <- meta_human$Row.names
meta_human$Row.names <-  NULL

data_meta <- merge(data, meta_human, by = "row.names")

data_nmds_meta <- merge(data_nmds, meta, by = "row.names")
data_nmds_meta$Fermentation <- as.factor(data_nmds_meta$Fermentation)
data_nmds_meta$Isolation.Site <- as.factor(data_nmds_meta$Isolation.Site)

tiff("pFBA_figures/human_ord_iso_site.tiff", width = 6.4, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
ggplot(data_nmds_meta, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=Isolation.Site), alpha = 0.8) +
  labs(color = "Isolation Site", x='NMDS axis 1', y='NMDS axis 2') +
  scale_color_discrete(labels=c("gut","oral","vaginal"))
dev.off()

# spec_names = c("acidophilus","brevis","casei","crispatus","delbrueckii","fermentum","gasseri","helveticus",
#                "iners","jensenii","johnsonii","paracasei","plantarum","reuteri","rhamnosus","salivarius")

palette <- brewer.pal(8,"Dark2")
palette2 <- rep(palette,each=2)

tiff("pFBA_figures/human_ord_species.tiff", width = 6.4, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
ggplot(data_nmds_meta, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=species,shape=species), alpha = 0.8) +
  labs(color = "Species", x='NMDS axis 1', y='NMDS axis 2') +
  scale_colour_manual(name = "Species",
                      values = palette2)+
  scale_shape_manual(name = "Species",
                     values = c(16,15,16,15,16,15,16,15,16,15,16,15,16,15,16,15))
dev.off()

#### Metabolomics Ordination Analysis Normalized Data ####
meta_human <- meta_human_global
data <- data_global
data <- data %>% select(-contains("1"))
data <- data %>% select(-contains("2"))
data$Row.names <- rownames(data)
data_meta <- merge(data, meta_human, by = "Row.names")
rownames(data_meta) <- data_meta$Row.names
data_meta_slim <- data_meta
data_meta_slim$Row.names <- NULL
data_meta_slim$Fermentation <- NULL
data_meta_slim$Isolation.Site <- NULL
data_meta_slim$species <- NULL
# View(iris)

# data_meta_alt <- data_meta_slim
data_meta_slim_norm <- as.data.frame(t(apply(data_meta_slim, 1, function(x) ((x)/sd(x)) )))
# rownames(data_meta_slim_norm) <- rownames(data_meta_slim)
data2 <- data_meta_slim_norm
data_dist <- vegdist(data2, method='bray')
# data_nmds <- as.data.frame(metaMDS(data_dist, k=2, trymax=100)$points)
data_nmds_df <- as.data.frame(cmdscale(data_dist, k=2, eig=TRUE)$points)
data_nmds <- cmdscale(data_dist, k=2, eig=TRUE)
# rownames(data_nmds) <- rownames(data2)
colnames(data_nmds_df) <- c('MDS1','MDS2')

# barplot(data_nmds$eig[1:10])
# data_nmds$eig[1]/sum(data_nmds$eig[1:10])
# data_nmds$eig[2]/sum(data_nmds$eig[1:10])

# data_pcoa <- pcoa(data_dist, correction = "cailliez")
# barplot(data_pcoa$values$Relative_eig[1:10])
# barplot(data_pcoa$values$Eigenvalues[1:10])
# biplot.pcoa(data_pcoa, data2)

# test_dist <- metaMDS(data_dist, k=2, trymax=100)

# biplot(data_nmds)
# pcoa.plot(data_nmds)
# ggord(data_nmds, as.factor(meta_human$Isolation.Site))

# which(data_nmds$MDS2 > 0.1)

## Add Meta-data, then plot
rownames(meta_human) <- meta_human$Row.names
meta_human$Row.names <-  NULL

# data_meta <- merge(data, meta_human, by = "row.names")

data_nmds_meta <- merge(data_nmds_df, meta_human, by = "row.names")
data_nmds_meta$Fermentation <- as.factor(data_nmds_meta$Fermentation)
data_nmds_meta$Isolation.Site <- as.factor(data_nmds_meta$Isolation.Site)

data_nmds_meta_means <- merge(data_nmds_meta,aggregate(cbind(mean.x=MDS1,mean.y=MDS2)~Isolation.Site,data_nmds_meta,mean),by="Isolation.Site")

contri_1 = as.numeric(round(100*data_nmds$eig[1]/sum(data_nmds$eig[1:20])))
contri_2 = as.numeric(round(100*data_nmds$eig[2]/sum(data_nmds$eig[1:20])))

tiff("pFBA_figures/human_ord_iso_site_norm.tiff", width = 6.4, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
p1 = ggplot(data_nmds_meta_means, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=factor(Isolation.Site)), alpha = 0.6, size = 2) +
  labs(color = "Isolation Site", 
       x=paste(sprintf('MDS axis 1 (%i',contri_1),'%)'), 
       y=paste(sprintf('MDS axis 2 (%i',contri_2),'%)')) +
  scale_color_discrete(labels=c("Intestinal","Oral","Vaginal")) + 
  #geom_point(aes(x=mean.x,y=mean.y,color=factor(Isolation.Site)),size=4) +
  geom_segment(aes(x=mean.x, y=mean.y, xend=MDS1, yend=MDS2,color=factor(Isolation.Site)), size = 1, alpha = 0.2) +
  theme(text = element_text(size=16))
dev.off()

## Adonis test for isolation site
sample_meta <- data_nmds_meta
rownames(sample_meta) <- sample_meta$Row.names
sample_meta <- sample_meta[,c("Isolation.Site","species")]
adonis(data_dist ~ Isolation.Site, data = sample_meta)
adonis(data_dist ~ species, data = sample_meta)
adonis(data_dist ~ Isolation.Site+species, data = sample_meta)

# spec_names = c("acidophilus","brevis","casei","crispatus","delbrueckii","fermentum","gasseri","helveticus",
#                "iners","jensenii","johnsonii","paracasei","plantarum","reuteri","rhamnosus","salivarius")

palette <- brewer.pal(8,"Dark2")
palette2 <- rep(palette,each=2)

data_nmds_meta_means2 <- merge(data_nmds_meta,aggregate(cbind(mean.x=MDS1,mean.y=MDS2)~species,data_nmds_meta,mean),by="species")

tiff("pFBA_figures/human_ord_species_norm.tiff", width = 6.4, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
p2 = ggplot(data_nmds_meta_means2, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=species,shape=species), alpha = 0.6, size = 2) +
  labs(color = "Species", 
       x=paste(sprintf('MDS axis 1 (%i',contri_1),'%)'), 
       y=paste(sprintf('MDS axis 2 (%i',contri_2),'%)')) +
  scale_colour_manual(name = "Species",
                      values = palette2) +
  #geom_point(aes(x=mean.x,y=mean.y,color=factor(species), shape=species), size=4) +
  scale_shape_manual(name = "Species",
                     values = c(19,17,19,17,19,17,19,17,19,17,19,17,19,17,19,17)) +
  geom_segment(aes(x=mean.x, y=mean.y, xend=MDS1, yend=MDS2,color=factor(species)), size = 1, alpha = 0.2) +
  theme(text = element_text(size=16))
dev.off()

# library(gridExtra)
# tiff("pFBA_figures/Figure_4.tiff", width = 13, height = 4.9, units = 'in', res = 600)
pdf("pFBA_figures/Figure_4.pdf", width = 13, height = 4.9)
grid.newpage()
print(p2, vp=viewport(x = .255, y = 0.5, width = .475, height = 0.95))
print(p1, vp = viewport(x = .75, y = 0.5, width = .475, height = 0.95))
print(grid.text("A", vp = viewport(x = 0.03, y = .95), gp=gpar(fontsize=18)))
print(grid.text("B", vp = viewport(x = 0.52, y = .95), gp=gpar(fontsize=18)))
dev.off()

#### PCA with Biplot ####
library(ggfortify)
data <- data_global
data <- data %>% select(-contains("1"))
data <- data %>% select(-contains("2"))
data_meta <- merge(data, meta_human, by = "row.names")
data_meta_slim <- data_meta
data_meta_slim$Row.names <- NULL
data_meta_slim$Fermentation <- NULL
data_meta_slim$Isolation.Site <- NULL
data_meta_slim$species <- NULL
# View(iris)

# data_meta_alt <- data_meta_slim
data_meta_slim_norm <- as.data.frame(t(apply(data_meta_slim, 1, function(x) ((x - mean(x))/sd(x)) )))
# data_meta <- cbind(data_meta_alt_norm, data_meta$Isolation.Site)
# colnames(data_meta)[length(colnames(data_meta))] <- 'Isolation.Site'

df_pca <- prcomp(data_meta_slim_norm)

df_out <- as.data.frame(df_pca$x)
df_out$group <- data_meta$Isolation.Site
# df_out$group <- sapply( strsplit(as.character(row.names(df)), "_"), "[[", 1 )
# head(df_out)

# percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
# percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

choices = 1:2
u.axis.labs <- paste("PC", choices, sep = "")
u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%%)", 100 * df_pca$sdev[choices]^2/sum(df_pca$sdev^2)))

# df_out$group <- factor(df_out$group, levels = c("virginica", "setosa", "versicolor"))

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group ))
p<-p+geom_point()+theme() + xlab(u.axis.labs[1]) + ylab(u.axis.labs[2])
p

df_out$group <- data_meta$species
p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group ))
p<-p+geom_point()+theme() + xlab(u.axis.labs[1]) + ylab(u.axis.labs[2])
p

library(devtools)
install_github("vqv/ggbiplot")
install.packages('scales')
library(ggbiplot)

ggbiplot(df_pca)

#### Boxplot for all products ####
meta <- meta_human_global
rownames(meta) <- meta$Row.names
meta$Row.names <- NULL

data <- data_global
data <- data %>% select(-contains("1"))
data <- data %>% select(-contains("2"))

data_meta <- merge(data, meta, by = "row.names")
data_meta$Fermentation <- NULL

data_meta$species <- NULL

rownames(data_meta) <- data_meta$Row.names 
data_meta$Row.names <- NULL

data_meta_long <- gather(data_meta, product, likelihood, 1:50, factor_key = TRUE)

tiff("pFBA_figures/products_boxplot.tiff", width = 15, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
ggplot(data = data_meta_long, aes(x = product, y = likelihood, color = Isolation.Site)) + 
  geom_boxplot(aes(x = product, y = likelihood, color = Isolation.Site), outlier.size = .25, position=position_dodge(width = 1)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  labs(x = "Products", y = "Average Pathway Likelihood", color = "Isolation Site") +
  scale_color_discrete(labels=c("gut","oral","vaginal"))
dev.off()

#### Boxplot for all products Normalized ####
meta <- meta_human_global
rownames(meta) <- meta$Row.names
meta$Row.names <- NULL

data <- data_global
data <- data %>% select(-contains("1"))
data <- data %>% select(-contains("2"))

data_meta <- merge(data, meta, by = "row.names")
data_meta$Fermentation <- NULL

data_meta$species <- NULL

rownames(data_meta) <- data_meta$Row.names 
data_meta$Row.names <- NULL

data_meta_alt <- data_meta
data_meta_alt$Isolation.Site <- NULL
data_meta_alt_norm <- as.data.frame(t(apply(data_meta_alt, 1, function(x) ((x - mean(x))/sd(x)) )))
data_meta <- cbind(data_meta_alt_norm, data_meta$Isolation.Site)
colnames(data_meta)[length(colnames(data_meta))] <- 'Isolation.Site'

data_meta_long <- gather(data_meta, product, likelihood, 1:50, factor_key = TRUE)

tiff("pFBA_figures/products_boxplot_norm.tiff", width = 15, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
ggplot(data = data_meta_long, aes(x = product, y = likelihood, color = Isolation.Site)) + 
  geom_boxplot(aes(x = product, y = likelihood, color = Isolation.Site), outlier.size = .25, position=position_dodge(width = 1)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  labs(x = "Products", y = "Average Pathway Likelihood", color = "Isolation Site") +
  scale_color_discrete(labels=c("gut","oral","vaginal"))
dev.off()

#### Lactate Boxplot with normalized data ####
meta <- meta_human_global
rownames(meta) <- meta$Row.names
meta$Row.names <- NULL

data <- data_global
data <- data %>% select(-contains("1"))
data <- data %>% select(-contains("2"))

data_meta <- merge(data, meta, by = "row.names")
data_meta_orig <- data_meta
data_meta$Fermentation <- NULL
data_meta$species <- NULL
data_meta$Isolation.Site <- NULL

rownames(data_meta) <- data_meta$Row.names 
data_meta$Row.names <- NULL

rownames(data_meta) <- data_meta$Row.names 
data_meta$Row.names <- NULL

data_meta_alt <- data_meta
data_meta_alt_norm <- as.data.frame(t(apply(data_meta_alt, 1, function(x) ((x - mean(x))/sd(x)) )))
data_meta_alt_norm <- data_meta_alt_norm %>% select(contains("Lactate"))
data_meta <- cbind(data_meta_alt_norm, data_meta_orig$species)
colnames(data_meta)[length(colnames(data_meta))] <- 'species'
keep <- c('crispatus','fermentum','gasseri','helveticus','iners','jensenii','johnsonii','plantarum','rhamnosus')
data_meta <- data_meta[data_meta$species %in% keep,]

data_meta_long <- gather(data_meta, product, likelihood, 1:2, factor_key = TRUE)

# tiff("pFBA_figures/Figure_6.tiff", width = 6, height = 5, units = 'in', res = 600)
pdf("pFBA_figures/Figure_6.pdf", width = 5, height = 4)
theme_set(theme_bw())
p1 = ggplot(data = data_meta_long, aes(x = species, y = likelihood, color = product)) + 
  geom_boxplot(aes(x = species, y = likelihood, color = product),
               outlier.size = .4, position=position_dodge(width = 1)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1), legend.position = c(.15,.85)) +
  labs(x = "Species", y = "Scaled Production Likelihood", color = "Product") +
  scale_color_discrete(labels=c("D-lactate","L-lactate"))
p1
dev.off()

### Process data to look into Jensenii Production of D-lactate ###
meta <- meta_human_global
rownames(meta) <- meta$Row.names
meta$Row.names <- NULL

data <- data_global

data_meta <- merge(data, meta, by = "row.names")
data_meta_orig <- data_meta
data_meta$Fermentation <- NULL
data_meta$species <- NULL
data_meta$Isolation.Site <- NULL

rownames(data_meta) <- data_meta$Row.names 
data_meta$Row.names <- NULL

rownames(data_meta) <- data_meta$Row.names 
data_meta$Row.names <- NULL

data_meta_alt <- data_meta
data_meta_alt_0 <- data_meta_alt %>% select(-contains("1"))
data_meta_alt_0 <- data_meta_alt_0 %>% select(-contains("2"))
data_meta_alt_norm_0 <- as.data.frame(t(apply(data_meta_alt_0, 1, function(x) ((x - mean(x))/sd(x)) )))
data_meta_alt_1 <- data_meta_alt %>% select(-contains("0"))
data_meta_alt_1 <- data_meta_alt_1 %>% select(-contains("2"))
data_meta_alt_norm_1 <- as.data.frame(t(apply(data_meta_alt_1, 1, function(x) ((x - mean(x))/sd(x)) )))
data_meta_alt_2 <- data_meta_alt %>% select(-contains("0"))
data_meta_alt_2 <- data_meta_alt_2 %>% select(-contains("1"))
data_meta_alt_norm_2 <- as.data.frame(t(apply(data_meta_alt_2, 1, function(x) ((x - mean(x))/sd(x)) )))
data_meta_alt_norm <- cbind(data_meta_alt_norm_0,data_meta_alt_norm_1,data_meta_alt_norm_2)
data_meta_alt_norm <- data_meta_alt_norm %>% select(contains("D_Lactate"))
data_meta <- cbind(data_meta_alt_norm, data_meta_orig$species)
colnames(data_meta)[length(colnames(data_meta))] <- 'species'
data_meta$Row.names <- data_meta_orig$Row.names
data_meta$Row.names <- as.factor(data_meta$Row.names)

data_meta_long <- gather(data_meta, product, likelihood, 1:3, factor_key = TRUE)

data_meta_long <- data_meta_long[data_meta_long$species == "jensenii", ] 

data_meta_long$Row.names <- as.factor(data_meta_long$Row.names)

palette <- brewer.pal(8,"Dark2")
palette2 <- c(palette,"Black","midnightblue")

theme_set(theme_bw())
p2 = ggplot(data = data_meta_long, aes(x = product, y = likelihood, group = Row.names, color = Row.names)) + 
  geom_point(aes(color = Row.names), 
    position = position_jitterdodge(jitter.width = NULL, jitter.height = 0,
                                    dodge.width = 0.7), 
             size = 2, alpha = 0.8) + 
  scale_color_manual(name = "L. jensenii\nGenome ID", values = palette2) +
  scale_x_discrete(labels = c("Glucose","Galactose","Mannose")) +
  labs(x = "Condition", y = "Scaled Pathway Likelihood", title = "D-lactate")
p2

tiff("pFBA_figures/Figure_6.tiff", width = 8, height = 4, units = 'in', res = 600)
grid.arrange(p1, p2, nrow = 1, widths=c(1,1))
dev.off()


#### AUCRF By Isolation Site Normalized ####

## Vaginal Products ##

meta <- read.csv('../Lactobacillus/Data/All_meta_data.csv', sep=",", header=TRUE, colClasses = "character", row.names=1)
meta <- meta[,c("Fermentation","Isolation.Site","species")]

data <- data_global
data <- data %>% select(-contains("1"))
data <- data %>% select(-contains("2"))

data_meta <- merge(data, meta, by = "row.names")
data_meta$Fermentation <- NULL
data_meta$species <- NULL

rownames(data_meta) <- data_meta$Row.names 
data_meta$Row.names <- NULL

data_meta_alt <- data_meta
data_meta_alt$Isolation.Site <- NULL
data_meta_alt_norm <- as.data.frame(t(apply(data_meta_alt, 1, function(x) ((x - mean(x))/sd(x)) )))
data_meta <- cbind(data_meta_alt_norm, data_meta$Isolation.Site)
colnames(data_meta)[length(colnames(data_meta))] <- 'Isolation.Site'

data_meta <- mutate(data_meta, Isolation.Site = ifelse(Isolation.Site == "vaginal_human",1,0))
data_meta$Isolation.Site <- as.factor(data_meta$Isolation.Site)

# Process data to aquire only testable mets #
data_meta_1 <- filter(data_meta, Isolation.Site == 1)
data_meta_1$Isolation.Site <- NULL
data_meta_0 <- filter(data_meta, Isolation.Site == 0)
data_meta_0$Isolation.Site <- NULL
medians_1 <- apply(data_meta_1, 2, function(x) median(x, na.rm=TRUE))
medians_0 <- apply(data_meta_0, 2, function(x) median(x, na.rm=TRUE))
met_keepers <- names(medians_1)[c(medians_1>medians_0 & medians_1>0)]
# met_keepers
data_meta_slim = data_meta[,met_keepers]
data_meta_slim$Isolation.Site = data_meta$Isolation.Site

set.seed(42)
fit <- AUCRF(Isolation.Site~., data = data_meta_slim, ranking=c("MDG","MDA"), ntree = 1000)
fitcv <- AUCRFcv(fit)
OptimalSet(fitcv)

# tiff("pFBA_figures/vaginal_aucrf_norm.tiff", width = 6, height = 5, units = 'in', res = 600)
plot(fitcv)
# dev.off()

# summary(fitcv)

# n = fitcv[["Kopt"]]
n = 7

product_names <- fit[["ranking"]][1:n] %>% names() %>% as.data.frame()
colnames(product_names) <- c("ids")

info_feats <- data_meta_slim[,c(as.character(product_names$ids),'Isolation.Site')]

info_feats_long <- gather(info_feats, product, likelihood, 1:n, factor_key = TRUE)

info_feats_long$product <- sub("0","", info_feats_long$product)
info_feats_long$product <- sub("_","-", info_feats_long$product)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 2
cols = gg_color_hue(n)

tiff("pFBA_figures/vaginal_products_boxplot_norm.tiff", width = 7, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
p1 = ggplot(data = info_feats_long, aes(x = product, y = likelihood, color = Isolation.Site)) + 
  geom_boxplot(aes(x = product, y = likelihood, color = Isolation.Site), outlier.size = .25, position=position_dodge(width = 1)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  scale_color_manual(labels=c("Intestinal/Oral","Vaginal"), values=c(cols[1],cols[2])) + 
  ylim(-2.2,3.8)+
  labs(x = "Products", y = "Scaled Production Likelihood", title = "Vaginal Features", color="Isolation Site")
p1
dev.off()

## Gut products ##

data <- data_global
data <- data %>% select(-contains("1"))
data <- data %>% select(-contains("2"))

data_meta <- merge(data, meta, by = "row.names")
data_meta$Fermentation <- NULL
data_meta$species <- NULL

rownames(data_meta) <- data_meta$Row.names 
data_meta$Row.names <- NULL

data_meta_alt <- data_meta
data_meta_alt$Isolation.Site <- NULL
data_meta_alt_norm <- as.data.frame(t(apply(data_meta_alt, 1, function(x) ((x - mean(x))/sd(x)) )))
data_meta <- cbind(data_meta_alt_norm, data_meta$Isolation.Site)
colnames(data_meta)[length(colnames(data_meta))] <- 'Isolation.Site'

data_meta <- mutate(data_meta, Isolation.Site = ifelse(Isolation.Site == "vaginal_human",1,0))
data_meta$Isolation.Site <- as.factor(data_meta$Isolation.Site)

# Process data to aquire only testable mets #
data_meta_1 <- filter(data_meta, Isolation.Site == 1)
data_meta_1$Isolation.Site <- NULL
data_meta_0 <- filter(data_meta, Isolation.Site == 0)
data_meta_0$Isolation.Site <- NULL
medians_1 <- apply(data_meta_1, 2, function(x) median(x, na.rm=TRUE))
medians_0 <- apply(data_meta_0, 2, function(x) median(x, na.rm=TRUE))
met_keepers <- names(medians_0)[c(medians_0>medians_1 & medians_0>0)]
# met_keepers
data_meta_slim = data_meta[,met_keepers]
data_meta_slim$Isolation.Site = data_meta$Isolation.Site

set.seed(20)
fit <- AUCRF(Isolation.Site~., data = data_meta_slim, ranking=c("MDG","MDA"), ntree = 1000)
fitcv <- AUCRFcv(fit)
# OptimalSet(fitcv)

# tiff("pFBA_figures/gut_aucrf_norm.tiff", width = 6, height = 5, units = 'in', res = 600)
plot(fitcv,xlim=c(0,25))
# dev.off()
fitcv[["AUCcurve"]][["AUC"]][6]

# summary(fitcv)

# n = fitcv[["Kopt"]]
n=8

product_names <- fit[["ranking"]][1:n] %>% names() %>% as.data.frame()
colnames(product_names) <- c("ids")

info_feats <- data_meta[,c(as.character(product_names$ids),'Isolation.Site')]

info_feats_long <- gather(info_feats, product, likelihood, 1:n, factor_key = TRUE)

info_feats_long$product <- sub("0","", info_feats_long$product)
info_feats_long$product <- sub("_","-", info_feats_long$product)

tiff("pFBA_figures/gut_products_boxplot_norm.tiff", width = 7, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
p2 = ggplot(data = info_feats_long, aes(x = product, y = likelihood, color = Isolation.Site)) + 
  geom_boxplot(aes(x = product, y = likelihood, color = Isolation.Site), outlier.size = .25, position=position_dodge(width = 1)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 55, hjust = 1)) +
  scale_color_manual(labels=c("Intestinal/Oral","Vaginal"), values=c(cols[1],cols[2])) +
  ylim(-2.2,3.8)+
  labs(x = "Products", y = "Scaled Production Likelihood", title = "Intestinal/Oral Features", color="Isolation Site")
p2
dev.off()

# tiff("pFBA_figures/Figure_5.tiff", width = 7, height = 4, units = 'in', res = 600)
pdf("pFBA_figures/Figure_5.pdf", width = 7.2, height = 4)
grid.newpage()
print(p2, vp = viewport(x = 1.62/7, y = 0.5, width = 3/7, height = 0.9))
print(p1, vp = viewport(x = 5.1/7, y = 0.5, width = 4/7, height = 0.9))
# grid.arrange(p2, p1, nrow = 1, widths=c(1.05,1.4))
print(grid.text("A", vp = viewport(x = 0.03, y = .95), gp=gpar(fontsize=16)))
print(grid.text("B", vp = viewport(x = 3.2/7, y = .95), gp=gpar(fontsize=16)))
dev.off()

# ## Oral products ##
# 
# data <- data_global
# data <- data %>% select(-contains("1"))
# data <- data %>% select(-contains("2"))
# 
# data_meta <- merge(data, meta, by = "row.names")
# data_meta$Fermentation <- NULL
# data_meta$species <- NULL
# 
# rownames(data_meta) <- data_meta$Row.names 
# data_meta$Row.names <- NULL
# 
# data_meta <- mutate(data_meta, Isolation.Site = ifelse(Isolation.Site == "oral_human",1,0))
# data_meta$Isolation.Site <- as.factor(data_meta$Isolation.Site)
# 
# # Process data to aquire only testable mets #
# data_meta_1 <- filter(data_meta, Isolation.Site == 1)
# data_meta_1$Isolation.Site <- NULL
# data_meta_0 <- filter(data_meta, Isolation.Site == 0)
# data_meta_0$Isolation.Site <- NULL
# medians_1 <- apply(data_meta_1, 2, function(x) median(x, na.rm=TRUE))
# medians_0 <- apply(data_meta_0, 2, function(x) median(x, na.rm=TRUE))
# met_keepers <- names(medians_1)[c(medians_1>medians_0 & medians_1>0)]
# # met_keepers
# data_meta_slim = data_meta[,met_keepers]
# data_meta_slim$Isolation.Site = data_meta$Isolation.Site
# 
# set.seed(42)
# fit <- AUCRF(Isolation.Site~., data = data_meta_slim, ranking=c("MDG","MDA"), ntree = 1000)
# fitcv <- AUCRFcv(fit)
# OptimalSet(fitcv)
# 
# tiff("pFBA_figures/oral_aucrf_norm.tiff", width = 6, height = 5, units = 'in', res = 600)
# plot(fitcv)
# dev.off()
# 
# # summary(fitcv)
# 
# # n = fitcv[["Kopt"]]
# n = 5
# 
# product_names <- fit[["ranking"]][1:n] %>% names() %>% as.data.frame()
# colnames(product_names) <- c("ids")
# 
# info_feats <- data_meta[,c(as.character(product_names$ids),'Isolation.Site')]
# 
# info_feats_long <- gather(info_feats, product, likelihood, 1:n, factor_key = TRUE)
# 
# tiff("pFBA_figures/oral_products_boxplot_norm.tiff", width = 7, height = 5, units = 'in', res = 600)
# theme_set(theme_bw())
# ggplot(data = info_feats_long, aes(x = product, y = likelihood, color = Isolation.Site)) + 
#   geom_boxplot(aes(x = product, y = likelihood, color = Isolation.Site), outlier.size = .25, position=position_dodge(width = 1)) +
#   theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
#   labs(x = "Products", y = "Likelihood", title = "Informative Oral Products")
# dev.off()
#### AUCRF By Isolation Site ####

# Vaginal Products

meta <- read.csv('../Lactobacillus/Data/All_meta_data.csv', sep=",", header=TRUE, colClasses = "character", row.names=1)
meta <- meta[,c("Fermentation","Isolation.Site","species")]

data <- data_global
data <- data %>% select(-contains("1"))
data <- data %>% select(-contains("2"))

data_meta <- merge(data, meta, by = "row.names")
data_meta$Fermentation <- NULL
data_meta$species <- NULL

rownames(data_meta) <- data_meta$Row.names 
data_meta$Row.names <- NULL

data_meta <- mutate(data_meta, Isolation.Site = ifelse(Isolation.Site == "vaginal_human",1,0))
data_meta$Isolation.Site <- as.factor(data_meta$Isolation.Site)

fit <- AUCRF(Isolation.Site~., data = data_meta, ranking=c("MDG","MDA"), ntree = 1000)
fitcv <- AUCRFcv(fit)
OptimalSet(fitcv)

tiff("pFBA_figures/vaginal_aucrf.tiff", width = 6, height = 5, units = 'in', res = 600)
plot(fitcv)
dev.off()

# summary(fitcv)

# n = fitcv[["Kopt"]]
n = 10

product_names <- fit[["ranking"]][1:n] %>% names() %>% as.data.frame()
colnames(product_names) <- c("ids")

info_feats <- data_meta[,c(as.character(product_names$ids),'Isolation.Site')]

info_feats_long <- gather(info_feats, product, likelihood, 1:n, factor_key = TRUE)

tiff("pFBA_figures/vaginal_products_boxplot.tiff", width = 7, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
ggplot(data = info_feats_long, aes(x = product, y = likelihood, color = Isolation.Site)) + 
  geom_boxplot(aes(x = product, y = likelihood, color = Isolation.Site), outlier.size = .25, position=position_dodge(width = 1)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  labs(x = "Products", y = "Likelihood", title = "Informative Vaginal Products")
dev.off()

# Gut products

data <- data_global
data <- data %>% select(-contains("1"))
data <- data %>% select(-contains("2"))

data_meta <- merge(data, meta, by = "row.names")
data_meta$Fermentation <- NULL
data_meta$species <- NULL

rownames(data_meta) <- data_meta$Row.names 
data_meta$Row.names <- NULL

data_meta <- mutate(data_meta, Isolation.Site = ifelse(Isolation.Site == "gut_human",1,0))
data_meta$Isolation.Site <- as.factor(data_meta$Isolation.Site)

fit <- AUCRF(Isolation.Site~., data = data_meta, ranking=c("MDG","MDA"), ntree = 1000)
fitcv <- AUCRFcv(fit)
OptimalSet(fitcv)

tiff("pFBA_figures/gut_aucrf.tiff", width = 6, height = 5, units = 'in', res = 600)
plot(fitcv)
dev.off()

# summary(fitcv)

# n = fitcv[["Kopt"]]
n=10

product_names <- fit[["ranking"]][1:n] %>% names() %>% as.data.frame()
colnames(product_names) <- c("ids")

info_feats <- data_meta[,c(as.character(product_names$ids),'Isolation.Site')]

info_feats_long <- gather(info_feats, product, likelihood, 1:n, factor_key = TRUE)

tiff("pFBA_figures/gut_products_boxplot.tiff", width = 7, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
ggplot(data = info_feats_long, aes(x = product, y = likelihood, color = Isolation.Site)) + 
  geom_boxplot(aes(x = product, y = likelihood, color = Isolation.Site), outlier.size = .25, position=position_dodge(width = 1)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  labs(x = "Products", y = "Likelihood", title = "Informative Gut Products")
dev.off()


# Oral products

data <- data_global
data <- data %>% select(-contains("1"))
data <- data %>% select(-contains("2"))

data_meta <- merge(data, meta, by = "row.names")
data_meta$Fermentation <- NULL
data_meta$species <- NULL

rownames(data_meta) <- data_meta$Row.names 
data_meta$Row.names <- NULL

data_meta <- mutate(data_meta, Isolation.Site = ifelse(Isolation.Site == "oral_human",1,0))
data_meta$Isolation.Site <- as.factor(data_meta$Isolation.Site)

fit <- AUCRF(Isolation.Site~., data = data_meta, ranking=c("MDG","MDA"), ntree = 1000)
fitcv <- AUCRFcv(fit)
OptimalSet(fitcv)

tiff("pFBA_figures/oral_aucrf.tiff", width = 6, height = 5, units = 'in', res = 600)
plot(fitcv)
dev.off()

# summary(fitcv)

n = fitcv[["Kopt"]]

product_names <- fit[["ranking"]][1:n] %>% names() %>% as.data.frame()
colnames(product_names) <- c("ids")

info_feats <- data_meta[,c(as.character(product_names$ids),'Isolation.Site')]

info_feats_long <- gather(info_feats, product, likelihood, 1:n, factor_key = TRUE)

tiff("pFBA_figures/oral_products_boxplot.tiff", width = 7, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
ggplot(data = info_feats_long, aes(x = product, y = likelihood, color = Isolation.Site)) + 
  geom_boxplot(aes(x = product, y = likelihood, color = Isolation.Site), outlier.size = .25, position=position_dodge(width = 1)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  labs(x = "Products", y = "Likelihood", title = "Informative Oral Products")
dev.off()

#### AUCRF By Species ####
species_aucrf_plots <- function(species_name){
  meta <- read.csv('PATRIC_genome_1508_with_meta.csv', sep=",", header=TRUE, colClasses = "character", row.names=1)
  meta <- meta[,c("Fermentation","Isolation.Site","species")]
  
  data <- data_global
  data <- data %>% select(-contains("1"))
  data <- data %>% select(-contains("2"))
  
  data_meta <- merge(data, meta, by = "row.names")
  data_meta$Fermentation <- NULL
  # data_meta$species <- NULL
  data_meta$Isolation.Site <- NULL
  
  rownames(data_meta) <- data_meta$Row.names 
  data_meta$Row.names <- NULL
  
  data_meta <- mutate(data_meta, species = ifelse(species == species_name,1,0))
  data_meta$species <- as.factor(data_meta$species)
  
  set.seed(42)
  fit <- AUCRF(species~., data = data_meta, ranking=c("MDG","MDA"), ntree = 1000)
  set.seed(42)
  fitcv <- AUCRFcv(fit)
  OptimalSet(fitcv)
  fitcv[["RFopt"]][["confusion"]]
  
  tiff(paste0("pFBA_figures/",species_name,"_aucrf.tiff"), width = 6, height = 5, units = 'in', res = 600)
  plot(fitcv, color="black" , ylim=c(0.6, 1.2))
  dev.off()
  
  # n = fitcv[["Kopt"]]
  n = 10
  
  product_names <- fit[["ranking"]][1:n] %>% names() %>% as.data.frame()
  colnames(product_names) <- c("ids")
  
  info_feats <- data_meta[,c(as.character(product_names$ids),'species')]
  
  info_feats_long <- gather(info_feats, product, likelihood, 1:n, factor_key = TRUE)
  
  tiff(paste0("pFBA_figures/",species_name,"_products_boxplot.tiff"), width = 7, height = 5, units = 'in', res = 600)
  theme_set(theme_bw())
  p = ggplot(data = info_feats_long, aes(x = product, y = likelihood, color = species)) + 
    geom_boxplot(aes(x = product, y = likelihood, color = species), outlier.size = .25, position=position_dodge(width = 1)) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
    # coord_cartesian(ylim = c(0,1)) +
    labs(x = "Products", y = "Likelihood", title = "")
  print(p)
  dev.off()
}

# Function Calls #

species_aucrf_plots("acidophilus")
species_aucrf_plots("agilis")
species_aucrf_plots("brevis")
species_aucrf_plots("casei")
species_aucrf_plots("crispatus")
species_aucrf_plots("delbrueckii")
species_aucrf_plots("fermentum")
species_aucrf_plots("gasseri")
species_aucrf_plots("helveticus")
species_aucrf_plots("iners")
species_aucrf_plots("jensenii")
species_aucrf_plots("johnsonii")
species_aucrf_plots("parabuchneri")
species_aucrf_plots("paracasei")
species_aucrf_plots("plantarum")
species_aucrf_plots("reuteri")
species_aucrf_plots("rhamnosus")
species_aucrf_plots("ruminis")
species_aucrf_plots("salivarius")

#### Species by metabolite heatmap ####
# meta <- read.csv('../Lactobacillus/Data/All_meta_data.csv', sep=",", header=TRUE, colClasses = "character", row.names=1)
# meta <- meta[,c("Fermentation","Isolation.Site","species")]

# meta <- meta_human_global
# rownames(meta) <- meta$Row.names
# meta$Row.names <- NULL
# 
# data <- data_global
# data <- data %>% select(-contains("1"))
# data <- data %>% select(-contains("2"))
# 
# data_meta <- merge(data, meta, by = "row.names")
# data_meta$Fermentation <- NULL
# # data_meta$species <- NULL
# data_meta$Isolation.Site <- NULL
# 
# rownames(data_meta) <- data_meta$Row.names 
# data_meta$Row.names <- NULL
# 
# data_meta_long <- gather(data_meta, product, likelihood, 1:50, factor_key = TRUE)
# 
# data_meta_long <- aggregate(likelihood ~ species + product, data_meta_long, median)
# 
# tiff("pFBA_figures/products_heatmap.tiff", width = 5, height = 7, units = 'in', res = 600)
# theme_set(theme_bw())
# ggplot(data = data_meta_long, aes(x = species, y = product)) + 
#   geom_tile(aes(fill = likelihood)) +
#   theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
#   scale_fill_gradient(low = "white", high = "black") +
#   labs(x = "Species", y = "Products", fill = "Median APL")
# dev.off()

#### Heatmap with dendrogram ####

# Read in data
meta <- meta_human_global
rownames(meta) <- meta$Row.names
meta$Row.names <- NULL

data <- data_global
data <- data %>% select(-contains("1"))
data <- data %>% select(-contains("2"))

data_meta <- merge(data, meta, by = "row.names")
data_meta_orig <- data_meta
data_meta$Fermentation <- NULL
data_meta$species <- NULL
data_meta$Isolation.Site <- NULL

rownames(data_meta) <- data_meta$Row.names 
data_meta$Row.names <- NULL

rownames(data_meta) <- data_meta$Row.names 
data_meta$Row.names <- NULL

data_meta_alt <- data_meta
data_meta_alt$Isolation.Site <- NULL
data_meta_alt_norm <- as.data.frame(t(apply(data_meta_alt, 1, function(x) ((x - mean(x))/sd(x)) )))
data_meta <- cbind(data_meta_alt_norm, data_meta_orig$species)
colnames(data_meta)[length(colnames(data_meta))] <- 'species'

data_meta_long2 <- tidyr::gather(data_meta, product, likelihood, 1:50, factor_key = TRUE)

# data_meta_long2_test <- subset(data_meta_long2, product == 'Adenine0' & species == 'plantarum',)
# median(data_meta_long2_test$likelihood)
# boxplot(data_meta_long2_test$likelihood)

data_meta_long <- aggregate(likelihood ~ species + product, data_meta_long2, median)

by_list <- list(data_meta$species)
data_meta$species <- NULL
data_meta_matrix <- aggregate(data_meta, by = by_list, FUN=median)
rownames(data_meta_matrix) <- data_meta_matrix$Group.1
data_meta_matrix$Group.1 <- NULL

# Run clustering
data_meta_matrix_2 <- as.matrix(data_meta_matrix)
rownames(data_meta_matrix_2) <- rownames(data_meta_matrix)
dendro <- as.dendrogram(hclust(d = dist(x = data_meta_matrix_2)))

# Create dendrogram plot
theme_set(default_theme)
dendro.plot <- ggdendrogram(data = dendro, rotate = TRUE) + 
                theme(axis.text.y = element_text(size = 9), 
                      axis.text.x = element_text(size = 0))

# Extract the order of the tips in the dendrogram
species.order <- order.dendrogram(dendro)
# Order the levels according to their position in the cluster
data_meta_long$species <- factor(x = data_meta_long$species,
                               levels = rownames(data_meta_matrix)[species.order], 
                               ordered = TRUE)

data_meta_long$product <- sub("0","", data_meta_long$product)
data_meta_long$product <- sub("_","-", data_meta_long$product)

# Create heatmap plot
theme_set(default_theme)
heatmap.plot <- ggplot(data = data_meta_long, aes(x = product, y = species)) +
                      geom_tile(aes(fill = likelihood)) +
                      theme(axis.text.y = element_blank(),
                            axis.title.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            legend.position = "top",
                            axis.text.x = element_text(angle = 55, hjust = 1),
                            plot.margin = ggplot2::margin(0, 4, 0, 20)) +
                      scale_fill_gradient(low = "white", high = "black") +
                      labs(x = "", y = "", fill = "Median\nSPL")

# theme_set(theme_bw())
data_meta_long2_subset <- subset(data_meta_long2, product == 'Adenine0' & species == 'rhamnosus',select=c('species','product','likelihood'))
p1 = ggplot(data = data_meta_long2_subset, aes(x = "reuteri", y = likelihood, color = "adenine0")) + 
  geom_boxplot(aes(x = "reuteri", y = likelihood, color = "adenine0"),
               outlier.size = .4, position=position_dodge(width = 1)) +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  labs(y = "Scaled Production Likelihood",x="Adenine",title = "L. rhamnosus") +
  scale_color_discrete(labels=c("Adenine"))

# All together
# tiff("pFBA_figures/Figure_3.tiff", width = 10, height = 4, units = 'in', res = 600)
pdf("pFBA_figures/Figure_3.pdf", width = 10, height = 4)
grid.newpage()
theme_set(default_theme)
print(heatmap.plot, vp = viewport(x = 0.47, y = 0.48, width = 0.7, height = 1))
print(dendro.plot, vp = viewport(x = 0.9, y = 0.515, width = 0.17, height = 0.60))
theme_set(theme_bw())
print(p1, vp = viewport(x = 0.06, y = 0.55, width = 0.1, height = 0.66))
print(grid.lines(x=unit(c(0.102,0.153),"npc"),y=unit(c(0.295,0.395),"npc"),gp=gpar(col="#F8766D")))
print(grid.lines(x=unit(c(0.102,0.153),"npc"),y=unit(c(0.80,0.425),"npc"),gp=gpar(col="#F8766D")))
dev.off()

#### Example reactions for likelihoods ####
library(tibble)
data_in <- read.table('../Lactobacillus/Data/reaction_likelihoods_all_genomes.tsv', sep='\t', header=TRUE, colClasses = "character", row.names=1)

data <- sapply(data_in, as.numeric)
rownames(data) <- rownames(data_in)
rm(data_in)
data <- as.data.frame(data)

data_medians <- apply(data, 2, median)
data_medians <- as.data.frame(data_medians)
rownames(data_medians) <- colnames(data)

data_medians <- rownames_to_column(data_medians, var = "rxn_id")
data_medians_arr = arrange(data_medians, desc(data_medians))

# Write CSV in R
write.csv(data_medians_arr, file = "reaction_medians.csv")

#### Reactome likelihood ####
data_in <- read.table('../Lactobacillus/Data/genome_likelihoods.tsv', sep='\t', header=TRUE, colClasses = "character", row.names=1)

data <- sapply(data_in, as.numeric)
rownames(data) <- rownames(data_in)
rm(data_in)
data <- as.data.frame(data)

meta <- meta_human_global
rownames(meta) <- meta$Row.names
meta$Row.names <- NULL

data_meta <- merge(data, meta, by = "row.names")
data_meta$Fermentation <- NULL
data_meta$species <- NULL

rownames(data_meta) <- data_meta$Row.names
data_meta$Row.names <- NULL
data_meta$median <- NULL

data_meta_long <- gather(data_meta, metric, likelihood, 1, factor_key = TRUE)

tiff("pFBA_figures/genome_likelihoods_boxplot.tiff", width = 3, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
ggplot(data = data_meta_long, aes(x = metric, y = likelihood, color = Isolation.Site)) +
  geom_boxplot(aes(x = metric, y = likelihood, color = Isolation.Site), outlier.size = .25, position=position_dodge(width = 1)) +
  theme(axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(0.55,0.68)) +
  labs(x = "Mean", y = "Average Reaction Likelihood", color = "Isolation Site") +
  scale_color_discrete(labels=c("gut","oral","vaginal"))
dev.off()

#### Reactome size ####
data_in <- read.table('../Lactobacillus/Data/genome_sizes.tsv', sep='\t', header=TRUE, colClasses = "character", row.names=1)

data <- sapply(data_in, as.numeric)
rownames(data) <- rownames(data_in)
rm(data_in)
data <- as.data.frame(data)

meta <- meta_human_global
rownames(meta) <- meta$Row.names
meta$Row.names <- NULL

data_meta <- merge(data, meta, by = "row.names")
data_meta$Fermentation <- NULL
data_meta$species <- NULL
data_meta$placeholder <- NULL

rownames(data_meta) <- data_meta$Row.names
data_meta$Row.names <- NULL

data_meta_long <- gather(data_meta, metric, likelihood, 1, factor_key = TRUE)

tiff("pFBA_figures/genome_size_boxplot.tiff", width = 3, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
ggplot(data = data_meta_long, aes(x = metric, y = likelihood, color = Isolation.Site)) +
  geom_boxplot(aes(x = metric, y = likelihood, color = Isolation.Site), outlier.size = .25, position=position_dodge(width = 1)) +
  theme(axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(300,600)) +
  labs(x = "Size", y = "Number of Reactions", color = "Isolation Site") +
  scale_color_discrete(labels=c("gut","oral","vaginal"))
dev.off()

#### Reaction Set Human Likelihoods Ordination ####
data_in <- read.table('../Lactobacillus/Data/reaction_probabilities.csv', sep=',', header=TRUE, colClasses = "character", row.names=1)

data <- sapply(data_in, as.numeric)
rownames(data) <- rownames(data_in)
rm(data_in)
data <- as.data.frame(data)

## Human Gut, Vaginal, Oral
meta <- read.csv('../Lactobacillus/Data/All_meta_data.csv', sep=",", header=TRUE, colClasses = "character", row.names=1)
meta <- meta[,c("Fermentation","Isolation.Site","species")]

target <- c("enviro","food_alcohol","other","unknown","waste","gut_other", "gut_swine", "gut_rodent","food_wheat", "food_milk", "food_veg","",
            "food_meat","food_other","gut_avian","gut_mammal", "blood_human", "bv_human", "other_human", "stomach_human", "milk_human", "gut_lab")
meta_human <- meta
meta_human$Row.names <- rownames(meta_human)
meta_human <- dplyr::filter(meta_human, !(Isolation.Site %in% target))

data <- data[c(meta_human$Row.names), ]
data <- data[complete.cases(data), ]

# Ordination analysis
data_dist <- vegdist(data, method='bray')
data_nmds <- as.data.frame(metaMDS(data_dist, k=2, trymax=100)$points)

## Add Meta-data, then plot
rownames(meta_human) <- meta_human$Row.names
meta_human$Row.names <-  NULL

data_meta <- merge(data, meta_human, by = "row.names")

data_nmds_meta <- merge(data_nmds, meta, by = "row.names")
data_nmds_meta$Fermentation <- as.factor(data_nmds_meta$Fermentation)
data_nmds_meta$Isolation.Site <- as.factor(data_nmds_meta$Isolation.Site)

tiff("pFBA_figures/likes_human_ord_iso_site.tiff", width = 5.9, height = 4.5, units = 'in', res = 600)
theme_set(theme_bw())
ggplot(data_nmds_meta, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=Isolation.Site), alpha = 0.5) +
  labs(color = "Isolation Site", x='NMDS axis 1', y='NMDS axis 2') +
  scale_color_discrete(labels=c("gut","oral","vaginal"))
dev.off()

palette <- brewer.pal(8,"Dark2")
palette2 <- rep(palette,each=2)

tiff("pFBA_figures/likes_human_ord_species.tiff", width = 5.9, height = 4.5, units = 'in', res = 600)
theme_set(theme_bw())
ggplot(data_nmds_meta, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=species, shape=species), alpha = 0.8) +
  labs(color = "Species", x='NMDS axis 1', y='NMDS axis 2') +
  scale_colour_manual(name = "Species",
                      values = palette2)+
  scale_shape_manual(name = "Species",
                     values = c(16,15,16,15,16,15,16,15,16,15,16,15,16,15,16,15))
dev.off()

#### Reaction Set All Likelihoods Ordination ####
data_in <- read.table('../Lactobacillus/Data/reaction_probabilities.csv', sep=',', header=TRUE, colClasses = "character", row.names=1)

data <- sapply(data_in, as.numeric)
rownames(data) <- rownames(data_in)
rm(data_in)
data <- as.data.frame(data)

## Human Gut, Vaginal, Oral
meta <- read.csv('../Lactobacillus/Data/All_meta_data.csv', sep=",", header=TRUE, colClasses = "character", row.names=1)
meta <- meta[,c("Fermentation","Isolation.Site","species")]

data <- data[complete.cases(data), ]

# Ordination analysis
data_dist <- vegdist(data, method='bray')
data_nmds <- as.data.frame(metaMDS(data_dist, k=2, trymax=1000)$points)

data_meta <- merge(data, meta, by = "row.names")

data_nmds_meta <- merge(data_nmds, meta, by = "row.names")
data_nmds_meta$Fermentation <- as.factor(data_nmds_meta$Fermentation)
data_nmds_meta$Isolation.Site <- as.factor(data_nmds_meta$Isolation.Site)

host_associated = c("blood_human","gut_human","gut_avian","gut_rodent","gut_swine","vaginal_human","milk_human","oral_human","gut_other","gut_mammal","other_human","bv_human")
fermented = c("food_veg","food_alcohol","food_milk","gut_lab","food_wheat","food_meat","food_other")
enviromental = c("enviro","waste")

data_nmds_meta = mutate(data_nmds_meta, 
      Isolation.Site = ifelse(Isolation.Site %in% host_associated,"host associated",
                              ifelse(Isolation.Site %in% fermented,"fermented food",
                                     ifelse(Isolation.Site %in% enviromental,"environmental","other/unknown"))))

tiff("pFBA_figures/likes_ord_iso_site.tiff", width = 6.1, height = 4.5, units = 'in', res = 600)
theme_set(theme_bw())
ggplot(data_nmds_meta, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=Isolation.Site), alpha = 0.5) +
  labs(color = "Isolation Site", x='NMDS axis 1', y='NMDS axis 2') 
dev.off()

palette <- brewer.pal(8,"Dark2")
palette2 <- rep(palette,each=2)

tiff("pFBA_figures/likes_ord_species.tiff", width = 5.9, height = 4.5, units = 'in', res = 600)
theme_set(theme_bw())
ggplot(data_nmds_meta, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=species, shape=species), alpha = 0.5) +
  labs(color = "Species", x='NMDS axis 1', y='NMDS axis 2') +
  scale_colour_manual(name = "Species",
                      values = palette2)+
  scale_shape_manual(name = "Species",
                     values = c(16,15,16,15,16,15,16,15,16,15,16,15,16,15,16,15))
dev.off()


#### PATRIC Gene Protein Families Ordinations and Core/Pan plot####
data_in <- read.table('../Lactobacillus/Data/GenomicFeatureFamilyMatrix.csv', sep=',', header=TRUE, colClasses = "character", row.names=1)

data <- as.data.frame(data_in)
data <- sapply(data, as.numeric)
data <- as.data.frame(data)
rownames(data) <- rownames(data_in)
data$Row.names <- as.factor(rownames(data))

meta_human <- meta_human_global
meta_human$Fermentation <- NULL
meta_human$Isolation.Site <- NULL
meta_human$species <- as.factor(meta_human$species)

data_meta <- merge(data, meta_human, by = "Row.names")
rownames(data_meta) <- data_meta$Row.names
data_meta$Row.names <- NULL
data_meta$species <- as.factor(data_meta$species)

data_meta_means <- aggregate(.~species,data_meta,mean)
rownames(data_meta_means) <- data_meta_means$species
data_meta_means$species <- NULL
data_meta_matrix <- as.matrix(data_meta_means)

# Dedrogram with table
# data_meta_long2 <- tidyr::gather(data_meta, product, likelihood, 1:50, factor_key = TRUE)
# 
# data_meta_long <- aggregate(likelihood ~ species + product, data_meta_long2, median)
# 
# by_list <- list(data_meta$species)
# data_meta$species <- NULL
# data_meta_matrix <- aggregate(data_meta, by = by_list, FUN=median)
# rownames(data_meta_matrix) <- data_meta_matrix$Group.1
# data_meta_matrix$Group.1 <- NULL

# Run clustering
# data_meta_matrix_2 <- as.matrix(data_meta_matrix)
# rownames(data_meta_matrix_2) <- rownames(data_meta_matrix)
dendro <- as.dendrogram(hclust(d = dist(x = data_meta_matrix)))

# Create dendrogram plot
theme_set(default_theme)
dendro.plot <- ggdendrogram(data = dendro, rotate = TRUE) + 
  theme(axis.text.y = element_text(size = 9), 
        axis.text.x = element_text(size = 0))

# Table with Species Information
meta <- meta_human_global

meta$Fermentation <- NULL
meta$Isolation.Site <- as.factor(meta$Isolation.Site)
freq_table <-  plyr::count(meta, vars = c('Isolation.Site','species'))
freq_table <- spread(freq_table, Isolation.Site, freq)
freq_table[is.na(freq_table)] <- 0
colnames(freq_table) <- c('Species','Intestinal','Oral','Vaginal')
freq_table <- transform(freq_table,Sum=(Intestinal+Oral+Vaginal)) %>% 
              transform(.,Intestinal=round((Intestinal/Sum)*100,0)) %>% 
              transform(.,Oral=round((Oral/Sum)*100,0)) %>% 
              transform(.,Vaginal=round((Vaginal/Sum)*100,0))
freq_table$Sum <- NULL
freq_table <- as.data.frame(sapply(freq_table,as.character))
freq_table[-1] <- lapply(freq_table[-1], function(x) paste(x,"%", sep=""))


data_in <- read.table('../Lactobacillus/Data/core_pan_data.csv', sep=',', header=TRUE, colClasses = "numeric", row.names=1)
data <- as.data.frame(data_in)
data$genomes <- rownames(data)
data_tall <- gather(data, "type","values", -genomes)
data_tall$type <- factor(data_tall$type, levels = c("pan","core"))
data_tall$genomes <- as.numeric(data_tall$genomes)
class(data_tall$genomes)

# tiff("pFBA_figures/Pan_core_plot.tiff", width = 4, height = 5, units = 'in', res = 600)
theme_set(theme_bw(base_size = 16))
p2 = ggplot(data_tall, aes(x=genomes, y=values, group=type, color=type)) +
  geom_line(size = 1) + 
  labs(y = 'Cross-Genus Protein Families', x='Number of Genomes') +
  scale_colour_manual(labels=c("Pan","Core"), values = c(palette[1],palette[2])) +
  coord_cartesian(ylim = c(50, 1550)) +
  scale_y_continuous(breaks = c(0,400,800,1200,1600)) +
  theme(legend.position = c(0.7, 0.5), plot.margin = ggplot2::margin(5, 10, 5, 5),
        legend.title=element_blank())
p2
# dev.off()

data_in <- read.table('../Lactobacillus/Data/filteredHumanAddedGenomeSizes.csv', sep=',', header=FALSE, colClasses = "numeric")
data <- as.data.frame(t(data_in))
theme_set(theme_bw(base_size = 16))
p1 = ggplot(data, aes(x='',y=V1))+
  geom_boxplot()+
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  coord_cartesian(ylim = c(50, 1550)) +
  scale_y_continuous(breaks = c(0,400,800,1200,1600)) +
  labs(y = "Cross-Genus Protein Families",x="",title = "")
p1

# p1 = ggplot(data = data_meta_long2_subset, aes(x = "reuteri", y = likelihood, color = "adenine0")) + 
#   geom_boxplot(aes(x = "reuteri", y = likelihood, color = "adenine0"),
#                outlier.size = .4, position=position_dodge(width = 1)) +
#   theme(axis.text.x = element_blank(),
#         legend.position = "none",
#         plot.title = element_text(hjust = .5)) +
#   labs(y = "Scaled Production Likelihood",x="Adenine",title = "L. rhamnosus") +
#   scale_color_discrete(labels=c("Adenine"))

colors_set <- c(c("black"),hue_pal()(3))

# tiff("pFBA_figures/Figure_1.tiff", width = 10, height = 5.5, units="in", res=600)
pdf("pFBA_figures/Figure_1.pdf", width = 11, height = 5.3)
grid.newpage()
print(p1, vp=viewport(x = 1.2/10, y = 0.55, width = 2/10, height = 0.89))
print(p2, vp = viewport(x = 4.7/10, y = 0.5, width = 3.6/8, height = 0.85))
print(grid.table(freq_table,vp = viewport(x = 8.4/10, y = .51), 
                 theme = ttheme_minimal(core=list(fg_params=list(hjust=0, x=0.1)),
                          colhead=list(fg_params=list(col=colors_set,hjust=0, x=0.1)),
                          rowhead=list(fg_params=list(col="white")))))
print(grid.text("A", vp = viewport(x = 0.02, y = .95), gp=gpar(fontsize=16)))
print(grid.text("B", vp = viewport(x = 2.4/10, y = .95), gp=gpar(fontsize=16)))
print(grid.text("C", vp = viewport(x = 7.1/10, y = .95), gp=gpar(fontsize=16)))
dev.off()

#### Core Genome ####
# data_in <- read.table('../Lactobacillus/Data/GenomicFeatureFamilyMatrix.csv', sep=',', header=TRUE, colClasses = "character", row.names=1)
# 
# data <- as.data.frame(data_in)
# data <- sapply(data, as.numeric)
# data <- as.data.frame(data)
# rownames(data) <- rownames(data_in)
# data$Row.names <- as.factor(rownames(data))
# 
# 0 %in% data[,1]

#### Old Ordinations and Pan/Core Line plot ####
data <- as.data.frame(data_in)
data <- sapply(data, as.numeric)
rownames(data) <- rownames(data_in)

data_dist <- vegdist(data, method='bray')
data_nmds_df <- as.data.frame(cmdscale(data_dist, k=2, eig=TRUE)$points)
data_nmds <- cmdscale(data_dist, k=2, eig=TRUE)
colnames(data_nmds_df) <- c('MDS1','MDS2')
data_nmds_df$Row.names <- rownames(data)

meta_human <- meta_human_global

data_nmds_meta <- merge(data_nmds_df, meta_human, by = "Row.names")
data_nmds_meta$Fermentation <- as.factor(data_nmds_meta$Fermentation)
data_nmds_meta$Isolation.Site <- as.factor(data_nmds_meta$Isolation.Site)

data_nmds_meta_means <- merge(data_nmds_meta,aggregate(cbind(mean.x=MDS1,mean.y=MDS2)~Isolation.Site,data_nmds_meta,mean),by="Isolation.Site")

contri_1 = as.numeric(round(100*data_nmds$eig[1]/sum(data_nmds$eig[1:20])))
contri_2 = as.numeric(round(100*data_nmds$eig[2]/sum(data_nmds$eig[1:20])))

tiff("pFBA_figures/PAT_GPF_iso_site.tiff", width = 6.4, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
p1 = ggplot(data_nmds_meta_means, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=factor(Isolation.Site)), alpha = 0.6, size = 2) +
  labs(color = "Isolation Site", 
       x=paste(sprintf('MDS axis 1 (%i',contri_1),'%)'), 
       y=paste(sprintf('MDS axis 2 (%i',contri_2),'%)')) +
  scale_color_discrete(labels=c("Gut","Oral","Vaginal")) + 
  geom_point(aes(x=mean.x,y=mean.y,color=factor(Isolation.Site)),size=4) +
  geom_segment(aes(x=mean.x, y=mean.y, xend=MDS1, yend=MDS2,color=factor(Isolation.Site)), size = 1, alpha = 0.2)
p1
dev.off()

palette <- brewer.pal(8,"Dark2")
palette2 <- rep(palette,each=2)

data_nmds_meta_means2 <- merge(data_nmds_meta,aggregate(cbind(mean.x=MDS1,mean.y=MDS2)~species,data_nmds_meta,mean),by="species")

tiff("pFBA_figures/PAT_GPF_species.tiff", width = 6.4, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
p2 = ggplot(data_nmds_meta_means2, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=species,shape=species), alpha = 0.6, size = 2) +
  labs(color = "Species", 
       x=paste(sprintf('MDS axis 1 (%i',contri_1),'%)'), 
       y=paste(sprintf('MDS axis 2 (%i',contri_2),'%)')) +
  scale_colour_manual(name = "Species",
                      values = palette2) +
  geom_point(aes(x=mean.x,y=mean.y,color=factor(species), shape=species), size=4) +
  scale_shape_manual(name = "Species",
                     values = c(19,17,19,17,19,17,19,17,19,17,19,17,19,17,19,17)) +
  geom_segment(aes(x=mean.x, y=mean.y, xend=MDS1, yend=MDS2,color=factor(species)), size = 1, alpha = 0.2)
p2
dev.off()

data_in <- read.table('../Lactobacillus/Data/core_pan_data.csv', sep=',', header=TRUE, colClasses = "numeric", row.names=1)
data <- as.data.frame(data_in)
data$genomes <- rownames(data)
data_tall <- gather(data, "type","values", -genomes)
data_tall$type <- factor(data_tall$type, levels = c("pan","core"))
data_tall$genomes <- as.numeric(data_tall$genomes)
class(data_tall$genomes)

tiff("pFBA_figures/Pan_core_plot.tiff", width = 2.5, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
p3 = ggplot(data_tall, aes(x=genomes, y=values, group=type, color=type)) +
  geom_line(size = 1) + 
  labs(y = 'Number of Cross-Genera Protein Families', x='Number of Genomes') +
  scale_colour_manual(labels=c("Pan","Core"), values = c(palette[1],palette[2])) +
  theme(legend.position = c(0.7, 0.5), plot.margin = ggplot2::margin(5, 10, 5, 5),legend.title=element_blank())
p3
dev.off()

# library(gridExtra)
tiff("pFBA_figures/Figure_1.tiff", width = 15.3, height = 5, units = 'in', res = 600)
grid.arrange(p3, p2, p1, nrow = 1, widths=c(1,2.56,2.56))
dev.off()

