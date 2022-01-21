###############################################################################
#
# One-Way-ANOVA + Tuckey Test for genome size, gene count and PGPT count analysis
# Info: (1) genome size and gene count is given in supplemental tables 2-4 (https://doi.org/10.1101/2021.12.13.472471)
#           Input format: size<tab>genes<tab>class (class means class labels for environments, plant sites, ... see below)
#       (2) PGPT count analysis is based on tabular PGPT annotation output achieved with pgpt_comp_fun_ascii_v2.py -if img -of tab
#           Input format: size<tab>class (size means PGPT count; class means class labels for environments, plant sites, ... see below)
###############################################################################
#install.packages("dplyr")
library(dplyr)
library(ggplot2)
#install.packages("ggpubr")
library(ggpubr)

#################
#ENVIRONMENTS
#################
#file input to choose (keep one option commented befor run)
my_data <- read.delim("ENV_GenomeSizeGeneCount_ANOVA.txt") # OPTION (1) mentioned above
my_data <- read.delim("PGPT_ENV_lev1-Direct_count_ANOVA.txt") #OPTION (2) mentioned above
set.seed(1234)
dplyr::sample_n(my_data, 10)
levels(my_data$env)
my_data$env <- ordered(my_data$env,levels = c("PA", "HA", "AA", "SA", "QA", "FA", "IA")) #class labels
#GENOME SIZE
group_by(my_data, env) %>% summarise(count = n(), mean = mean(size, na.rm = TRUE), sd = sd(size, na.rm = TRUE))
ggboxplot(my_data, x = "env", y = "size", color = "env", palette = c("#7EAB55", "#E9B48A", "#F6E4D6", "#79582C", "#436FBE", "#FAE6A3", "#A5A6A6"),order = c("PA", "HA", "AA", "SA", "QA", "FA", "IA"),ylab = "Genome Size", xlab = "Environment")
#ggline(my_data, x = "env", y = "size", color = "env", palette = c("#7EAB55", "#E9B48A", "#F6E4D6", "#79582C", "#436FBE", "#FAE6A3", "#A5A6A6"),order = c("PA", "HA", "AA", "SA", "QA", "FA", "IA"),ylab = "Genome Size", xlab = "Environment")
res.aov <- aov(size ~ env, data = my_data)
summary(res.aov)
tt <- TukeyHSD(res.aov)
str(tt)
tt$env[,"p adj"]
print(tt,digits=22)
#GENES
group_by(my_data, env) %>% summarise(count = n(), mean = mean(genes, na.rm = TRUE), sd = sd(genes, na.rm = TRUE))
ggboxplot(my_data, x = "env", y = "size", color = "env", palette = c("#7EAB55", "#E9B48A", "#F6E4D6", "#79582C", "#436FBE", "#FAE6A3", "#A5A6A6"),order = c("PA", "HA", "AA", "SA", "QA", "FA", "IA"),ylab = "Gene Count", xlab = "Environment")
#ggline(my_data, x = "env", y = "size", color = "env", palette = c("#7EAB55", "#E9B48A", "#F6E4D6", "#79582C", "#436FBE", "#FAE6A3", "#A5A6A6"),order = c("PA", "HA", "AA", "SA", "QA", "FA", "IA"),ylab = "Genome Size", xlab = "Environment")
res.aov <- aov(genes ~ env, data = my_data)
summary(res.aov)
tt <- TukeyHSD(res.aov)
str(tt)
tt$env[,"p adj"]
print(tt,digits=22)

#################
#PSITES/PSPHERES
#################
#file input to choose (keep one option commented befor run)
my_data <- read.delim("PSITE_GenomeSizeGeneCount_ANOVA.txt") # OPTION (1) mentioned above
my_data <- read.delim("PGPT_PSITE_lev1-Direct_count_ANOVA.txt") #OPTION (2) mentioned aboveset.seed(1234)
dplyr::sample_n(my_data, 10)
levels(my_data$psite)
my_data$psite <- ordered(my_data$psite,levels = c("PA_RHIZ", "PA_PHYL", "PA_FLOW", "PA_FRUIT", "PA_SEED"))
#GENOME SIZE
group_by(my_data, psite) %>% summarise(count = n(), mean = mean(size, na.rm = TRUE), sd = sd(size, na.rm = TRUE))
ggboxplot(my_data, x = "psite", y = "size", color = "psite", palette = c("#A3795A", "#768967", "#9471B8", "#C35552", "#A19056"),order = c("PA_RHIZ", "PA_PHYL", "PA_FLOW", "PA_FRUIT", "PA_SEED"),ylab = "Genome Size", xlab = "Plant Sphere")
res.aov <- aov(size ~ psite, data = my_data)
summary(res.aov)
tt <- TukeyHSD(res.aov)
str(tt)
tt$psite[,"p adj"]
print(tt,digits=22)
#GENES
group_by(my_data, psite) %>% summarise(count = n(), mean = mean(genes, na.rm = TRUE), sd = sd(size, na.rm = TRUE))
ggboxplot(my_data, x = "psite", y = "genes", color = "psite", palette = c("#A3795A", "#768967", "#9471B8", "#C35552", "#A19056"),order = c("PA_RHIZ", "PA_PHYL", "PA_FLOW", "PA_FRUIT", "PA_SEED"),ylab = "Gene Count", xlab = "Plant Sphere")
res.aov <- aov(genes ~ psite, data = my_data)
summary(res.aov)
tt <- TukeyHSD(res.aov)
str(tt)
tt$psite[,"p adj"]
print(tt,digits=22)

#################
#PPHEN
#################
#file input to choose (keep one option commented befor run)
my_data <- read.delim("PPHEN_GenomeSizeGeneCount_ANOVA.txt") # OPTION (1) mentioned above
my_data <- read.delim("PGPT_PPHEN_lev1-Direct_count_ANOVA.txt") #OPTION (2) mentioned aboveset.seed(1234)
set.seed(1234)
dplyr::sample_n(my_data, 10)
levels(my_data$pphen)
my_data$pphen <- ordered(my_data$pphen,levels = c("PA_SYMB", "PA_PATH"))
#GENOME SIZE
group_by(my_data, pphen) %>% summarise(count = n(), mean = mean(size, na.rm = TRUE), sd = sd(size, na.rm = TRUE))
ggboxplot(my_data, x = "pphen", y = "size", color = "pphen", palette = c("#00B04F", "#F38D10"),order = c("PA_SYMB", "PA_PATH"),ylab = "Genome Size", xlab = "Plant Phenotype")
res.aov <- aov(size ~ pphen, data = my_data)
summary(res.aov)
tt <- TukeyHSD(res.aov)
str(tt)
tt$pphen[,"p adj"]
print(tt,digits=22)
#GENES
group_by(my_data, pphen) %>% summarise(count = n(), mean = mean(genes, na.rm = TRUE), sd = sd(genes, na.rm = TRUE))
ggboxplot(my_data, x = "pphen", y = "genes", color = "pphen", palette = c("#00B04F", "#F38D10"),order = c("PA_SYMB", "PA_PATH"),ylab = "Gene Count", xlab = "Plant Phenotype")
res.aov <- aov(genes ~ pphen, data = my_data)
summary(res.aov)
tt <- TukeyHSD(res.aov)
str(tt)
tt$phen[,"p adj"]
print(tt,digits=22)


###############################################################################
#
# SCOARY Result-based HEATMAP
# Info: Takes as input the output of Scoary (supplemental tables 5-6 https://doi.org/10.1101/2021.12.13.472471), that has then to be converted into the following format:
#       Input format: CLASS<tab>Class1<tab>Class2 ...ClassN
#                     PGPT1<tab>value<tab>value ...
#                     PGPT2 ...
#      (Class1 to ClassN states for the respective class labels for environments, plant sites, ... see above)
#      (Valuse are log Odds Ratios) 
#     (Scoary input files achieved by pgpt_comp_fun_ascii_v2.py if img -of roary and supplemental tables 2-4 https://doi.org/10.1101/2021.12.13.472471)
#
###############################################################################
library("RColorBrewer")
#install.packages("gplots")
library("gplots")
#install.packages("varhandle")
library("varhandle")
getwd()

col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)

tab = read.csv("PPHEN_LEV3_lorZscale.txt",header = TRUE, sep = "\t", dec = ".") #exemplarily
t <- as.matrix(tab[, -1])
storage.mode(t) <- "numeric"
rownames(t) <- tab$CLASS
#1-minus Pearson correlation distance with average linkage
heatmap.2(t, distfun = function(x) as.dist(1-cor(t(x))), hclustfun = function(x) hclust(x, method="average"), 
          scale = "none", col = bluered(100), trace = "none", density.info = "density",cexRow=0.7,cexCol=1,
          margins = c(7,20),lhei=c(4,15),lwid=c(2,7),keysize=0.75, key.par = list(cex=0.5))

#Euclidean distance with Ward's linkage
heatmap.2(t, distfun = function(x) dist(x, method="euclidean"), hclustfun = function(x) hclust(x, method="ward.D"), 
          scale = "none", col = bluered(100), trace = "none", density.info = "density",cexRow=0.7,cexCol=1,
          margins = c(7,20),lhei=c(4,15),lwid=c(2,7),keysize=0.75, key.par = list(cex=0.5))

###############################################################################
#
# PCOA of PGPT counts of Ontology level 2 for PPHEN (classes PA_SYMB, PA_PATH)
# Info: Takes as input the output of pgpt_comp_fun_ascii_v2.py -if img -of tab and and of supplemental tables 2-4 https://doi.org/10.1101/2021.12.13.472471in converted to format:
#       Input format: IMG_SAMPLE_ID<tab>GENOME_SIZE><tab>Level2_BIO_FERT...<tab>Level2_BIO-REMED...<tab>...<tab>CLASS
#                     XYXYXYXY<tab>6000000<tab>238<tab>167<tab>...<tab>PA_SYMB
#                     XXXXXXXX<tab> ...
#       (Values are given as pure PGPT count per PGPT subclass)
#
###############################################################################
library("RColorBrewer")
#install.packages("gplots")
library("gplots")
#install.packages("varhandle")
library("varhandle")
#install.packages("FactoMineR")
library("FactoMineR")
#install.packages("factoextra")
library("factoextra")
#install.packages("corrplot")
library("corrplot")
#install.packages("pca3d")
library("pca3d")
#install.packages("vegan")
library("vegan")
#install.packages("ecodist")
library("ecodist")
#install.packages("ape")
library("ape")
#install.packages("cluster")
library("cluster")
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

#Data acquisition and formatting
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)

tab = read.csv("PGPT-COUNT_LEV2_PPHEN.txt",header = TRUE, sep = "\t", dec = ".") #exemplarily
t <- as.matrix(tab[,-1])
storage.mode(t) <- "numeric"
rownames(t) <- tab$IMG_SAMPLE_ID

res.pca <- PCA(tab[3:10], graph = FALSE) #tab[3:10] might be extended for other ontology levels, as done above
print(res.pca)
#Eigen values
eig.val <- get_eigenvalue(res.pca)
eig.val
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 80))
var <- get_pca_var(res.pca)
var
corrplot(var$cos2, is.corr=FALSE)
#PCA plotting
pca <- prcomp(tab[3:10], scale.=TRUE, center.=TRUE) #tab[3:10] might be extended for other ontology levels, as done above
gr <- factor(tab[,11]) # should call last column in file, so if changes done above, adapt 
summary(gr)
#3d (uncomment next line)
#pca3d(pca, group=gr, show.plane=FALSE)
#2d
#Option 1
pca2d(pca, group=gr, legend="topleft", biplot=TRUE, biplot.vars=8, 
      palette=c("#F38D10", "#00B04F"),show.centroids = TRUE,
      show.shadows = FALSE, radius = 0.7,
      show.ellipses = TRUE, ellipse.ci = 0.95)
#Option 2
ggbiplot(pca,ellipse=TRUE,  labels=rownames(tab[3:10]), groups=gr, 
         choices=c(1,2),obs.scale = 0.5, var.scale = 0.5,var.axes=TRUE)
