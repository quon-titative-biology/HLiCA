library(Seurat)
setwd("project/Liver_Cell_Atlas/Analysis_02/")


Mesenchyme <- read.csv("02_Lineages_itr01/01_Mesenchyme/Mesenchyme_annotations.csv", row.names = 1)
Endothelia <- read.csv("02_Lineages_itr01/02_Endothelia/Endothelia_annotations.csv", row.names = 1)
Cholangiocyte <- read.csv("02_Lineages_itr01/03_Cholangiocytes/Cholangiocyte_annotations.csv", row.names = 1)
Myeloid_cell <- read.csv("02_Lineages_itr01/04_Myeloid_Cell/Myeloid_Cell_annotations.csv", row.names = 1)
Lymphocyte <- read.csv("02_Lineages_itr01/05_Lymphocyte/Lymphocyte_annotations.csv", row.names = 1)
Hepatocyte <- read.csv("02_Lineages_itr01/06_Hepatocyte/Hepatocyte_annotations.csv", row.names = 1)


beta.annotations <- rbind(Mesenchyme, Endothelia, Cholangiocyte, Myeloid_cell, Lymphocyte, Hepatocyte)
alpha.annotations <- read.csv("01_Alpha_Annotations.csv", row.names = 1)

df <- merge(alpha.annotations, beta.annotations, by = 'row.names', all = TRUE) 
rownames(df) <- df$Row.names

# Test

# Mesenchyme
mes_preclean <- readRDS("02_Lineages_itr01/01_Mesenchyme/01_PreClean/chosen_parameters/Mesenchyme_unclean.rds")

test <- mes_preclean$seu
test@meta.data$Alpha.Annotations <- NULL
test@meta.data$Beta.Annotation.Lineage <- NULL
test@meta.data$Potential.Doublets <- NULL

test <- AddMetaData(test, df)

all.equal(as.factor(test$Alpha.Annotations), mes_preclean$seu$Alpha.Annotations)
all.equal(as.factor(test$Potential.Doublets), mes_preclean$seu$Potential.Doublets)



mes_clean <- readRDS("02_Lineages_itr01/01_Mesenchyme/02_Clean/chosen_parameters/Mesenchyme_clean.rds")

test <- mes_clean$seu
test@meta.data$Alpha.Annotations <- NULL
test@meta.data$Beta.Annotation.Lineage <- NULL
test@meta.data$Potential.Doublets <- NULL

test <- AddMetaData(test, df)

all.equal(as.factor(test$Alpha.Annotations), mes_clean$seu$Alpha.Annotations)
all.equal(as.factor(test$Potential.Doublets), mes_clean$seu$Potential.Doublets)
all.equal(factor(test$Beta.Annotation.Lineage, levels = levels(mes_clean$seu$Beta.Annotation.Lineage)), mes_clean$seu$Beta.Annotation.Lineage)
all.equal(factor(test$Beta.Annotation.SubLineage, levels = levels(mes_clean$seu$Beta.Annotation.SubLineage)), mes_clean$seu$Beta.Annotation.SubLineage)
all.equal(as.factor(test$Potential.Doublets.2), mes_clean$seu$Potential.Doublets.2)



# Endothelia
endo_preclean <- readRDS("02_Lineages_itr01/02_Endothelia/01_PreClean/chosen_parameters/Endothelia_unclean.rds")

test <- endo_preclean$seu
test@meta.data$Alpha.Annotations <- NULL
test@meta.data$Beta.Annotation.Lineage <- NULL
test@meta.data$Potential.Doublets <- NULL

test <- AddMetaData(test, df)

all.equal(as.factor(test$Alpha.Annotations), endo_preclean$seu$Alpha.Annotations)
all.equal(as.factor(test$Potential.Doublets), endo_preclean$seu$Potential.Doublets)



endo_clean <- readRDS("02_Lineages_itr01/02_Endothelia/02_Clean/chosen_parameters/Endothelia_clean.rds")

test <- endo_clean$seu
test@meta.data$Alpha.Annotations <- NULL
test@meta.data$Beta.Annotation.Lineage <- NULL
test@meta.data$Potential.Doublets <- NULL

test <- AddMetaData(test, df)

all.equal(as.factor(test$Alpha.Annotations), endo_clean$seu$Alpha.Annotations)
all.equal(as.factor(test$Potential.Doublets), endo_clean$seu$Potential.Doublets)
all.equal(factor(test$Beta.Annotation.Lineage, levels = levels(endo_clean$seu$Beta.Annotation.Lineage)), endo_clean$seu$Beta.Annotation.Lineage)
all.equal(factor(test$Beta.Annotation.SubLineage, levels = levels(endo_clean$seu$Beta.Annotation.SubLineage)), endo_clean$seu$Beta.Annotation.SubLineage)
all.equal(as.factor(test$Potential.Doublets.2), endo_clean$seu$Potential.Doublets.2)


# Cholangiocytes
cho_preclean <- readRDS("02_Lineages_itr01/03_Cholangiocytes/01_PreClean/chosen_parameters/Cholangiocyte_unclean.rds")

test <- cho_preclean$seu
test@meta.data$Alpha.Annotations <- NULL
test@meta.data$Beta.Annotation.Lineage <- NULL
test@meta.data$Potential.Doublets <- NULL

test <- AddMetaData(test, df)

all.equal(as.factor(test$Alpha.Annotations), cho_preclean$seu$Alpha.Annotations)
all.equal(as.factor(test$Potential.Doublets), cho_preclean$seu$Potential.Doublets)



cho_clean <- readRDS("02_Lineages_itr01/03_Cholangiocytes/02_Clean/chosen_parameters/Cholangiocyte_clean.rds")

test <- cho_clean$seu
test@meta.data$Alpha.Annotations <- NULL
test@meta.data$Beta.Annotation.Lineage <- NULL
test@meta.data$Potential.Doublets <- NULL

test <- AddMetaData(test, df)

all.equal(as.factor(test$Alpha.Annotations), cho_clean$seu$Alpha.Annotations)
all.equal(as.factor(test$Potential.Doublets), cho_clean$seu$Potential.Doublets)
all.equal(factor(test$Beta.Annotation.Lineage, levels = levels(cho_clean$seu$Beta.Annotation.Lineage)), cho_clean$seu$Beta.Annotation.Lineage)
all.equal(factor(test$Beta.Annotation.SubLineage, levels = levels(cho_clean$seu$Beta.Annotation.SubLineage)), cho_clean$seu$Beta.Annotation.SubLineage)
all.equal(factor(test$Potential.Doublets.2, levels = levels(cho_clean$seu$Potential.Doublets.2)), cho_clean$seu$Potential.Doublets.2)


# Myeloid Cells
myeloid_preclean <- readRDS("02_Lineages_itr01/04_Myeloid_Cell/01_PreClean/chosen_parameters/Myeloid_Cells_unclean.rds")

test <- myeloid_preclean$seu
test@meta.data$Alpha.Annotations <- NULL
test@meta.data$Beta.Annotation.Lineage <- NULL
test@meta.data$Potential.Doublets <- NULL

test <- AddMetaData(test, df)

all.equal(as.factor(test$Alpha.Annotations), myeloid_preclean$seu$Alpha.Annotations)
all.equal(as.factor(test$Potential.Doublets), myeloid_preclean$seu$Potential.Doublets)



myeloid_clean <- readRDS("02_Lineages_itr01/04_Myeloid_Cell/02_Clean/chosen_parameters/Myeloid_Cell_clean.rds")

test <- myeloid_clean$seu
test@meta.data$Alpha.Annotations <- NULL
test@meta.data$Beta.Annotation.Lineage <- NULL
test@meta.data$Potential.Doublets <- NULL

test <- AddMetaData(test, df)

all.equal(as.factor(test$Alpha.Annotations), myeloid_clean$seu$Alpha.Annotations)
all.equal(as.factor(test$Potential.Doublets), myeloid_clean$seu$Potential.Doublets)
all.equal(factor(test$Beta.Annotation.Lineage, levels = levels(myeloid_clean$seu$Beta.Annotation.Lineage)), myeloid_clean$seu$Beta.Annotation.Lineage)
all.equal(factor(test$Beta.Annotation.SubLineage, levels = levels(myeloid_clean$seu$Beta.Annotation.SubLineage)), myeloid_clean$seu$Beta.Annotation.SubLineage)
all.equal(factor(test$Potential.Doublets.2, levels = levels(myeloid_clean$seu$Potential.Doublets.2)), myeloid_clean$seu$Potential.Doublets.2)


# Lmyphocyte 
lmp_preclean <- readRDS("02_Lineages_itr01/05_Lymphocyte/01_PreClean/chosen_parameters/Lymphocyte_unclean.rds")

test <- lmp_preclean$seu
test@meta.data$Alpha.Annotations <- NULL
test@meta.data$Beta.Annotation.Lineage <- NULL
test@meta.data$Potential.Doublets <- NULL

test <- AddMetaData(test, df)

all.equal(as.factor(test$Alpha.Annotations), lmp_preclean$seu$Alpha.Annotations)
all.equal(as.factor(test$Potential.Doublets), lmp_preclean$seu$Potential.Doublets)



lmp_clean <- readRDS("02_Lineages_itr01/05_Lymphocyte/02_Clean/chosen_parameters/Lymphocyte_clean.rds")

test <- lmp_clean$seu
test@meta.data$Alpha.Annotations <- NULL
test@meta.data$Beta.Annotation.Lineage <- NULL
test@meta.data$Potential.Doublets <- NULL

test <- AddMetaData(test, df)

all.equal(as.factor(test$Alpha.Annotations), lmp_clean$seu$Alpha.Annotations)
all.equal(as.factor(test$Potential.Doublets), lmp_clean$seu$Potential.Doublets)
all.equal(factor(test$Beta.Annotation.Lineage, levels = levels(lmp_clean$seu$Beta.Annotation.Lineage)), lmp_clean$seu$Beta.Annotation.Lineage)
all.equal(factor(test$Beta.Annotation.SubLineage, levels = levels(lmp_clean$seu$Beta.Annotation.SubLineage)), lmp_clean$seu$Beta.Annotation.SubLineage)
all.equal(factor(test$Potential.Doublets.2, levels = levels(lmp_clean$seu$Potential.Doublets.2)), lmp_clean$seu$Potential.Doublets.2)


# Hpeatocytes 
hep_preclean <- readRDS("02_Lineages_itr01/06_Hepatocyte/01_PreClean/chosen_parameters/Hepatocyte_unclean.rds")

test <- hep_preclean$seu
test@meta.data$Alpha.Annotations <- NULL
test@meta.data$Beta.Annotation.Lineage <- NULL
test@meta.data$Potential.Doublets <- NULL

test <- AddMetaData(test, df)

all.equal(as.factor(test$Alpha.Annotations), hep_preclean$seu$Alpha.Annotations)



hep_preclean2 <- readRDS("02_Lineages_itr01/06_Hepatocyte/02_PreClean_2/chosen_parameters/Hepatocyte_unclean_2.rds")

test <- hep_preclean2$seu
test@meta.data$Alpha.Annotations <- NULL
test@meta.data$Beta.Annotation.Lineage <- NULL
test@meta.data$Potential.Doublets <- NULL

test <- AddMetaData(test, df)

all.equal(as.factor(test$Alpha.Annotations), hep_preclean2$seu$Alpha.Annotations)
all.equal(as.factor(test$Potential.Doublets), factor(hep_preclean2$seu$Potential.Doublets))
all.equal(factor(test$Beta.Annotation.Lineage), factor(hep_preclean2$seu$Beta.Annotation.Lineage))



hep_clean <- readRDS("02_Lineages_itr01/06_Hepatocyte/03_Clean/chosen_parameters/Hepatocyte_clean.rds")

test <- hep_clean$seu
test@meta.data$Alpha.Annotations <- NULL
test@meta.data$Beta.Annotation.Lineage <- NULL
test@meta.data$Potential.Doublets <- NULL

test <- AddMetaData(test, df)

all.equal(as.factor(test$Alpha.Annotations), hep_clean$seu$Alpha.Annotations)
all.equal(factor(test$Potential.Doublets), factor(hep_clean$seu$Potential.Doublets))
all.equal(factor(test$Beta.Annotation.Lineage, levels = levels(hep_clean$seu$Beta.Annotation.Lineage)), hep_clean$seu$Beta.Annotation.Lineage)
all.equal(factor(test$Beta.Annotation.SubLineage, levels = levels(hep_clean$seu$Beta.Annotation.SubLineage)), hep_clean$seu$Beta.Annotation.SubLineage)
all.equal(factor(test$Potential.Doublets.2, levels = levels(hep_clean$seu$Potential.Doublets.2)), hep_clean$seu$Potential.Doublets.2)



table(df$Alpha.Annotations)
table(df$Beta.Annotation.Lineage)
table(df$Potential.Doublets)
table(df$Beta.Annotation.SubLineage)
table(df$Potential.Doublets.2)


sum(duplicated(table(df$Alpha.Annotations)))
sum(duplicated(table(df$Beta.Annotation.Lineage)))
sum(duplicated(table(df$Potential.Doublets)))
sum(duplicated(table(df$Beta.Annotation.SubLineage)))
sum(duplicated(table(df$Potential.Doublets.2)))

seu <- readRDS("01_Broad_Cell_Annotations/seu.rds")
dim(df)[1] == length(colnames(seu))


write.csv(df, "~/project/Liver_Cell_Atlas/Analysis_02/02_Beta_Annotations.csv")
