
## Load Libraries
library(here)

library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

library(sp)
library(scales)
library(viridis)
library(colorspace)

library(igraph)
library(tiff)
library(tidyr)
library(Seurat)

library(scatterpie)




source("scripts/00_pretty_plots.R")


#############
## DEG between cholangiocyte types
#############

annotation_C107_8um<-read.csv(here("data/annotation_C107_8um_additional.csv"))

localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/reseq/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

head(object)

rownames(annotation_C107_8um)<-annotation_C107_8um$bin
annotation_C107_8um<-annotation_C107_8um[match(colnames(object), rownames(annotation_C107_8um)),]
identical(colnames(object),annotation_C107_8um$bin)

object<-AddMetaData(object,annotation_C107_8um)

cholangiocyte_bins<-subset(object, subset = final_anno == "intrahepatic cholangiocyte")


### resize the ducts
bin_count <- table(cholangiocyte_bins$bile_duct_cluster)

cholangiocyte_bins$bile_duct_size <- sapply(cholangiocyte_bins$bile_duct_cluster, function(cluster) {
  bin_num <- bin_count[as.character(cluster)]
  if (bin_num <= 3) {
    return("Bile ductules")
  } else if (bin_num < 30) {
    return("Interlobular duct")
  } else if (bin_num >= 30) {
    return("Septal duct")
  }
})

## how does that translate to mean diameter?
chol_coords <- as.data.frame(cholangiocyte_bins@images$slice1$centroids@coords)
rownames(chol_coords)<-colnames(cholangiocyte_bins)

# plot bile ducts annotated by size
chol_coords$bin<-rownames(chol_coords)
plt_cholangiocyte_spatial<-merge(chol_coords,cholangiocyte_bins@meta.data, by="bin")

## only ducts with more than 3 - for those less assign 16um since it will be two bins max
hulls <- plt_cholangiocyte_spatial %>%
  group_by(bile_duct_cluster, bile_duct_size) %>%
  filter(n() >= 3) %>%  slice(chull(x, y)) %>%  ungroup()

hull_diameters <- hulls %>%  group_by(bile_duct_cluster,bile_duct_size) %>%
  summarise(diameter = if(n() > 1) max(dist(cbind(x, y))) else NA_real_,.groups = "drop")

plt_cholangiocyte_spatial<-merge(plt_cholangiocyte_spatial, hull_diameters[,c("bile_duct_cluster","diameter")], by="bile_duct_cluster", all.x=T)

## duct 86 is a square and the diameter is 19.02984 but should be 16 since it is 8um bins
#So scale all diameters by this factor
16/19.02984
plt_cholangiocyte_spatial$diameter_adjusted<-plt_cholangiocyte_spatial$diameter*(16/19.02984)

## add back in ducts with only 2-3 bins (diameter = 16)
#2 bins
bins2<-names(bin_count[which(bin_count<3 & bin_count>1)])
plt_cholangiocyte_spatial$diameter_adjusted[which(plt_cholangiocyte_spatial$bile_duct_cluster%in%bins2)]<-16

mean_diameter_by_size <- plt_cholangiocyte_spatial %>%
  group_by(bile_duct_size) %>%
  summarise(min_diameter = min(diameter_adjusted, na.rm = TRUE),
            mean_diameter = mean(diameter_adjusted, na.rm = TRUE),
            max_diameter = max(diameter_adjusted, na.rm = TRUE),
            n = n(),.groups = "drop")
mean_diameter_by_size


ggplot(plt_cholangiocyte_spatial, aes(x,y, color=diameter_adjusted))+geom_point()
ggplot(plt_cholangiocyte_spatial, aes(bile_duct_size, diameter_adjusted))+geom_boxplot()


plt_cholangiocyte_spatial$bile_duct_size_by_diameter <- sapply(1:nrow(plt_cholangiocyte_spatial), function(x) {
  bin_dia <- plt_cholangiocyte_spatial$diameter_adjusted[x]
  if(is.na(bin_dia)){NA}else{
    if (bin_dia <= 20) {"Bile ductules"  
    } else if (bin_dia > 20 & bin_dia < 100 ) { "Interlobular duct"
    } else if (bin_dia >= 100) { "Septal duct"  }
  }
  })

save(plt_cholangiocyte_spatial, file="data/HLiCA_duct_size_diameter.RData")

ggplot(plt_cholangiocyte_spatial, aes(bile_duct_size_by_diameter, diameter_adjusted))+geom_boxplot()

plt_cholangiocyte_spatial<-plt_cholangiocyte_spatial[which(!(is.na(plt_cholangiocyte_spatial$bile_duct_size_by_diameter))),]

mean_diameter_by_size <- plt_cholangiocyte_spatial %>%
  group_by(bile_duct_size_by_diameter) %>%
  summarise(min_diameter = min(diameter_adjusted, na.rm = TRUE),
            mean_diameter = mean(diameter_adjusted, na.rm = TRUE),
            max_diameter = max(diameter_adjusted, na.rm = TRUE),
            n = n(),.groups = "drop")
mean_diameter_by_size

# ## adjust si
# bin_size <- 13.4543
# unit_per_um <- bin_size / 8
# scale_length_um <- 100
# scale_length_plot <- scale_length_um * unit_per_um
# 
# 
# y_start<-7000
# x_start<-(-3000)

plot_bins<-ggplot() +
  geom_point(data = plt_cholangiocyte_spatial,aes(x = y, y = -x, color = bile_duct_size_by_diameter),size = 0.25) +
  geom_rect(aes(xmin = min(plt_cholangiocyte_spatial$y), xmax = max(plt_cholangiocyte_spatial$y),
                ymin = (-min(plt_cholangiocyte_spatial$x)), ymax = (-max(plt_cholangiocyte_spatial$x))),
            fill = NA,color = "black",linewidth = 1) +
  coord_fixed() +   theme_void()
plot_bins




## mucus secreting distribution
load(here("data/mucus_secreting.RData"))

plt_cholangiocyte_spatial$type<-"non-mucus secreating"
plt_cholangiocyte_spatial$type[which(plt_cholangiocyte_spatial$bin%in%MUC_cholangiocyte)]<-"mucus secreting"

df_per<-as.data.frame(table(plt_cholangiocyte_spatial$type, plt_cholangiocyte_spatial$bile_duct_size_by_diameter))

df_plot <- df_per %>%
  group_by(Var2) %>%
  mutate(prop = Freq / sum(Freq)) %>%
  filter(Var1 == "mucus secreting")

ggplot(df_plot, aes(x = Var2, y = prop)) +
  geom_col(width = 0.6, fill = "steelblue") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Secreting cholangiocytes (%)") +
  theme_classic()


df_plot<-df_per %>%  group_by(Var2) %>% mutate(percent = Freq / sum(Freq) * 100)

df_plot$Var1<-factor(df_plot$Var1, levels = c("non-mucus secreating","mucus secreting" ))

ducts_mucus_count<-ggplot(df_plot, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill", color="black") +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Percent of bins")+theme_bw()+scale_fill_manual(values=c("#cccc16","#eb150c"))
 # geom_text(data = totals,  aes(x = Var2, y = 1.05, label = label),  inherit.aes = FALSE)
ducts_mucus_count
save_plts(ducts_mucus_count, "HLiCA_mucus_duct_size", w=5,h=6)

tab <- table(plt_cholangiocyte_spatial$type, plt_cholangiocyte_spatial$bile_duct_size_by_diameter)
chisq.test(tab)


#### pie plot
df <- plt_cholangiocyte_spatial %>%
  count(bile_duct_size_by_diameter, type) %>%
  tidyr::pivot_wider(names_from = type, values_from = n, values_fill = 0)

df_noNA<-df[which(!(is.na(df$bile_duct_size_by_diameter))),]



mean_diameter_by_size
duct_counts <- plt_cholangiocyte_spatial %>%
  distinct(bile_duct_cluster, bile_duct_size_by_diameter) %>%
  count(bile_duct_size_by_diameter)
duct_counts

df_noNA$x <- c(1,1.5,3)
df_noNA$y <- 1
df_noNA$r <- c(14/147, 54/147, 1.0)  # biologically defined relative sizes

muscus_scatter<-ggplot() +
  geom_scatterpie(data = df_noNA, aes(x = x, y = y, r = r), cols = c("non-mucus secreating","mucus secreting")) +
  coord_equal() +
  scale_x_continuous( breaks = df$x,labels = df$bile_duct_size_by_diameter) +
  theme_classic() +scale_fill_manual(values=c("#cccc16","#eb150c"))+
  labs(x = "Bile duct size", y = NULL)
save_plts(muscus_scatter, "HLiCA_mucus_piechart_duct_size", w=6,h=4)





# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# head(cho_meta)
# 
# 
# 
# #### as continuous
# percent_df <- cho_meta %>%
#   count(bile_duct_cluster, type) %>% 
#   group_by(bile_duct_cluster) %>%
#   mutate(percent = n / sum(n) * 100)
# 
# continuous_size_ducts<-read.csv(file=here("data/metadata_for_C107_bileductsize.csv"))
# bin_count <- as.data.frame(table(continuous_size_ducts$cholangiocyte_cluster))
# colnames(bin_count)<-c("bile_duct_cluster","size_bile_duct_cluster")
# 
# percent_df_size<-merge(percent_df, bin_count, by="bile_duct_cluster")
# 
# percent_df_size_secreting<-percent_df_size[which(percent_df_size$type=="mucus secreting"),]
# ggplot(percent_df_size_secreting, aes(size_bile_duct_cluster, percent, color=type))+geom_point()+stat_smooth(method="lm")
# 
# 
# ## convert to 8um^2
# cho_meta<-merge(cho_meta, bin_count, by="bile_duct_cluster")
# cho_meta$um_size<-cho_meta$size_bile_duct_cluster*64
# 
# levels(df_plot$Var2)<-c("Large duct\n(>1280µm²)","Small duct\n(<1280µm²)", "tehcnical")
# ggplot(df_plot, aes(x = Var2, y = prop)) +
#   geom_col(width = 0.6, fill = "steelblue") +
#   scale_y_continuous(labels = scales::percent) +
#   labs(x = "", y = "Secreting cholangiocytes (%)") +
#   theme_classic()
# 
# 
# ## counts for ducts non-technical
# duct_count<-cho_meta[,c("bile_duct_size","bile_duct_cluster")]
# duct_count<-duct_count[!duplicated(duct_count),]
# 
# totals <- df_per %>%
#   group_by(Var2) %>%
#   summarise(total = sum(Freq))
# 
# totals$ducts<-c(30, 218)
# totals$label<-paste(totals$ducts, " ducts\n(",totals$total,") bins",sep="" )
# 
# df_per$Var1 <- factor(df_per$Var1, levels = c("non-mucus secreating", "mucus secreting"))
# levels(df_per$Var2)<-c("Large duct\n(>1280µm²)","Small duct\n(<1280µm²)", "tehcnical")
# levels(totals$Var2)<-c("Large duct\n(>1280µm²)","Small duct\n(<1280µm²)", "tehcnical")
# 
# ducts_mucus_count<-ggplot(df_per, aes(x = Var2, y = Freq, fill = Var1)) +
#   geom_bar(stat = "identity", position = "fill", color="black") +
#   scale_y_continuous(labels = scales::percent) +
#   labs(y = "Percent of bins")+theme_bw()+scale_fill_manual(values=c("#cccc16","#eb150c"))+
#   geom_text(data = totals,  aes(x = Var2, y = 1.05, label = label),  inherit.aes = FALSE) 
# save_plts(ducts_mucus_count, "HLiCA_mucus_duct_size", w=5,h=6)
# 
