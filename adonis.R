#adonis analysis
otu = read.delim("merged_abundance_table_species.txt", row.names = 1)
metadata = read.csv("MIA2_zscore.csv", row.names = 1)
metadata = metadata[metadata$SID != 0, ]
rownames(metadata) = metadata$SID
metadata = metadata[rownames(otu), ]
library(vegan)
#install.packages("vegan")
library(ggpubr)
otu = as.data.frame(t(otu))
bray_dist = vegdist(otu, "bray")
bray_dist = as.matrix(bray_dist)
pcoa = cmdscale(bray_dist, k=3, eig=T)

#sex
sub_metadata = metadata[metadata$sex == "m", ]
sub_bray = bray_dist[rownames(sub_metadata), rownames(sub_metadata)]
pcoa = cmdscale(sub_bray, k=3, eig=T)
library(ggplot2)
# get coordinate string, format to dataframme
points = as.data.frame(pcoa$points) 
eig = pcoa$eig
points$group = sub_metadata$group
colnames(points) = c("PC1", "PC2", "PC3","group") 
points$names = rownames(points)
p = ggplot(points, aes(x=PC1, y=PC2, color=group)) + geom_point(alpha=.7, size=2) +
  #scale_color_manual(values=mycolors)+
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  stat_ellipse(level=0.95)+
  theme(plot.margin=unit(rep(0.5,4), 'lines'),
        panel.background=element_rect(fill='transparent', color='black'),
        panel.border=element_rect(fill='transparent', color='transparent'),
        panel.grid=element_blank())
  #geom_text(aes(label = names), fontface = "bold", hjust = .5, vjust = -.25) 
p
ado = adonis2(bray_dist ~ sex+genetype+treatment+genetype*treatment, metadata, permutations = 9999)#permanova
ado

library(ggplot2)
ado$new_name = rownames(ado)
#reorder(new_name, -LDA)
ado$label = ifelse(ado$`Pr(>F)`<0.05, "*", "ns")
ado = na.omit(ado)
ggplot(ado, aes(x = reorder(new_name, R2), y = R2) )+
  geom_bar(stat = "identity")+
  geom_text(aes(label=label), hjust = 0,color="red")+
  theme_classic()+
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  coord_flip()

