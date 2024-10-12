#ggclusterNet
library(BiocManager)
install("tidyfst")
install("igraph")
install("sna")
install("phyloseq")
install("ggalluvial")
install("ggraph")#error
install("WGCNA")
install("ggnewscale")
install("pulsar")
install("patchwork")
remotes::install_github("taowenmicro/EasyStat")
remotes::install_github("taowenmicro/ggClusterNet")
library(tidyverse)
library("tidyfst")
library("igraph")
library("sna")
library("phyloseq")
library("ggalluvial")
library("ggraph")#error
library("WGCNA")
library("ggnewscale")
library("pulsar")
library("patchwork")
library(ggClusterNet)
library(xlsx)
metadata = read.csv("MIA2_zscore.csv", row.names = 1)
otutab = read.delim("merged_abundance_table_genus.txt", row.names = 1)
metadata = metadata[metadata$SID %in% colnames(otutab), ]
rownames(metadata) = metadata$SID
result = cor_Big_micro2(otutab)

taxonomy = read.xlsx("taxonomy.xlsx", sheetIndex = 1, header = T)
taxonomy = taxonomy[!duplicated(taxonomy$Genus), ]
rownames(taxonomy)= taxonomy$Species
rownames(taxonomy)= gsub("g__", "",taxonomy$Genus)
taxonomy = taxonomy[rownames(otutab), ]
#subgroup network
sub_metadata1 = metadata[metadata$sex == "f" & metadata$group == "WT", ]
sub_otutab = otutab[,colnames(otutab) %in% rownames(sub_metadata1)]
ps = phyloseq(sample_data(sub_metadata1), 
              otu_table(as.matrix(sub_otutab), taxa_are_rows = T),
              tax_table((as.matrix(taxonomy))))##tax and otutab's rownames are the same
result = corMicro(ps = ps, N = 50, 
                  r.threshold = 0.4,
                  p.threshold = 0.05,
                  method = "spearman")
cor = result[[1]]
igraph = graph_from_adjacency_matrix(cor, diag = F, mode = "undirected", weighted = T)
dat = net_properties.4(igraph, n.hub = F)
nodepro = node_properties(igraph)
nodepro = as.data.frame(nodepro)
#write.csv(dat, "ggclusternet/fNK_net.txt")
#write.csv(nodepro, "ggclusternet/fNK_node.csv")
ps_net = result[[3]]
otu_table = ps_net %>%
  vegan_otu() %>%
  t() %>%
  as.data.frame()
#gp = data.frame(ID = rownames(otu_table), group = sample(1:3, 92, replace = T))
#layout = PolygonClusterG(cor = cor, nodeGroup = gp)
layout = model_maptree2(cor = cor, method = "cluster_fast_greedy", seed = 12)
node = layout[[1]]
tax_table = ps_net %>%
  vegan_tax()%>%
  as.data.frame()
nodes = nodeadd(plotcord = node, otu_table = otu_table, tax_table = tax_table)
edge = edgeBuild(cor = cor, node = node)
write.csv(edge, "ggclusternet/edge_fwt.csv")
#nodes$level = sample(c("enriched", "depleted","others"), 92, replace = T)
grp2 = read.xlsx("difference.p0.05.genus_D  A.xlsx", sheetIndex = 1)
grp2$enrichment = ifelse(grp2$enrichment == "D<A", "enriched", "depleted")
rownames(grp2) = grp2$Row.names
nodes$level = "others"
for (i in 1:50){
  if (rownames(nodes)[i] %in% grp2$Row.names){
    nodes[i, 'level'] = grp2[rownames(nodes)[i], "enrichment"]
  }
}
nodes$closeness = 0
for (i in 1:50){
  if (rownames(nodes)[i] %in% rownames(nodepro)){
    nodes[i, 'closeness'] = nodepro[rownames(nodes)[i], "igraph.closeness"]
  }
}
#nodes$showpoint = ifelse(nodes$level == 'others', FALSE,TRUE)
# 开始绘图
p1 <- ggplot() + 
  geom_segment(data = edge, aes(x = X1, y = Y1, xend = X2, yend = Y2,color = cor), size = 0.4) +
  geom_point(data = nodes, shape = 21, colour = "black", stroke = 0.4,  aes(X1, X2, size = mean,fill = level))+
  scale_fill_manual(values = c('enriched' = "#c62320", 'depleted'= "#132b69", 'others' = "white"))+
  scale_color_manual(values = c("+" = "pink","-" = "#9ab8d4"))+
  scale_size_continuous(range = c(2, 5)) + 
  scale_x_continuous(breaks = c(-6,-3,0,3,6)) + 
  scale_y_continuous(breaks = c(-6,-3,0,3,6)) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA), legend.text = element_text(size = 7)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p1
ggsave("ggclusternet/fNKMIA.pdf", p1, width = 50*3, height = 32*3, units = 'mm') # 保存图片
nodes_mwt = nodes
write.csv(nodes, "ggclusternet/nodes_fnkmia.csv")
dif1 = nodesnodes_mwtdif1 = nodes[nodes$level == 'enriched', ]
nodepro[rownames(dif1), ]

mwt_hub = c("Erysipelotrichaceae_unclassified")
mnkmia_hub = c("Anaerotruncus")
fwt_hub = c("Pseudoflavonifractor")
fnkmia_hub = c("Faecalibaculum")



##closeness of key bacteria
all_hub = c(mwt_hub, mnkmia_hub, fwt_hub, fnkmia_hub)
nodepro1 = nodepro#mwt
nodepro2 = nodepro#mnkmia
nodepro3 = nodepro#fwt
nodepro4 = nodepro#fnkmia
mwt_closeness = as.data.frame(nodepro1)[all_hub,]
mnkmia_closeness = as.data.frame(nodepro2)[all_hub,]
fwt_closeness = as.data.frame(nodepro3)[all_hub, ]
fnkmia_closeness = as.data.frame(nodepro4)[all_hub, ]
all_closeness = cbind(mwt_closeness, mnkmia_closeness,fwt_closeness, fnkmia_closeness)
rownames(all_closeness) = all_hub
all_closeness = all_closeness[,c(2,6,10,14)]
all_closeness[is.na(all_closeness)] = 0
colnames(all_closeness) = c("mWT", "mNKMIA","fWT", "fNKMIA")
library(pheatmap)
p2 = pheatmap(all_closeness,scale = "column",
         cluster_cols = F, cluster_rows = T,
         border_color = "grey",
         #breaks = c(seq(0,0.3,0.1),seq(0.4,0.6,0.1)),
         #display_numbers = data_mark,frontsize_number = 15,
         color = c(colorRampPalette(c( "#6996e3","white"))(length(seq(-1,0,0.1))), colorRampPalette(c("white", "#d8527c"))(length(seq(0.1,1,0.1))))
)
p2
ggsave("ggclusternet/key_bac.pdf", p2, width = 41*3, height = 41*3, units = 'mm')
##closeness of all bacteria

overlap = intersect(rownames(nodepro1), rownames(nodepro2))
overlap = intersect(overlap, rownames(nodepro3))
overlap = intersect(overlap, rownames(nodepro4))
nodepro10 = as.data.frame(nodepro1)[overlap, ]
nodepro20 = as.data.frame(nodepro2)[overlap, ]
nodepro30 = as.data.frame(nodepro3)[overlap, ]
nodepro40 = as.data.frame(nodepro4)[overlap, ]
all_closeness = cbind(nodepro10, nodepro20,nodepro30, nodepro40)
all_closeness = all_closeness[,c(2,6,10,14)]
library(dplyr)
all_closeness0 <- mutate_all(all_closeness, funs(. * (!is.na(.) & . != 1)))
colnames(all_closeness0) = c("mWT", "mNKMIA","fWT", "fNKMIA")
all_closeness0 = as.data.frame(t(all_closeness0))
all_closeness0 = all_closeness0[,-c(47,48,50)]
library(pheatmap)
p3 = pheatmap(all_closeness0,scale = "column",
              cluster_cols = T, cluster_rows = F,cellwidth = 5,
              show_colnames = F,border_color = NA,
              #display_numbers = data_mark,frontsize_number = 15,
              color = c(colorRampPalette(c( "white", "#6996e3"))(length(seq(0,0.2,0.01))), colorRampPalette(c("#6996e3", "#d8527c"))(length(seq(0.3,0.7,0.001))))
)
p3
ggsave("ggclusternet/closeness.pdf", p3, width = 70*2, height = 27*2, units = 'mm')
##venn diagram of edge
edge1m = paste(edge$OTU_1, edge$OTU_2, sep = "_")#mwt
edge2m = paste(edge$OTU_1, edge$OTU_2, sep = "_")#mnkmia
edge3m = paste(edge$OTU_1, edge$OTU_2, sep = "_")#fwt
edge4m = paste(edge$OTU_1, edge$OTU_2, sep = "_")#fnkmia

length(intersect(edge1m, edge2m))
length(intersect(edge3m, edge4m))
length(intersect(edge1m, edge3m))
length(intersect(edge2m, edge4m))


edge1f = paste(edge$OTU_1, edge$OTU_2, sep = "_")#wt
edge2f = paste(edge$OTU_1, edge$OTU_2, sep = "_")#mia
edge3f = paste(edge$OTU_1, edge$OTU_2, sep = "_")#nk
edge4f = paste(edge$OTU_1, edge$OTU_2, sep = "_")#nkmia

length(intersect(edge1f, edge2f))
length(intersect(edge1f, edge3f))
length(intersect(edge1f, edge4f))

length(intersect(edge1f, edge1m))
length(intersect(edge2f, edge2m))
length(intersect(edge3f, edge3m))
length(intersect(edge4f, edge4m))


#鲁棒性：系统在网络结构发生扰动的情况下，外部干扰抵御能力的保持能力
res = Robustness.Random.removal(ps = ps_net, Top = 50)
res = Robustness.Targeted.removal(ps = ps_net, Top = 50)
#负相关比例：越大代表网络抵御扰动能力越强
res = negative.correlation.ratio(ps = ps)
#组成稳定性：评估两个网络之间变化多大
res = community.stability(ps = ps, time = F)
res = natural.con.micro(ps = ps, con.method = 'pulsar')
#网络抗毁性，网络中的节点发生自然失效或遭受故意攻击的条件下，网络维持其功能的能力

