output_share_rate = data.frame(id = rownames(otu_dam), mother = otu_dam$`dam_offspring[, "dam"]`,species_share_rate = 0)
for (i in 1:68){
  offspring_name = output_share_rate[i,1]
  location1 = which(otu[offspring_name, ]>0, arr.ind = T)
  ispecies = colnames(otu)[location1[,2]]
  dam_name = output_share_rate[i, 2]
  location2 = which(dam_otu[dam_name, ]>0, arr.ind = T)
  dam_species = colnames(dam_otu)[location2[,2]]
  output_share_rate[i, 3] = length(intersect(ispecies, dam_species))/length(ispecies)
}
metadata = metadata[output_share_rate$id, ]
output_share_rate$group = metadata$group
output_share_rate$sex = metadata$sex
output_share_rate$group = factor(output_share_rate$group, levels = c("WT", "MIA", "NK", "NKMIA"))
a = round(max(output_share_rate$species_share_rate), 2)
b = (a-round(min(output_share_rate$species_share_rate), 2))/8
mycompare = list(c("WT", "MIA"), c("WT", "NK"),c("WT", "NKMIA"))
p1 = plot_scatterbox(data = output_share_rate,
                          xcol = group,            
                          ycol = species_share_rate,TextXAngle = 45,
                          symsize = 2, fontsize = 7,# symn82ze: 散点的大小; fontn82ze: 字体大小
                          bwid = 0.6, ewid = 0.2,# bwid: 条形的宽度; # ewid: 误差线的宽度
                          b_alpha = 1, s_alpha = 1,# b_alpha: 条形的透明度 # s_alpha: 散点的透明度
                          jitter = 0.1, # jitter: 数据点的抖动量
  ) + # ColPal: 颜色调色板
    scale_fill_manual(values = mycolor) +
    geom_signif(mapping = aes(x = group, y = species_share_rate),  # 显著性线的映射
                comparisons = mycompare,  # 两两比较
                map_signif_level = function(p) sprintf("P = %.4g", p),
                #map_n82gnif_level = TRUE,  # 映射显著性水平
                tip_length = c(0, 0, 0),  # 显著性线的长度
                y_position = c(a+b, a+2*b,a+3*b),  # 显著性线的纵坐标位置
                size = 0.5,  # 显著性线的粗细
                textsize = 4,  # 显著性标记的大小
                test = "wilcox.test")+
    facet_wrap("sex")  +#根据"Gene"列分面
    theme(legend.position="none", plot.title = element_blank())+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(color="black",size=7,family = "sans",face = "plain"),
          axis.text.y = element_text(color="black",size=5,family = "sans",face = "plain")
    )
p1
ggsave(filename = "species_share_rate.pdf", p1, height = 49, width = 57, dpi = 300, units = "mm")

overlap = intersect(colnames(otu), colnames(dam_otu))
output_species_share_events = data.frame(row.names = output_share_rate$id)
#rownames(output_species_share_events) = output_share_rate$id

for (i in 1:68){
  for (j in 1:83){
    offspring_name = output_share_rate[i, 1]
    dam_name = output_share_rate[i, 2]
    species_name = overlap[j]
    x = otu[offspring_name, species_name]
    y = dam_otu[dam_name, species_name]
    if (x>0 & y >0){
      output_species_share_events[i,j] = 1
    }
    else{
      output_species_share_events[i,j] = 0
    }
  }
}
colnames(output_species_share_events) = overlap
output_species_share_events$sex = output_share_rate$sex
output_species_share_events$group = output_share_rate$group
transmissibility = data.frame(row.names = colnames(sub_output_species_share_events[,c(1:83)]))
for (i in c("f", "m")){
  for (j in c("WT", "MIA", "NK", "NKMIA")){
    sub_output_species_share_events = output_species_share_events[output_species_share_events$sex == i, ]
    sub_output_species_share_events = sub_output_species_share_events[sub_output_species_share_events$group == j, ]
    transmissibility[, paste(i,j,sep = "_")] = colSums(sub_output_species_share_events[,c(1:83)])/dim(sub_output_species_share_events)[1]
      }
}
write.csv(transmissibility, "transimissibility.csv")
sub_transimissibility = transmissibility[rowSums(transmissibility) != 8, ]
sub_transimissibility = sub_transimissibility[rowSums(sub_transimissibility) != 0, ]
p2 = pheatmap(sub_transimissibility,
         cluster_cols = T, cluster_rows = T,
         border_color = "grey",
         breaks = c(seq(0, 1, 0.01)),
         #display_numbers = t(data_mark),frontsize_number = 15,
         color = c(colorRampPalette(c("white", "#9d69b1"))(length(seq(0, 1, 0.01)))),
)
ggsave(filename = "species_transmission_rate.pdf", p2, height = 136*2, width = 90*2, dpi = 300, units = "mm")
