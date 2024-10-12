oft = read.csv("MIA2_zscore.csv", row.names = 1)
#oft = na.omit(oft)
control = oft[oft$sex == "m" & oft$group == "WT", ]
oft$group = factor(oft$group, labels = label_group)
oft$Ztc = (mean(na.omit(control$Time_in_center_100))-oft$Time_in_center_100)/sd(na.omit(control$Time_in_center_100))
oft$Zdr = (mean(na.omit(control$Distance_in_center_100))-oft$Distance_in_center_100)/sd(na.omit(control$Distance_in_center_100))
oft$Zof = (oft$Ztc+oft$Zdr)/2
oft$Ztoa = (mean(na.omit(control$Time_in_open_100))-oft$Time_in_open_100)/sd(na.omit(control$Time_in_open_100))
oft$Zdor = (mean(na.omit(control$Distance_in_open_100))-oft$Distance_in_open_100)/sd(na.omit(control$Distance_in_open_100))
oft$Zepm = (oft$Ztoa+oft$Zdor)/2

oft$Ztst = (oft$TST_immobile_time_100-mean(control$TST_immobile_time_100))/sd(control$TST_immobile_time_100)

control_ppi = oft[oft$sex == "m" & oft$group == "WT", ]
control_ppi = control_ppi[complete.cases(control_ppi$SI), ]
oft$Ztcst = (mean(control_ppi$SI)-oft$SI)/sd(control_ppi$SI)


control_ppi = oft[oft$sex == "m" & oft$group == "WT", ]
control_ppi = control_ppi[complete.cases(control_ppi$n74, control_ppi$n78, control_ppi$n82, control_ppi$n86, control_ppi$n90), ]
oft$Z78 = (mean(na.omit(control_ppi$n78))-oft$n78)/sd(na.omit(control_ppi$n78))
oft$Z82 = (mean(na.omit(control_ppi$n82))-oft$n82)/sd(na.omit(control_ppi$n82))
oft$Z86 = (mean(na.omit(control_ppi$n86))-oft$n86)/sd(na.omit(control_ppi$n86))
oft$Zppi = (oft$Z78+oft$Z82+oft$Z86)/3
write.csv(oft, "MIA2_zscore.csv")
