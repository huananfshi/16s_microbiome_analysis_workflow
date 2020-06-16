#loading required packages
library(vegan)
library(ALDEx2)
library(ComplexHeatmap)
library(RColorBrewer)
metadata=read.csv('~/.../metadata_all.csv',row.names = 1)
asv=read.csv('~/.../asv_grouped.csv',header = T,row.names = 1)
# ALDEx2 clr, with MC =16, using Kruskal–Wallis with BH FDR correction
samples_30 = rownames(metadata[metadata$Time=='30',])
asv_run_30= asv[,samples_30]
con_30 = metadata[samples_30,]$Y.A
names(con_30)=samples_30
asv_run.clr_30 = aldex.clr(asv_run_30,con_30,mc.samples = 16,denom="all", verbose=FALSE)
asv_run.kw_30 = aldex.kw(asv_run.clr_30)
asv_run.effect_30=aldex.effect(asv_run.clr_30,verbose = F)
asv_run.all_30=data.frame(asv_run.kw_30,asv_run.effect_30)
# par(mfrow=c(1,2))
# aldex.plot(asv_run.all_30,type = 'MA',test="kw",xlab="Log-ratio abundance",ylab="Difference")
# aldex.plot(asv_run.all_30, type="MW", test="kw", xlab="Dispersion",ylab="Difference")
# dev.off()
write.csv(asv_run.kw_30,file = '~/.../30days_asv_stats.csv')

#reformat clr transformed taxa table for beta diversity: Aitchson distance and PCA visualization
#using the median of 16 MC instances of each sample
asv_30.clr = list()
for(i in samples_30){
  asv_30.clr[[i]]=apply(data.frame(asv_run.clr_30@analysisData[[i]]),1,median)
}
asv_30.clr = t(data.frame(asv_30.clr))
asv_30.dist=vegdist(asv_30.clr, method = "euclidean")
asv_30.stat = adonis(asv_30.dist~con_30,method = "euclidean")
PCA_30 = betadisper(asv_30.dist,con_30)

#calculate taxa proportion/relative abundance for alpha diversity
asv_30_prop={}
for(i in samples_30){
  asv_30_prop[[i]]= asv_run_30[i]/apply(asv_run_30,2,sum)[i]

}
asv_30_prop = t(data.frame(asv_30_prop))
asv_30_prop = cbind(as.character(metadata[samples_30,]$Y.A),asv_30_prop)
asv_30_prop = data.frame(asv_30_prop)

shannon_30=diversity(t(asv_run_30), index = "shannon")
shannon_30.t=t.test(shannon_30~metadata[samples_30,]$Y.A)
simpson_30=diversity(t(asv_run_30), index = "simpson")
simpson_30.t=t.test(simpson_30~metadata[samples_30,]$Y.A)
chao1_30=estimateR(t(asv_run_30),index='chao')['S.chao1',]
chao1_30.t=t.test(chao1_30~metadata[samples_30,]$Y.A)
S_30 = apply(t(asv_run_30)[,-1]>0,1,sum)
evenness_30=shannon_30/log(S_30)
evenness_30.t=t.test(evenness_30~metadata[samples_30,]$Y.A)

#60 days and 90 days are the same as 30 days
samples_60 = rownames(metadata[metadata$Time=='60',])
asv_run_60= asv[,samples_60]
con_60 = metadata[samples_60,]$Y.A
names(con_60)=samples_60
asv_run.clr_60 = aldex.clr(asv_run_60,con_60,mc.samples = 16,denom="all", verbose=FALSE)
asv_run.kw_60 = aldex.kw(asv_run.clr_60)
asv_run.effect_60=aldex.effect(asv_run.clr_60,verbose = F)
asv_run.all_60=data.frame(asv_run.kw_60,asv_run.effect_60)
# par(mfrow=c(1,2))
# aldex.plot(asv_run.all_60,type = 'MA',test="kw",xlab="Log-ratio abundance",ylab="Difference")
# aldex.plot(asv_run.all_60, type="MW", test="kw", xlab="Dispersion",ylab="Difference")
#dev.off()
write.csv(asv_run.kw_60,file = '~/.../60days_asv_stats.csv')

asv_60.clr = list()
for(i in samples_60){
  asv_60.clr[[i]]=apply(data.frame(asv_run.clr_60@analysisData[[i]]),1,median)
}
asv_60.clr = t(data.frame(asv_60.clr))
asv_60.dist=vegdist(asv_60.clr, method = "euclidean")
asv_60.stat = adonis(asv_60.dist~con_60,method = "euclidean")
PCA_60 = betadisper(asv_60.dist,con_60)

asv_60_prop={}
for(i in samples_60){
  asv_60_prop[[i]]= asv_run_60[i]/apply(asv_run_60,2,sum)[i]

}
asv_60_prop = t(data.frame(asv_60_prop))
asv_60_prop = cbind(as.character(metadata[samples_60,]$Y.A),asv_60_prop)
asv_60_prop = data.frame(asv_60_prop)

shannon_60=diversity(t(asv_run_60), index = "shannon")
shannon_60.t=t.test(shannon_60~metadata[samples_60,]$Y.A)
simpson_60=diversity(t(asv_run_60), index = "simpson")
simpson_60.t=t.test(simpson_60~metadata[samples_60,]$Y.A)
chao1_60=estimateR(t(asv_run_60),index='chao')['S.chao1',]
chao1_60.t=t.test(chao1_60~metadata[samples_60,]$Y.A)
S_60 = apply(t(asv_run_60)[,-1]>0,1,sum)
evenness_60=shannon_60/log(S_60)
evenness_60.t=t.test(evenness_60~metadata[samples_60,]$Y.A)


samples_90 = rownames(metadata[metadata$Time=='90'&metadata$Location=='Fecal',])
asv_run_90= asv[,samples_90]
con_90 = metadata[samples_90,]$Y.A
names(con_90)=samples_90
asv_run.clr_90 = aldex.clr(asv_run_90,con_90,mc.samples = 16,denom="all", verbose=FALSE)
asv_run.kw_90 = aldex.kw(asv_run.clr_90)
asv_run.effect_90=aldex.effect(asv_run.clr_90,verbose = F)
asv_run.all_90=data.frame(asv_run.kw_90,asv_run.effect_90)
# par(mfrow=c(1,2))
# aldex.plot(asv_run.all_90,type = 'MA',test="kw",xlab="Log-ratio abundance",ylab="Difference")
# aldex.plot(asv_run.all_90, type="MW", test="kw", xlab="Dispersion",ylab="Difference")
#dev.off()
write.csv(asv_run.kw_90,file = '~/.../90days_asv_stats.csv')

asv_90.clr = list()
for(i in samples_90){
  asv_90.clr[[i]]=apply(data.frame(asv_run.clr_90@analysisData[[i]]),1,median)
}
asv_90.clr = t(data.frame(asv_90.clr))
asv_90.dist=vegdist(asv_90.clr, method = "euclidean")
asv_90.stat = adonis(asv_90.dist~con_90,method = "euclidean")
PCA_90 = betadisper(asv_90.dist,con_90)
asv_90_prop={}
for(i in samples_90){
  asv_90_prop[[i]]= asv_run_90[i]/apply(asv_run_90,2,sum)[i]

}
asv_90_prop = t(data.frame(asv_90_prop))
asv_90_prop = cbind(as.character(metadata[samples_90,]$Y.A),asv_90_prop)
asv_90_prop = data.frame(asv_90_prop)

shannon_90=diversity(t(asv_run_90), index = "shannon")
shannon_90.t=t.test(shannon_90~metadata[samples_90,]$Y.A)
simpson_90=diversity(t(asv_run_90), index = "simpson")
simpson_90.t=t.test(simpson_90~metadata[samples_90,]$Y.A)
chao1_90=estimateR(t(asv_run_90),index='chao')['S.chao1',]
chao1_90.t=t.test(chao1_90~metadata[samples_90,]$Y.A)
S_90 = apply(t(asv_run_90)[,-1]>0,1,sum)
evenness_90=shannon_90/log(S_90)
evenness_90.t=t.test(evenness_90~metadata[samples_90,]$Y.A)


# a proportion summary table for stacked bar chart (generated in GraphPad)
asv_prop = rbind(asv_30_prop,asv_60_prop,asv_90_prop)
asv_taxa = colnames(asv_prop)
# get genus information only
asv_taxa=strsplit(asv_taxa,";")
asv_taxa=unique(asv_taxa)
for(i in 1:length(asv_taxa)){
  asv_taxa[[i]]=unique(asv_taxa[[i]])
}

asv_genus=list()
for(i in 1:length(asv_taxa)){
  l=length(asv_taxa[[i]])
  asv_genus[[i]] = paste(asv_taxa[[i]][l-1],"_",asv_taxa[[i]][l])
}
asv_genus=unlist(asv_genus)
asv_genus=gsub(' _ ',"_",asv_genus)
asv_genus[1]='Y_A'
colnames(asv_prop)=asv_genus
write.csv(asv_prop,"~/.../asv_prop.csv")

# ploting alpha diversity
# !!!! check the level of metadata i.e. levels(metadata[samples_30,]$Y.A) to determine the x axis labels order :names =levels(metadata[samples_30,]$Y.A)

par(cex.lab=1.5,cex.axis=1.5,mfrow=c(1,3))
boxplot(shannon_30~metadata[samples_30,]$Y.A, xlab = "", ylab = "",ylim=c(1.5,2.6),outline=F,names=c('Aged','Young'))
points(shannon_30~metadata[samples_30,]$Y.A)
mtext(paste('p=',as.character(round(shannon_30.t$p.value,5))),cex=1)
title('Shannon')

boxplot(simpson_30~metadata[samples_30,]$Y.A, xlab = "", ylab = "",ylim=c(0.5,1),outline=F,names=c('Aged','Young'))
points(simpson_30~metadata[samples_30,]$Y.A)
mtext(paste('p=',as.character(round(simpson_30.t$p.value,4))),cex=1)
title('Simpson')

boxplot(chao1_30~metadata[samples_30,]$Y.A, xlab = "", ylab = "",ylim=c(15,50),outline=F,names=c('Aged','Young'))
points(chao1_30~metadata[samples_30,]$Y.A)
mtext(paste('p=',as.character(round(chao1_30.t$p.value,4))),cex=1)
title('Chao1')


par(cex.lab=1.5,cex.axis=1.5,mfrow=c(1,3))
boxplot(shannon_60~metadata[samples_60,]$Y.A, xlab = "", ylab = "",ylim=c(0.5,2.8),outline=F,names=c('Aged','Young'))
points(shannon_60~metadata[samples_60,]$Y.A)
mtext(paste('p=',as.character(round(shannon_60.t$p.value,4))),cex=1)
title('Shannon')

boxplot(simpson_60~metadata[samples_60,]$Y.A, xlab = "", ylab = "",ylim=c(0,1),outline=F,names=c('Aged','Young'))
points(simpson_60~metadata[samples_60,]$Y.A)
mtext(paste('p=',as.character(round(simpson_60.t$p.value,4))),cex=1)
title('Simpson')

boxplot(chao1_60~metadata[samples_60,]$Y.A, xlab = "", ylab = "",ylim=c(20,55),outline=F,names=c('Aged','Young'))
points(chao1_60~metadata[samples_60,]$Y.A)
mtext(paste('p=',as.character(round(chao1_60.t$p.value,4))),cex=1)
title('Chao1')

par(cex.lab=1.5,cex.axis=1.5,mfrow=c(1,3))
boxplot(shannon_90~metadata[samples_90,]$Y.A, xlab = "", ylab = "",ylim=c(1.6,2.6),outline=F,names=c('Aged','Young'))
points(shannon_90~metadata[samples_90,]$Y.A)
mtext(paste('p=',as.character(round(shannon_90.t$p.value,5))),cex=1)
title('Shannon')

boxplot(simpson_90~metadata[samples_90,]$Y.A, xlab = "", ylab = "",ylim=c(0.6,1),outline=F,names=c('Aged','Young'))
points(simpson_90~metadata[samples_90,]$Y.A)
mtext(paste('p=',as.character(round(simpson_90.t$p.value,4))),cex=1)
title('Simpson')

boxplot(chao1_90~metadata[samples_90,]$Y.A, xlab = "", ylab = "",ylim=c(25,55),outline=F,names=c('Aged','Young'))
points(chao1_90~metadata[samples_90,]$Y.A)
mtext(paste('p=',as.character(round(chao1_90.t$p.value,4))),cex=1)
title('Chao1')


# plotting PCAs

par(cex.lab=1.5,cex.axis=1.5,mfrow=c(1,1))
pca.30.fig = ordiplot(PCA_30, type = "none",
                      xlab = paste("PC 1 (", round(PCA_30$eig[1]/sum(PCA_30$eig)*100,2), "% explained)"),
                      ylab = paste("PC 2 (", round(PCA_30$eig[2]/sum(PCA_30$eig)*100,2), "% explained)"), xlim = c(-20,40))

points(pca.30.fig, "sites", pch = 17, col = "black",  select =con_30 == "A")
points(pca.30.fig, "sites", pch = 19, col = "black", select = con_30 == "Y")
ordispider(PCA_30,con_30,col= c("black","black"))
legend("topright",legend = c('Aged','Young'), pch = c(17,19),bty='n')
mtext(paste("p=", as.character(asv_30.stat$aov.tab$`Pr(>F)`[1]),'\n',"F=",as.character(round(asv_30.stat$aov.tab$F.Model[1],2))),cex=1.5,at=30)
title('30 days post FTG',cex=1.5)

par(cex.lab=1.5,cex.axis=1.5,mfrow=c(1,1))
pca.60.fig = ordiplot(PCA_60, type = "none",
                      xlab = paste("PC 1 (", round(PCA_60$eig[1]/sum(PCA_60$eig)*100,2), "% explained)"),
                      ylab = paste("PC 2 (", round(PCA_60$eig[2]/sum(PCA_60$eig)*100,2), "% explained)"), xlim = c(-20,25))

points(pca.60.fig, "sites", pch = 17, col = "black",  select =con_60 == "A")
points(pca.60.fig, "sites", pch = 19, col = "black", select = con_60 == "Y")
ordispider(PCA_60,con_60,col= c("black","black"))
mtext(paste("p=", as.character(asv_60.stat$aov.tab$`Pr(>F)`[1]),'\n',"F=",as.character(round(asv_60.stat$aov.tab$F.Model[1],2))),cex=1.5,at=30)
title('60 days post FTG',cex=1.5)

par(cex.lab=1.5,cex.axis=1.5,mfrow=c(1,1))
pca.90.fig = ordiplot(PCA_90, type = "none",
                      xlab = paste("PC 1 (", round(PCA_90$eig[1]/sum(PCA_90$eig)*100,2), "% explained)"),
                      ylab = paste("PC 2 (", round(PCA_90$eig[2]/sum(PCA_90$eig)*100,2), "% explained)"), xlim = c(-20,30))

points(pca.90.fig, "sites", pch = 17, col = "black",  select =con_90 == "A")
points(pca.90.fig, "sites", pch = 19, col = "black", select = con_90 == "Y")
ordispider(PCA_90,con_90,col= c("black","black"))
mtext(paste("p=", as.character(asv_90.stat$aov.tab$`Pr(>F)`[1]),'\n',"F=",as.character(round(asv_90.stat$aov.tab$F.Model[1],2))),cex=1.5,at=30)
title('90 days post FTG',cex=1.5)
dev.off()

# generating effect size heatmap
genus_30 = rownames(asv_run.all_30)[asv_run.all_30$kw.eBH<=0.2]
genus_60 = rownames(asv_run.all_60)[asv_run.all_60$kw.eBH<=0.2]
genus_90 = rownames(asv_run.all_90)[asv_run.all_90$kw.eBH<=0.2]

genus=c(genus_30,genus_60,genus_90)
genus =unique(genus)

asv_30.filter=asv_run.all_30[genus,c('kw.eBH','effect')]
asv_60.filter=asv_run.all_60[genus,c('kw.eBH','effect')]
asv_90.filter=asv_run.all_90[genus,c('kw.eBH','effect')]
rownames(asv_90.filter) = genus

asv_30.filter[is.na(asv_30.filter)]=0
asv_60.filter[is.na(asv_60.filter)]=0
asv_90.filter[is.na(asv_90.filter)]=0

genus=strsplit(genus,";")
genus=unique(genus)
for(i in 1:length(genus)){
  genus[[i]]=unique(genus[[i]])
}

genus_sample=list()
for(i in 1:length(genus)){
  l=length(genus[[i]])
  genus_sample[[i]] = paste(genus[[i]][l-1],"_",genus[[i]][l])
}
genus_sample=unlist(genus_sample)
genus_sample=gsub(' _ ',"_",genus_sample)
asv = data.frame(asv_30.filter$effect,asv_60.filter$effect,asv_90.filter$effect)
colnames(asv)=c("30 days","60 days","90 days")
rownames(asv)=genus_sample

asv_mat = as.matrix(asv)

H_map = Heatmap(asv_mat,
                col=colorRampPalette(brewer.pal(8, "RdBu"))(25),
                show_row_dend = F,
                show_column_dend = FALSE,
                show_column_names=T,
                column_names_centered = T,
                column_names_rot = 0,
                column_names_side = "top",
                show_row_names=TRUE,
                row_names_side = "left",
                column_names_gp = gpar(fontsize = 10),
                row_names_gp = gpar(fontsize = 10),
                column_title="Effect of Young Microbiome",
                column_title_side = "top",
                column_title_gp = gpar(fontsize = 14),
                cluster_rows = F,
                cluster_columns = F,
                show_heatmap_legend = T,
                heatmap_legend_param = list(at = c(-3,0,3),
                                            labels = c("-3","0","3"),
                                            color_bar = "continuous",
                                            legend_direction="vertical",
                                            labels_gp = gpar(fontsize = 8),
                                            legend_height= unit(4, "cm"),
                                            title_position = "leftcenter-rot",
                                            title = "Effect"
                ),
                width = unit(7,'cm')
)

draw(H_map,
     heatmap_legend_side = "right",
     row_title="Genus",
     row_title_side="left",
     column_title = "",
     column_title_side = "bottom",
     column_title_gp = gpar(fontsize = 8),
     row_title_gp = gpar(fontsize = 15),

)
