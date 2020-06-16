# an example of how p.adjust function works
statssss = read.csv("~/.../scfa_enzymes_taxa_stats.csv",header = T,sep = ',')
head(statssss)
colnames(statssss)

statssss['p_adjust_30']=p.adjust(statssss$A_30_Fecal_Y_30_Fecal,method = 'fdr')
statssss['p_adjust_60']=p.adjust(statssss$A_60_Fecal_Y_60_Fecal,method = 'fdr')
statssss['p_adjust_90']=p.adjust(statssss$A_90_Fecal_Y_90_Fecal,method = 'fdr')
write.csv(statssss,"~/.../scfa_enzymes_taxa_stats_adjusted.csv",row.names = F)
