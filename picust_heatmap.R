library(ComplexHeatmap)
enzyme30 = read.csv("~/.../30days_enzymes.csv",header = T,row.names = 1)
enzyme30_matrix = as.matrix(enzyme30)
enzyme60 = read.csv("~/...s/60days_enzymes.csv",header = T,row.names = 1)
enzyme60_matrix = as.matrix(enzyme60)
enzyme90 = read.csv("~/.../90days_enzymes.csv",header = T,row.names = 1)
enzyme90_matrix = as.matrix(enzyme90)
enzyme = cbind(enzyme30_matrix,enzyme60_matrix,enzyme90_matrix)
colnames(enzyme)=c('30_Y','30_A','60_Y','60_A','90_Y','90_A')
Heatmap(enzyme,
        column_split = c('30_Y','30_A','60_Y','60_A','90_Y','90_A'),
        column_gap = unit(c(0,2,0,2,0),'mm'),
        column_title=c('30 days',"",'60 days',"",'90 days',""),
        show_row_dend = F,
        show_column_dend = F,
        cluster_columns = F,
        row_names_side = "left",
        show_column_names =T,
        column_names_side = "top",
        column_names_rot = 0,
        heatmap_legend_param = list(at = c(-2,2),
                                    labels = c("low","high"),
                                    color_bar = "continuous",
                                    legend_direction="vertical",
                                    labels_gp = gpar(fontsize = 8),
                                    legend_height= unit(4, "cm"),
                                    title_position = "leftcenter-rot",
                                    title = "Relative abundance"
        ),
        width = unit(10,'cm'))
dev.off()
