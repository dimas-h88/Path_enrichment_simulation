camera <- read.csv('summary_table_camera_pval10.csv')
fgsea <- read.csv('summary_table_fgsea_pval10.csv')
fora <- read.csv('summary_table_fora_pval10.csv')
fgseaSE <- read.csv('summary_table_SE_pval10.csv')



par(mfrow = c(1, 5), mar = c(5, 4, 10, 2))

hist(camera[camera$k == 0.05,]$pval, breaks = 20, main = 'CAMERA', xlab ='')
hist(fgsea[fgsea$k == 0.05,]$pval, breaks = 20, main = 'fgsea', xlab ='')
hist(fgseaSE[fgseaSE$k == 0.05,]$pval, breaks = 20, main = 'fgsea_SE', xlab ='')
hist(fora[fora$k == 0.05,]$pval, breaks = 20, main = 'fora', xlab ='')

mtext("Distribution of p-values", side = 3, line = -2.5, cex = 2, outer = TRUE)




frac <- ggplot(summary_report, aes(x=method, y=Fraction_of_differentially_expressed_pathways, fill = factor(k) )) +  
  geom_boxplot() + theme_bw() + labs(title="% of pathways which are not diffexpressef, but have padj <0.05")

frac




