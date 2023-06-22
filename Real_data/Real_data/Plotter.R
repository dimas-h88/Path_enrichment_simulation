library(ggplot2)
library(data.table)

sum_plot <- function(path, val, a) {
  table <- data.table(seed = numeric(),
                      Fraq = numeric(),
                      Method = character(),
                      N = character())
  
  files <- list.files(path = path, full.names = TRUE)
  for (file_path in files) {
    var <- fread(file_path)
    if (val == 'pval') {
      var <- var[, list(Fraq = mean(pval < 0.05, na.rm = TRUE)), by = .(seed, N / 2, Method)]
    } else if (val == 'padj') {
      var <- var[, list(Fraq = mean(padj < 0.05, na.rm = TRUE)), by = .(seed, N / 2, Method)]
    }
    table <- rbind(table, var)
  }
  if (val == 'pval') {
    yl <- 'Fraction'
    hline <- geom_hline(yintercept = 0.05, lty = 3)
  } else if (val == 'padj') {
    yl <- 'Fraction'
    hline <- NULL  
  }
  if (a == 'M0') {
    tit <- paste0(val, " < 0.05 in M0 vs M0 ")
  } else if (a == 'M1') {
    tit <- paste0(val, " < 0.05 in M0 vs M1 ")
  }
  table$N <- factor(table$N)
  table$N <- reorder(table$N, as.numeric(as.character(table$N))) # reorder N levels by increasing order
  plot <- ggplot(table, aes(x = Method, y = Fraq, fill = N)) +
    geom_boxplot() +
    theme_classic(base_size = 20) +
    hline +  # Add geom_hline conditionally
    labs(title = tit, y = yl) +
    ylim(0, 0.35)
  
  # Save the plot as a PNG image
  ggsave(paste0(a, "_", val, ".png"), plot = plot, width = 8, height = 6, dpi = 300)
}


sum_plot('./M1', 'pval', 'M1')


camera <- fread('./M0/summary_table_camera_m0.tsv')
fgsea <- fread('./M0/summary_table_fgsea_m0.tsv')
fgseaSE <- fread('./M0/summary_table_fgsea_SE_m0.tsv')
fora <- fread('./M0/summary_table_fora_m0.tsv')




par(mfrow = c(2, 4), mar = c(5, 4, 10, 2))

for (a in c(4, 40)) {
  hist(camera[camera$N == a,]$pval, breaks = 20, main = 'CAMERA', xlab = '', freq = FALSE)
  hist(fgsea[fgsea$N == a,]$pval, breaks = 20, main = 'fgsea', xlab = '', freq = FALSE)
  hist(fgseaSE[fgseaSE$N == a,]$pval, breaks = 20, main = 'fgsea_SE', xlab = '', freq = FALSE)
  hist(fora[fora$N == a,]$pval, breaks = 20, main = 'fora', xlab = '', freq = FALSE)
}

mtext("2 resamples in group", side = 3, line = -2.5, cex = 2.5, outer = TRUE)
mtext("20 resamples in group", side = 3, line = -32, cex = 2.5, outer = TRUE)



