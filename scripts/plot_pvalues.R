library(ggplot2)
library(plotly)

plot_pvalues <- function(filtered_table, p_title, column, num_pairs){
  # keep combinations with p-value under threshold
  
  p <- ggplot(data = filtered_table, mapping = aes(x = c(1:dim(filtered_table)[1]),
                                                   y = log(filtered_table[, column]))) +
    geom_count(size = 0.3) +
    xlab("Index") + 
    ylab("Log p-value") + 
    ggtitle(p_title, subtitle = "Only positions with p-value under threshhold") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    geom_hline(yintercept=log(0.05/num_pairs), color = "green")
  
  return(p)
}

p_values <- read.table("final_output.tsv", header = T, sep = '\t')

# plotting for pval1
filtered_table <- p_values[p_values[,8] < (0.05/dim(p_values)[1]),]
discarded <- length(filtered_table[,8][filtered_table[,8] == as.double(0)])
kept <- dim(filtered_table)[1] - discarded
title <- paste("Log pval1, ", kept, " from ", dim(filtered_table)[1])
p <- plot_pvalues(filtered_table, title, 8, nrow(p_values))
ggplotly(p)

filtered_table <- filtered_table[filtered_table[,8] > 0,]
write.table(filtered_table, "h1n1_filtered_pval1.tsv" ,row.names=FALSE, sep="\t", quote = FALSE)

# plotting for pval2
filtered_table <- p_values[p_values[,9] < (0.05/dim(p_values)[1]),]
discarded <- length(filtered_table[,9][filtered_table[,9] == as.double(0)])
kept <- dim(filtered_table)[1] - discarded
title <- paste("Log pval2, ", kept, " from ", dim(filtered_table)[1])
filtered_table <- filtered_table[order(filtered_table$pos2),]
p <- plot_pvalues(filtered_table, title, 9, nrow(p_values))
write.table(filtered_table, "hiv1_filtered_pval2.tsv" ,row.names=FALSE, sep="\t", quote = FALSE)

filtered_table <- filtered_table[filtered_table[,9] > 0,]
write.table(filtered_table, "h1n1_filtered_pval2.tsv" ,row.names=FALSE, sep="\t", quote = FALSE)
ggplotly(p)


# trying to take the minimum p-value between the two but that gave the same results
# as taking pval2
min_col = matrix(0,dim(p_values)[1],ncol=1)
p_values$min_pval = min_col[1,]
for (i in c(1:dim(p_values)[1])) {
  p_values[i,10] <- min(p_values[i,8], p_values[i,9])
}

# remove above p-value threshold
filtered_table <- p_values[p_values[,10] < (0.05/dim(p_values)[1]),]
# how many with 0 p-value
discarded <- length(filtered_table[,10][filtered_table[,10] == as.double(0)])
kept <- dim(filtered_table)[1] - discarded
title <- paste("Log minimum pval, ", kept, " from ", dim(filtered_table)[1])
p <- plot_pvalues(filtered_table, title, 9)
write.table(filtered_table, "filtered_pval2.tsv" ,row.names=FALSE, sep="\t", quote = FALSE)
ggplotly(p)

