dummy=read.csv("dummy.txt")
library(ggplot2)
## Create list based on group information

split_and_list <- function(data, group_column, value_column) {
    unique_groups <- unique(data[[group_column]])
    grouped_data <- split(data[[value_column]], data[[group_column]])
    grouped_list <- lapply(unique_groups, function(group) {
        sublist <- grouped_data[[as.character(group)]]
        return(sublist)
    })
    names(grouped_list) <- unique_groups
    return(grouped_list)
}
grouped_lists <- split_and_list(dummy, group_column = "Cluster", value_column = "adj.p.value")


##  Function for generating multiple qq-plot in one

gg_qqplot <- function(ps_list, group_labels, colors) {
  df <- data.frame()
  for (i in seq_along(ps_list)) {
    n <- length(ps_list[[i]])
    tmp <- data.frame(
      observed = -log10(sort(ps_list[[i]])),
      expected = -log10(ppoints(n)),
      group = group_labels[i]
    )
    df <- rbind(df, tmp)
  }
  
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  
  ggplot(df) +
    geom_point(aes(expected, observed, color = group, shape = group), size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    scale_color_manual(values = colors) +
    xlab(log10Pe) +
    ylab(log10Po) 
    
}


group_labels <- c("Cluster 1", "Cluster 2", "Cluster 3","Cluster 4","Cluster 5")
colors <- c("#009999", "#8B4513", "#56B4E9","#009E73","#CC79A7")

# Plotting multiple QQ plots 
gg_qqplot(grouped_lists, group_labels, colors)
