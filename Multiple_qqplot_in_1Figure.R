library(ggpubr)
library(ggplot2)
dummy=read.csv("dummy.txt")
#Function
gg_qqplot <- function(ps, c,title) {
    n  <- length(ps)
    df <- data.frame(
        observed = -log10(sort(ps)),
        expected = -log10(ppoints(n))
    )
    log10Pe <- expression(paste("Expected -log"[10], plain(P)))
    log10Po <- expression(paste("Observed -log"[10], plain(P)))
    ggplot(df) +
        geom_point(aes(expected, observed), shape = 1, size = 3,color=c) +
        geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
        xlab(log10Pe) + ylab(log10Po)+
         ggtitle (title) + scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 2))
}

c1=gg_qqplot(dummy$adj.p.value,"#009999","cluster1") +
theme_bw(base_size = 16) +
theme(
axis.ticks = element_line(size = 0.5),
panel.grid = element_blank()
# panel.grid = element_line(size = 0.5, color = "grey80")
)
c2=gg_qqplot(dummy$adj.p.value,"#8B4513","cluster2") +
theme_bw(base_size = 16) +
theme(
axis.ticks = element_line(size = 0.5),
panel.grid = element_blank()
# panel.grid = element_line(size = 0.5, color = "grey80")
)
c3=gg_qqplot(dummy$adj.p.value,"#56B4E9","cluster3") +
theme_bw(base_size = 16) +
theme(
axis.ticks = element_line(size = 0.5),
panel.grid = element_blank()
# panel.grid = element_line(size = 0.5, color = "grey80")
)
c4=gg_qqplot(dummy$adj.p.value,"#009E73","cluster4") +
theme_bw(base_size = 16) +
theme(
axis.ticks = element_line(size = 0.5),
panel.grid = element_blank()
# panel.grid = element_line(size = 0.5, color = "grey80")
)
c5=gg_qqplot(dummy$adj.p.value,"#CC79A7","cluster5") +
theme_bw(base_size = 16) +
theme(
axis.ticks = element_line(size = 0.5),
panel.grid = element_blank()
# panel.grid = element_line(size = 0.5, color = "grey80")
)
ggarrange(c1,c2,c3,c4,c5,
labels = c("A", "B","C","D","E"),
ncol = 2, nrow = 3)
