library(cowplot)
library(ggplot2)
library(dplyr)

result.dir <- "/Users/hirak/Documents/Papers/papers/mercury-paper/dataframes/swim/"
s <- readr::read_csv(file.path(result.dir,"salmon_truth_1"))
merc <- readr::read_csv(file.path(result.dir,"merc_truth_1_1"))


combineDF <- function(a,b,names=c("A","B")) {
  data.frame(y=c(a$salmon_1, b$mercury_1),
             x=c(a$sample_01, b$sample_01),
             method=rep(names,c(nrow(a),nrow(b))))
}

dat <- combineDF(s,merc,names=c("Salmon","Terminus"))

make3XHistogram <- function(dat, thresh=.5, 
                            alsoShowGood=TRUE, alpha=.25) {
  dat2 <- dat %>% filter(x > 0 | y > 0) %>%
    mutate(x = log10(x + 1), y = log10(y + 1))
  if (!alsoShowGood) {
    dat2 <- dat2 %>% filter(abs(x - y) > thresh)
  }
  dat2 <- dat2 %>% 
    mutate(bad=ifelse(abs(x - y) > thresh,"red","black"))
  g0 <- ggplot(dat2, aes(x,y,color=bad)) + 
    geom_point(size=.5, alpha=alpha, show.legend=FALSE) + 
    geom_abline(intercept=0, slope=1, col="blue") +
    geom_abline(intercept=.5, slope=1, col="red") +
    geom_abline(intercept=-.5, slope=1, col="red") +
    facet_wrap(~method) + 
    xlab("true") + ylab("Estimate") +
    theme_bw() +
    scale_colour_manual(values=c("black","red"))
  g1 <- dat %>% filter(x == 0 & y >= 0.1) %>%
    mutate(x = log10(x + 1), y = log10(y + 1)) %>%
    ggplot(aes(y, fill=method)) + 
    geom_histogram(position="dodge") + 
    theme_bw() +
    ggtitle("Transcripts for which true = 0 and estimate >= 0.1") + xlab("log10(estimate + 1)") +
    theme(plot.title = element_text(size = 10))
  g2 <- dat %>% filter(x > 1 & y < 1) %>%
    mutate(x = log10(x + 1), y = log10(y + 1)) %>%
    ggplot(aes(x, fill=method)) + 
    geom_histogram(position="dodge") + 
    theme_bw() +
    ggtitle("Transcripts for which estimate < 1") + xlab("log10(true + 1)") + 
    theme(plot.title = element_text(size = 10))
  g3 <- dat %>% filter(x > 0 & y > 0) %>%
    mutate(x = log10(x + 1), y = log10(y + 1)) %>%
    filter(abs(x - y) > thresh) %>%
    mutate(residual = y - x) %>%
    ggplot(aes(residual, fill=method)) + 
    geom_histogram(position="dodge") +
    theme_bw() +
    ggtitle(paste0("Transcripts for which true > 0, estimate > 0 and residual > ",thresh)) +
    theme(plot.title = element_text(size = 10))
  cowplot::plot_grid(g0,g1,g2,g3)
}

png(filename="/Users/hirak/Documents/Papers/papers/mercury-paper/swim_salmon_vs_terminus.png")
make3XHistogram(dat)
