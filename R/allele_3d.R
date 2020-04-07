library(cowplot)
library(ggplot2)
library(dplyr)

result.dir <- "/Users/hirak/Documents/Papers/papers/mercury-paper/dataframes/swim/"
s <- readr::read_csv(file.path(result.dir,"salmon_allele"))
merc <- readr::read_csv(file.path(result.dir,"mercury_allele"))
mm <- readr::read_csv(file.path(result.dir,"mmc_allele.2"))


combineDF <- function(a,b,c,names=c("A","B","C")) {
  data.frame(y=c(a$NumReads, b$NumReads, c$mmcollapse),
             x=c(a$count, b$count, c$count),
             method=rep(names,c(nrow(a),nrow(b),nrow(c))))
}

dat <- combineDF(s,merc,mm,names=c("Salmon","Terminus","mmcollapse"))

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
    xlab("log10(true + 1)") + ylab("log10(estimate + 1)") +
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
    ggtitle(paste0("Transcripts for which true > 0, estimate > 0 and absolute residual > ",thresh)) +
    theme(plot.title = element_text(size = 10))
  cowplot::plot_grid(g0,g1,g2,g3)
  
  
}

make3XHistogram(dat)
