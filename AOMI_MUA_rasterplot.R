library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(Rmisc)
library(reshape2)
library(ggbio)

pval_df <- read.csv("AOMI_MUA_C4.csv")  # read csv file 
datapoints = as.numeric(ncol(pval_df))
pval_df <- melt(pval_df)

t = seq(from = -3.500, to = 8.900, length.out = 100) %>%
  round(digits=1) %>% 
  rep(times = datapoints) %>%
  data.frame()
names(t) <- "Time"

r <- cbind(t, pval_df)
names(r)[2] <- "Sessions"
names(r)[3] <- "p_value"

cols <- rainbow(7)[-7]

rasterplot <- ggplot(r, aes(x=Time, y=Sessions, fill=p_value)) +
  geom_raster() +
  scale_fill_gradientn(colours = cols) +
  coord_fixed(ratio = 1/8) +
  xlim(0,6.5) +
  scale_fill_gradientn(
    colours=c("black", "red", "orange", "yellow", "white"),
    values = c(0.00, 0.05, 0.25, 0.75, 1.00))
rasterplot + labs(title = "C4, Mass Univariate Analysis, N=17") + 
  theme_cowplot(12)
