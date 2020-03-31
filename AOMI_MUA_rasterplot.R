library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(Rmisc)
library(reshape2)

channel = "C4"

pval_df <- read.csv(sprintf("AOMI_MUA_pval_%s.csv", channel))  # read csv file 
p_datapoints = as.numeric(ncol(pval_df))
pval_df <- melt(pval_df)

p_t = seq(from = -3.500, to = 8.900, length.out = 100) %>%
  round(digits=1) %>% 
  rep(times = p_datapoints) %>%
  data.frame()
names(p_t) <- "Time"

p_r <- cbind(p_t, pval_df)
names(p_r)[2] <- "Sessions"
names(p_r)[3] <- "p_value"

rasterplot <- ggplot(p_r, aes(x=Time, y=Sessions, fill=p_value)) +
  geom_raster() +
  scale_fill_gradientn(colours = cols) +
  coord_fixed(ratio = 1/8) +
  xlim(0,6.5) +
  scale_fill_gradientn(
    colours=c("black", "red", "orange", "yellow", "white"),
    values = c(0.00, 0.05, 0.25, 0.75, 1.00))
rasterplot + labs(title = sprintf("%s, Mass Univariate Analysis (p-value), N=17", channel)) + 
  theme_cowplot(12)

tstat_df <- read.csv(sprintf("AOMI_MUA_tstat_%s.csv", channel))  # read csv file 
t_datapoints = as.numeric(ncol(tstat_df))
tstat_df <- melt(tstat_df)

t_t = seq(from = -3.500, to = 8.900, length.out = 100) %>%
  round(digits=1) %>% 
  rep(times = t_datapoints) %>%
  data.frame()
names(t_t) <- "Time"

t_r <- cbind(t_t, tstat_df)
names(t_r)[2] <- "Sessions"
names(t_r)[3] <- "tstat"

cols <- rainbow(7)[-7]

rasterplot <- ggplot(t_r, aes(x=Time, y=Sessions, fill=tstat)) +
  geom_raster() +
  scale_fill_gradientn(colours = cols) +
  coord_fixed(ratio = 1/8) +
  xlim(0,6.5)
  # scale_fill_gradientn(
  #   colours=c("black", "red", "orange", "yellow", "white"),
  #   values = c(0.00, 0.05, 0.25, 0.75, 1.00))
rasterplot + labs(title = sprintf("%s, Mass Univariate Analysis (t-statistic), N=17", channel)) + 
  theme_cowplot(12)

