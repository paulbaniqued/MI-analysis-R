# Motor Imagery ERP Analysis
# Initial Packages Used
library(edf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(eegkit)
library(Rmisc)

channel_name = "C3"
class1 = "Left"
class2 = "Right"

# Set empty data frame for left and right classes
cum_erds_1 <- data.frame()[1:150,]
cum_erds_2 <- data.frame()[1:150,]

erds.plot <- ggplot() +
  geom_rect(aes(xmin = -6, xmax = -3.5, ymin = -Inf, ymax = Inf), fill = 'black', alpha = 0.40) +
  geom_rect(aes(xmin = 6.5, xmax = 9, ymin = -Inf, ymax = Inf), fill = 'black', alpha = 0.40) +
  geom_rect(aes(xmin = -3.5, xmax = 0, ymin = -Inf, ymax = Inf), fill = 'black', alpha = 0.20) +
  geom_vline(aes(xintercept = 0, y = NULL, size = 0.20, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
  geom_vline(aes(xintercept = -3.5, y = NULL, size = 0.15, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
  geom_vline(aes(xintercept = 6.5, y = NULL, size = 0.15, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
  scale_x_continuous(breaks=c(-3.5, 0, 6.5)) +
  ylim(-100, 200) +                                       # CHANGE THIS TO CHANGE Y AXIS
  xlab("seconds") + ylab("%ERD/ERS") +
  annotate("text", x = -4.75, y = -90, label = "Rest") +
  annotate("text", x = -1.75, y = -90, label = "Baseline") +
  annotate("text", x = 3.25, y = -90, label = "Trial") +
  annotate("text", x = 7.75, y = -90, label = "Rest") +
  theme_cowplot(12)
erds.plot

# Load data

folder_dir = "C:/Users/mnpdeb/"                      # Ada Laptop
pax_no <- readline(prompt = "How many participants?: ")

for (i in 1:pax_no)                                
{
  # Start here
  p_code <- readline(prompt = "Enter participant code (1-17): ")
  filename <- paste(folder_dir, "AOMI_", p_code, ".csv", sep = "")
  import_eeg <- read.csv(filename)  # read csv file 
  
  # Session selector
  sesh_no <- readline(prompt = "How many sessions?: ")
  
  for (i in 1:sesh_no)
  {
    day <- as.numeric(readline(prompt = "Which session?: "))
    eeg <- dplyr::filter(import_eeg, Session == day)
    
    Fs = 100 # sampling frequency
    
    Time <- eeg[1:1501,] %>% select(Time)
    
    eeg_left <- eeg %>% 
      dplyr::filter(Class == "left") %>%
      select(Trial, Time, C3) %>%
      mutate(eegfilter(C3, Fs = Fs, lower = 8.5, upper = 11.5, 
                                                 method = "butter", order = 4, forwardreverse = TRUE, scale = FALSE, plot = FALSE)) %>%
      spread(Trial, C3) %>%
      select(-Time) %>%
      sapply(function(x) x^2)
    eeg_left <- cbind(Time, eeg_left)

    eeg_right <- eeg %>% 
      dplyr::filter(Class == "right") %>%
      select(Trial, Time, C3) %>%
      mutate(eegfilter(C3, Fs = Fs, lower = 8.5, upper = 11.5, 
                       method = "butter", order = 4, forwardreverse = TRUE, scale = FALSE, plot = FALSE)) %>%
      spread(Trial, C3) %>%
      select(-Time) %>%
      sapply(function(x) x^2)
    eeg_right <- cbind(Time, eeg_right)
    
    # Insert Trial Selection Here for eeg_left and eeg_right
    
    # Averaging epochs + stats
    
    #Initialisation
    t_epoch_duration = 1 #seconds
    t_epoch_interval = 0.1 #seconds
    t_epoch_t = seq(from = -6.000, to = 9.000, length.out = 150)
    t_epoch_t <- round(t_epoch_t, digits=1)
    t_epoch_t <- data.frame(t_epoch_t)
    names(t_epoch_t) <- "Time"
    epoch_counter = 1
    
    t_epoch_sig1 <- data.frame()
    t_epoch_sig2 <- data.frame()
    
    for (i in 1:151)
    {
      t_epoch_point <- t_epoch_t[epoch_counter,]
      t_epoch_start = t_epoch_point - t_epoch_duration*(0.5)
      t_epoch_end = t_epoch_point + t_epoch_duration*(0.5)
      eeg_bp_use1 <- dplyr::filter(eeg_left, Time >= t_epoch_start & Time <= t_epoch_end)
      eeg_bp_use1 <- eeg_bp_use1[2:21]
      eeg_bp_use2 <- dplyr::filter(eeg_right, Time >= t_epoch_start & Time <= t_epoch_end)
      eeg_bp_use2 <- eeg_bp_use2[2:21]
      t_epoch_sig1 <- rbind(t_epoch_sig1, colMeans(eeg_bp_use1))
      t_epoch_sig2 <- rbind(t_epoch_sig2, colMeans(eeg_bp_use2))
      epoch_counter = epoch_counter + 1
    }
    
    for (epochi in 1:20)
    {
      names(t_epoch_sig1)[epochi] <- epochi
      names(t_epoch_sig2)[epochi] <- epochi
    }
    
    t_epoch_sig1 <- slice(t_epoch_sig1, 1:150)
    t_epoch_sig1_b <- cbind(Time = t_epoch_t, t_epoch_sig1)
    t_epoch_sig2 <- slice(t_epoch_sig2, 1:150)
    t_epoch_sig2_b <- cbind(t_epoch_t, t_epoch_sig2)
    
    # Baseline Correction (Absolute) by B.Herrmann
    baseline1 <- dplyr::filter(t_epoch_sig1_b, between(Time, -3.5, 0)) %>% 
      select(-Time) %>%
      colMeans()
    baseline_ave1 <- data.frame(t(baseline1))
    baseline_ave1 <- baseline_ave1[rep(1,each=150),]
    
    baseline2 <- dplyr::filter(t_epoch_sig2_b, between(Time, -3.5, 0)) %>% 
      select(-Time) %>%
      colMeans()
    baseline_ave2 <- data.frame(t(baseline2))
    baseline_ave2 <- baseline_ave2[rep(1,each=150),]
    
    eeg_basecorr1 <- t_epoch_sig1 - baseline_ave1 %>% data.frame()
    eeg_basecorr2 <- t_epoch_sig2 - baseline_ave2 %>% data.frame()
    
    # %ERD/ERS
    erd_pct1 <- ((eeg_basecorr1)/baseline_ave1)*100
    erd_pct2 <- ((eeg_basecorr2)/baseline_ave2)*100

    sqd_ave1 <- data.frame()
    sqd_ave1 <- rowMeans(erd_pct1)
    cum_erds_1 <- cbind(cum_erds_1, sqd_ave1)
    
    sqd_ave2 <- data.frame()
    sqd_ave2 <- rowMeans(erd_pct2)
    cum_erds_2 <- cbind(cum_erds_2, sqd_ave2)
    
    sqd_ave_t <- cbind(t_e = t_epoch_t, sqd_ave1, sqd_ave2)
    
    stdErr <- function(x) {sd(x)/ sqrt(length(x))}
    conf <- function(x) {CI(x, ci = 0.95)}
    
    sd1 <- apply(erd_pct1, 1, sd)
    sem1 <- apply(erd_pct1, 1, stdErr)
    CI_1 <- apply(erd_pct1, 1, conf)
    CI_1_upper <- CI_1[1,]
    CI_1_lower <- CI_1[3,]
    sqd_ave1 <- cbind(t_epoch_t, sqd_ave1, sd1, sem1, CI_1_upper, CI_1_lower)
    
    sd2 <- apply(erd_pct2, 1, sd)
    sem2 <- apply(erd_pct2, 1, stdErr)
    CI_2 <- apply(erd_pct2, 1, conf)
    CI_2_upper <- CI_2[1,]
    CI_2_lower <- CI_2[3,]
    sqd_ave2 <- cbind(t_epoch_t, sqd_ave2, sd2, sem2, CI_2_upper, CI_2_lower)
    
    
    # Individual Plots
    erds.plot +
      geom_path(data = sqd_ave_t, aes(Time, sqd_ave2, colour=class2), size= 2, span = 0.1, alpha=0.9) +
      geom_path(data = sqd_ave_t, aes(Time, sqd_ave1, colour=class1), size= 2, span = 0.1, alpha=0.9) +
      geom_ribbon(data = sqd_ave2, aes(Time, ymin=CI_2_lower, ymax=CI_2_upper, fill=class2), alpha=0.4) +
      geom_ribbon(data = sqd_ave1, aes(Time, ymin=CI_1_lower, ymax=CI_1_upper, fill=class1), alpha=0.4) +
      scale_colour_manual(name="Legend", values=c(class2, class1)) + 
      ylim(-100, 850)

    # Plot for cumulative plots
    indplot <- cbind(t_epoch_t, sqd_ave1, sqd_ave2)
    erds.plot <- erds.plot + 
      geom_path(data = indplot, aes(t_epoch_t, sqd_ave2, colour=channel2), size= 1, span = 0.1, alpha=0.5) +
      geom_path(data = indplot, aes(t_epoch_t, sqd_ave1, colour=channel1), size= 1, span = 0.1, alpha=0.5) +
      scale_colour_manual(name="Legend", values=c("purple3", "forestgreen"))
    print(erds.plot)
  
  }
  

  
  
  
  


}

cum_ave1 <- data.frame()
cum_ave1 <- rowMeans(cum_erds_1)
cum_ave1 <- data.frame(cum_ave1)

cum_ave2 <- data.frame()
cum_ave2 <- rowMeans(cum_erds_2)
cum_ave2 <- data.frame(cum_ave2)


sqd_ave_t <- cbind(t_e = t_epoch_t, cum_ave1, cum_ave2)
sqd_ave_t <- data.frame(sqd_ave_t)

stdErr <- function(x) {sd(x)/ sqrt(length(x))}
conf <- function(x) {CI(x, ci = 0.95)}

sd1 <- apply(cum_erds_1, 1, sd)
sem1 <- apply(cum_erds_1, 1, stdErr)
CI_1 <- apply(cum_erds_1, 1, conf)
CI_1_upper <- CI_1[1,]
CI_1_lower <- CI_1[3,]
cum_ave1 <- cbind(t_epoch_t, cum_ave1, sd1, sem1, CI_1_upper, CI_1_lower)

sd2 <- apply(cum_erds_2, 1, sd)
sem2 <- apply(cum_erds_2, 1, stdErr)
CI_2 <- apply(cum_erds_2, 1, conf)
CI_2_upper <- CI_2[1,]
CI_2_lower <- CI_2[3,]
cum_ave2 <- cbind(t_epoch_t, cum_ave2, sd2, sem2, CI_2_upper, CI_2_lower)

erds.plot + 
  geom_ribbon(data = cum_ave2, aes(t_epoch_t, ymin=CI_2_lower, ymax=CI_2_upper, fill=channel2), alpha=0.3) +
  geom_ribbon(data = cum_ave1, aes(t_epoch_t, ymin=CI_1_lower, ymax=CI_1_upper, fill=channel1), alpha=0.3) +
  geom_smooth(data = sqd_ave_t, aes(t_epoch_t, cum_ave2, colour=channel2), size= 1.5, span = 0.1) +
  geom_smooth(data = sqd_ave_t, aes(t_epoch_t, cum_ave1, colour=channel1), size= 1.5, span = 0.1) +
  ylim(-100, 150) +
  labs(title = sprintf("Grand Average of ERD/ERS in %s Across All Participants' %s and %s Trials",  channel_name, channel1, channel2), subtitle = expression(paste("Relative ", mu, " Power % Change at 9-11 Hz, (S-B)/B x100"))) +
  scale_colour_manual(name="Legend", values=c("purple3", "forestgreen")) +
  scale_fill_manual(name="Legend", values=c("purple3", "forestgreen")) +
  theme_cowplot(12)

# Mass univariate T-Tests (c/o Matt Craddock)
runningT <- data.frame()[1:100,]
pvals <- data.frame()[1:100,]
rT_counter = 1

for (i in 1:100)
{
  t_x = cum_erds_1[rT_counter,]
  t_y = cum_erds_2[rT_counter,]
  ttest = t.test(t_x,t_y)
  runningT <- rbind(runningT, ttest$statistic)
  pvals <- rbind(pvals, ttest$p.value)
  rT_counter = rT_counter + 1
}

pvals <- cbind(t_epoch_t, pvals)
names(pvals)[2] <- "pval"
pval.plot <- ggplot(data = pvals, aes(t_epoch_t, pval)) + geom_point()
pval.plot
