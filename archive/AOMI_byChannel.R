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
cum_erds_1 <- data.frame()[1:100,]
cum_erds_2 <- data.frame()[1:100,]

erds.plot <- ggplot() +
  geom_rect(aes(xmin = 6.5, xmax = 9, ymin = -Inf, ymax = Inf), fill = 'black', alpha = 0.40) +
  geom_rect(aes(xmin = -3.5, xmax = 0, ymin = -Inf, ymax = Inf), fill = 'black', alpha = 0.20) +
  geom_vline(aes(xintercept = 0, y = NULL, size = 0.15, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
  scale_x_continuous(breaks=c(-3.5, 0, 6.5)) +
  ylim(-200, 200) +                                       # CHANGE THIS TO CHANGE Y AXIS
  xlab("seconds") + ylab("%ERD/ERS") +
  annotate("text", x = -1.75, y = -190, label = "Baseline") +
  annotate("text", x = 3.25, y = -190, label = "Trial") +
  annotate("text", x = 7.75, y = -190, label = "Rest") +
  theme_cowplot(12)
erds.plot

# Load data

folder_dir = "C:/Users/mnpdeb/"                      # Ada Laptop
pax_no <- readline(prompt = "How many participants?: ")

for (i in 1:pax_no)                                
{
  # Start here
  p_code <- readline(prompt = "Enter participant code (1-17): ")
  filename <- paste(folder_dir, "AOMI_", p_code, "_org.csv", sep = "")
  import_eeg <- read.csv(filename)  # read csv file 
  
  # Session selector
  #sesh_no <- readline(prompt = "How many sessions?: ")
  
  for (i in 1:sesh_no)
  {
    day <- as.numeric(readline(prompt = "Which session?: "))
    eeg <- dplyr::filter(import_eeg, Session == day)
    
    Fs = 500 # sampling frequency
    
    eeg_ftd <- eegfilter(eeg[7:13], Fs = Fs, lower = 8.5, upper = 11.5, 
                         method = "butter", order = 4, forwardreverse = TRUE, scale = FALSE, plot = FALSE)
    eeg <- eeg %>% select(Trial, Class, Time)
    eeg <- cbind(eeg, eeg_ftd) %>%
      arrange(Trial, Class)
    
    Time <- eeg[1:6200,] %>% select(Time) 
    
    eeg_left <- eeg %>% 
      dplyr::filter(Class == "Left") %>%
      select(Trial, Time, C3) %>%
      spread(Trial, C3) %>%
      select(-Time)
    
    eeg_left_bp <- (eeg_left)^2
    eeg_left_bp <- cbind(Time, eeg_left_bp)

    eeg_right <- eeg %>% 
      dplyr::filter(Class == "Right") %>%
      select(Trial, Time, C3) %>%
      spread(Trial, C3) %>%
      select(-Time)
    
    eeg_right_bp <- (eeg_right)^2
    eeg_right_bp <- cbind(Time, eeg_right_bp)
    
    # Insert Trial Selection Here for eeg_left and eeg_right
    
    # Averaging epochs + stats
    
    #Initialisation
    t_epoch_duration = 2 #seconds
    t_epoch_interval = 0.1 #seconds
    t_epoch_t = seq(from = -3.500, to = 8.900, length.out = 100)
    t_epoch_t <- round(t_epoch_t, digits=1)
    t_epoch_t <- data.frame(t_epoch_t)
    names(t_epoch_t) <- "Time"
    epoch_counter = 1
    
    t_epoch_sig1 <- data.frame()
    t_epoch_sig2 <- data.frame()
    
    for (i in 1:100)
    {
      t_epoch_point <- t_epoch_t[epoch_counter,]
      t_epoch_start = t_epoch_point - t_epoch_duration*(0.5)
      t_epoch_end = t_epoch_point + t_epoch_duration*(0.5)
      eeg_bp_use1 <- dplyr::filter(eeg_left_bp, Time >= t_epoch_start & Time <= t_epoch_end)
      eeg_bp_use1 <- eeg_bp_use1[2:21]
      eeg_bp_use2 <- dplyr::filter(eeg_right_bp, Time >= t_epoch_start & Time <= t_epoch_end)
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
    
    t_epoch_sig1 <- slice(t_epoch_sig1, 1:100)
    t_epoch_sig1_b <- cbind(Time = t_epoch_t, t_epoch_sig1)
    t_epoch_sig2 <- slice(t_epoch_sig2, 1:100)
    t_epoch_sig2_b <- cbind(t_epoch_t, t_epoch_sig2)
    
    # Baseline Correction (Absolute) by B.Herrmann
    baseline1 <- dplyr::filter(t_epoch_sig1_b, between(Time, -3.5, 0)) %>% 
      select(-Time) %>%
      colMeans()
    baseline_ave1 <- data.frame(t(baseline1))
    baseline_ave1 <- baseline_ave1[rep(1,each=100),]
    
    baseline2 <- dplyr::filter(t_epoch_sig2_b, between(Time, -3.5, 0)) %>% 
      select(-Time) %>%
      colMeans()
    baseline_ave2 <- data.frame(t(baseline2))
    baseline_ave2 <- baseline_ave2[rep(1,each=100),]
    
    eeg_basecorr1 <- t_epoch_sig1 - baseline_ave1 %>% data.frame()
    eeg_basecorr2 <- t_epoch_sig2 - baseline_ave2 %>% data.frame()
    
    # %ERD/ERS
    erd_pct1 <- ((eeg_basecorr1)/baseline_ave1)*100
    erd_pct2 <- ((eeg_basecorr2)/baseline_ave2)*100

    sqd_ave1 <- data.frame(rowMeans(erd_pct1))
    names(sqd_ave1) <- "Signal"
    cum_erds_1 <- cbind(cum_erds_1, sqd_ave1)
    
    sqd_ave2 <- data.frame(rowMeans(erd_pct2))
    names(sqd_ave2) <- "Signal"
    cum_erds_2 <- cbind(cum_erds_2, sqd_ave2)
    
    stdErr <- function(x) {sd(x)/ sqrt(length(x))}
    conf <- function(x) {CI(x, ci = 0.95)}
    
    sd1 <- apply(erd_pct1, 1, sd)
    sem1 <- apply(erd_pct1, 1, stdErr)
    CI_1 <- apply(erd_pct1, 1, conf)
    CI_upper <- CI_1[1,]
    CI_lower <- CI_1[3,]
    class_col1 <- data.frame(rep.int(class1, 100)) 
    names(class_col1) <- "Class"
    sqd_ave1 <- cbind(t_epoch_t, sqd_ave1, CI_upper, CI_lower, class_col1)
    
    sd2 <- apply(erd_pct2, 1, sd)
    sem2 <- apply(erd_pct2, 1, stdErr)
    CI_2 <- apply(erd_pct2, 1, conf)
    CI_upper <- CI_2[1,]
    CI_lower <- CI_2[3,]
    class_col2 <- data.frame(rep.int(class2, 100))
    names(class_col2) <- "Class"
    sqd_ave2 <- cbind(t_epoch_t, sqd_ave2, CI_upper, CI_lower, class_col2)
    
    sqd_ave_plot <- rbind(sqd_ave1, sqd_ave2)
    
    # Individual Plots
    indplot <- erds.plot +
      geom_path(data = sqd_ave_plot, aes(Time, Signal, colour=Class), size= 2, span = 0.1, alpha=0.9) +
      geom_ribbon(data = sqd_ave_plot, aes(Time, ymin=CI_lower, ymax=CI_upper, fill=Class), alpha=0.4) +
      labs(title = sprintf("%s of Participant %s (20 Trials each class)",  channel_name, p_code), subtitle = expression(paste("Relative ", mu, " Power % Change at 9-11 Hz, (S-B)/B x100"))) +
      ylim(-150, 150)
    indplot
    
    ggsave(sprintf("%s_%s_%s.png", p_code, day, channel_name), indplot)

    # # Plot for cumulative plots
    # indplot <- cbind(t_epoch_t, sqd_ave1, sqd_ave2)
    # erds.plot <- erds.plot + 
    #   geom_path(data = indplot, aes(t_epoch_t, sqd_ave2, colour=channel2), size= 1, span = 0.1, alpha=0.5) +
    #   geom_path(data = indplot, aes(t_epoch_t, sqd_ave1, colour=channel1), size= 1, span = 0.1, alpha=0.5) +
    #   scale_colour_manual(name="Legend", values=c("purple3", "forestgreen"))
    # print(erds.plot)
    # 
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
