# Motor Imagery ERD/S Analysis
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
cum_erds <- data.frame()[1:100,]


mua_p <- data.frame()[1:100,]
mua_p_BH <- data.frame()[1:100,]
mua_t <- data.frame()[1:100,]

# erds.plot <- ggplot() +
#   geom_rect(aes(xmin = 6.5, xmax = 9, ymin = -Inf, ymax = Inf), fill = 'black', alpha = 0.40) +
#   geom_rect(aes(xmin = -3.5, xmax = 0, ymin = -Inf, ymax = Inf), fill = 'black', alpha = 0.20) +
#   geom_vline(aes(xintercept = 0, y = NULL, size = 0.15, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
#   scale_x_continuous(breaks=c(-3.5, 0, 6.5)) +
#   ylim(-200, 200) +                                       # CHANGE THIS TO CHANGE Y AXIS
#   xlab("seconds") + ylab("%ERD/ERS") +
#   annotate("text", x = -1.75, y = -190, label = "Baseline") +
#   annotate("text", x = 3.25, y = -190, label = "Trial") +
#   annotate("text", x = 7.75, y = -190, label = "Rest") +
#   theme_cowplot(12)
# erds.plot

# Load data

folder_dir = "E:"                      # Ada Laptop
# pax_no <- readline(prompt = "How many participants?: ")
pax_no = 15 #15 with complete files
pax = list(1,2,3,4,5,7,8,9,10,11,12,13,15,16,17) # List with complete files
p_code_count = 1 #access participant in the list
dp_count = 1 # overall sessions counter for all participants x days


for (i in 1:pax_no)                                
{
  # Start here
  # p_code <- readline(prompt = "Enter participant code (1-17): ")
  p_code = pax[p_code_count]
  filename <- paste(folder_dir, "AOMI_", p_code, "_org.csv", sep = "")
  import_eeg <- read.csv(filename)  # read csv file 
  
  # Session selector
  #sesh_no <- readline(prompt = "How many sessions?: ")
  
  for (day in 1:5)
  {
    # day <- as.numeric(readline(prompt = "Which session?: "))
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
    
    sqd_ave2 <- data.frame(rowMeans(erd_pct2))
    names(sqd_ave2) <- "Signal"

    # Statistics
    
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
    
    erds.plot <- ggplot() +
      #geom_rect(aes(xmin = 6.5, xmax = 9, ymin = -Inf, ymax = Inf), fill = 'black', alpha = 0.10) +
      geom_rect(aes(xmin = -3.5, xmax = 0, ymin = -Inf, ymax = Inf), fill = 'black', alpha = 0.10) +
      geom_vline(aes(xintercept = 0, y = NULL, alpha = 0.6), size = 1, linetype = "dashed", show.legend = FALSE) +
      geom_vline(aes(xintercept = 6.5, y = NULL, alpha = 0.6), size = 1, linetype = "dashed", show.legend = FALSE) +
      geom_hline(aes(x = NULL, yintercept = 0,  alpha = 0.6), size = 1, linetype = "dashed", show.legend = FALSE) +
      xlim(-3.5, 9.8) +
      scale_x_continuous(breaks=c(-3.0, -2.0, -1.0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9.0)) +
      ylim(-200, 200) +                                       # CHANGE THIS TO CHANGE Y AXIS+
      xlab("time (seconds)") + ylab("power change (%)") +
      # annotate("text", x = -1.75, y = 190, label = "Baseline") +
      # annotate("text", x = 3.25, y = 190, label = "Trial") +
      # annotate("text", x = 7.75, y = 190, label = "Rest") +
      theme_cowplot(12)
    erds.plot
    
    ind.plot <- erds.plot +
      geom_path(data = sqd_ave_plot, aes(Time, Signal, colour=Class), size= 1.5, span = 0.1, alpha=0.9) +
      geom_ribbon(data = sqd_ave_plot, aes(Time, ymin=CI_lower, ymax=CI_upper, fill=Class), alpha=0.4) +
      labs(title = sprintf("%s, Participant %s Session %s", channel_name, p_code, day), subtitle = expression(paste("20 trials per class: Relative ", mu, " Power % at 9-11 Hz, (S-B)/B x100"))) +
      ylim(-200, 200)
    ind.plot
    
    # Cumulative
    current_session <- data.frame(list(rep(sprintf("p%s-s%s", p_code, day), 100)))
    names(current_session)[1] <- "Session"
    
    cum_erds_1_use <- cbind(current_session, sqd_ave1)
    cum_erds <- rbind(cum_erds, cum_erds_1_use)
    cum_erds_2_use <- cbind(current_session, sqd_ave2)
    cum_erds <- rbind(cum_erds, cum_erds_2_use)
    
    # # Plot for cumulative plots
    # indplot <- cbind(t_epoch_t, sqd_ave1, sqd_ave2)
    # erds.plot <- erds.plot + 
    #   geom_path(data = indplot, aes(t_epoch_t, sqd_ave2, colour=channel2), size= 1, span = 0.1, alpha=0.5) +
    #   geom_path(data = indplot, aes(t_epoch_t, sqd_ave1, colour=channel1), size= 1, span = 0.1, alpha=0.5) +
    #   scale_colour_manual(name="Legend", values=c("purple3", "forestgreen"))
    # print(erds.plot)
    
    # Mass univariate T-Tests (c/o Matt Craddock)
    runningT <- data.frame()[1:100,]
    pvals <- data.frame()[1:100,]
    crit <- data.frame()[1:100,]
    crit_BH <- data.frame()[1:100,]
    rT_counter = 1
    
    for (i in 1:100)
    {
      t_x = erd_pct1[rT_counter,]
      t_y = erd_pct2[rT_counter,]
      ttest = t.test(t_x,t_y)
      runningT <- rbind(runningT, ttest$statistic)
      pvals <- rbind(pvals, ttest$p.value)
      rT_counter = rT_counter + 1
    }
    
    #Benjamini-Hochberg Correction
    pvals <- unlist(pvals) %>% as.numeric()
    pvals_BH <- p.adjust(pvals, method="BH") %>% data.frame()
    pvals <- data.frame(pvals)
    names(pvals_BH) <- "pval"
    names(pvals) <- "pval"
    
    mua_p <- cbind(mua_p, pvals)
    mua_p_BH <- cbind(mua_p_BH, pvals_BH)
    names(mua_p)[dp_count] <- sprintf("p%s-s%s", p_code, day)
    names(mua_p_BH)[dp_count] <- sprintf("p%s-s%s_BH", p_code, day)
    
    crit <- cbind(crit, (0 + (pvals$pval <= .05)))
    names(crit) <- "crit"
    crit[crit == 0] <- NA
    crit[1:28,] <- NA
    crit[82:100,] <- NA
    
    crit_BH <- cbind(crit_BH, (0 + (pvals_BH$pval <= .05)))
    names(crit_BH) <- "crit_BH"
    crit_BH[crit_BH == 0] <- NA
    crit_BH[1:28,] <- NA
    crit_BH[82:100,] <- NA
    
    pvals <- cbind(t_epoch_t, pvals, crit)
    pvals_BH <- cbind(t_epoch_t, pvals_BH, crit_BH)
    
    sigind.plot <- ind.plot + 
      geom_line(data = pvals, aes(x = Time, y = crit-170),
                na.rm = TRUE, size = 2, colour="black") +
      annotate("text", x = 7.5, y = -170, label = "p < 0.05") +
      geom_line(data = pvals_BH, aes(x = Time, y = crit_BH-190),
                na.rm = TRUE, size = 2, colour="blue") + 
      annotate("text", x = 7.5, y = -190, label = "BH-corrected")
    sigind.plot
    
    #ggsave(sprintf("%s_%s_%s_erds.png", p_code, day, channel_name), sigind.plot, width=8.5, height=6, dpi=150, units="in", device='png')
    
    pval.plot <- ggplot(data = pvals, aes(Time, pval)) + geom_point(alpha = 0.2) + geom_point(data = pvals_BH, aes(Time, pval), colour = "blue") +
      geom_hline(aes(yintercept = 0.05, color = "red"), linetype = "dashed", show.legend = FALSE) +
      geom_vline(aes(xintercept = 0, y = NULL, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
      geom_vline(aes(xintercept = 6.5, y = NULL, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
      scale_y_continuous(breaks=c(0.05, 0.25, 0.50, 0.75, 1.00)) +
      scale_x_continuous(breaks=c(-3.5, 0.0, 6.5)) +
      xlab("Time (seconds)") + ylab("p-value") + 
      labs(title = sprintf("p-values: Left vs. Right in %s, Participant %s Session %s",  channel_name, p_code, day), subtitle = expression(paste("Benjamini-Hochberg Corrected p-values"))) + 
      theme_cowplot(12)
    pval.plot
    
    #ggsave(sprintf("%s_%s_%s_pvals.png", p_code, day, channel_name), pval.plot, width=8.5, height=6, dpi=150, units="in", device='png')
    
    mua_t <- cbind(mua_t, runningT)
    names(mua_t)[dp_count] <- sprintf("p%s-s%s", p_code, day)
    
    runningT <- cbind(t_epoch_t, runningT)
    names(runningT)[2] <- "tstat"
    runningT.plot <- ggplot(data = runningT, aes(Time, tstat)) + geom_point() + 
      geom_hline(aes(yintercept = 0.0, color = "red"), linetype = "dashed", show.legend = FALSE) +
      geom_vline(aes(xintercept = 0, y = NULL, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
      geom_vline(aes(xintercept = 6.5, y = NULL, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
      scale_x_continuous(breaks=c(-3.5, 0.0, 6.5)) +
      xlab("Time (seconds)") + ylab("t-statistic") + 
      labs(title = sprintf("t-statistic: Left vs. Right in %s, Participant %s Session %s",  channel_name, p_code, day)) + 
      theme_cowplot(12)
    runningT.plot
    
    #ggsave(sprintf("%s_%s_%s_runningt.png", p_code, day, channel_name), runningT.plot, width=8.5, height=6, dpi=150, units="in", device='png')
    
    
    dp_count = dp_count + 1
  }
  
  p_code_count = p_code_count + 1

}

# Save EEG data per participant as CSV file
write.table(mua_p, file = sprintf("AOMI_MUA_pval_%s.csv", channel_name), sep = ",", row.names = FALSE)
write.table(mua_p_BH, file = sprintf("AOMI_MUA_pval_BH_%s.csv", channel_name), sep = ",", row.names = FALSE)
write.table(mua_t, file = sprintf("AOMI_MUA_tstat_%s.csv", channel_name), sep = ",", row.names = FALSE)

write.table(cum_erds, file = sprintf("AOMI_ERDS_%s.csv", channel_name), sep = ",", row.names = FALSE)

# cum_ave1 <- data.frame()
# cum_ave1 <- rowMeans(cum_erds_1)
# cum_ave1 <- data.frame(cum_ave1)
# 
# cum_ave2 <- data.frame()
# cum_ave2 <- rowMeans(cum_erds_2)
# cum_ave2 <- data.frame(cum_ave2)
# 
# 
# sqd_ave_t <- cbind(t_e = t_epoch_t, cum_ave1, cum_ave2)
# sqd_ave_t <- data.frame(sqd_ave_t)
# 
# stdErr <- function(x) {sd(x)/ sqrt(length(x))}
# conf <- function(x) {CI(x, ci = 0.95)}
# 
# sd1 <- apply(cum_erds_1, 1, sd)
# sem1 <- apply(cum_erds_1, 1, stdErr)
# CI_1 <- apply(cum_erds_1, 1, conf)
# CI_1_upper <- CI_1[1,]
# CI_1_lower <- CI_1[3,]
# cum_ave1 <- cbind(t_epoch_t, cum_ave1, sd1, sem1, CI_1_upper, CI_1_lower)
# 
# sd2 <- apply(cum_erds_2, 1, sd)
# sem2 <- apply(cum_erds_2, 1, stdErr)
# CI_2 <- apply(cum_erds_2, 1, conf)
# CI_2_upper <- CI_2[1,]
# CI_2_lower <- CI_2[3,]
# cum_ave2 <- cbind(t_epoch_t, cum_ave2, sd2, sem2, CI_2_upper, CI_2_lower)
# 
# erds.plot + 
#   geom_ribbon(data = cum_ave2, aes(t_epoch_t, ymin=CI_2_lower, ymax=CI_2_upper, fill=channel2), alpha=0.3) +
#   geom_ribbon(data = cum_ave1, aes(t_epoch_t, ymin=CI_1_lower, ymax=CI_1_upper, fill=channel1), alpha=0.3) +
#   geom_smooth(data = sqd_ave_t, aes(t_epoch_t, cum_ave2, colour=channel2), size= 1.5, span = 0.1) +
#   geom_smooth(data = sqd_ave_t, aes(t_epoch_t, cum_ave1, colour=channel1), size= 1.5, span = 0.1) +
#   ylim(-100, 150) +
#   labs(title = sprintf("Grand Average of ERD/ERS in %s Across All Participants' %s and %s Trials",  channel_name, channel1, channel2), subtitle = expression(paste("Relative ", mu, " Power % Change at 9-11 Hz, (S-B)/B x100"))) +
#   scale_colour_manual(name="Legend", values=c("purple3", "forestgreen")) +
#   scale_fill_manual(name="Legend", values=c("purple3", "forestgreen")) +
#   theme_cowplot(12)
# 
# # Mass univariate T-Tests (c/o Matt Craddock)
# runningT <- data.frame()[1:100,]
# pvals <- data.frame()[1:100,]
# rT_counter = 1
# 
# for (i in 1:100)
# {
#   t_x = cum_erds_1[rT_counter,]
#   t_y = cum_erds_2[rT_counter,]
#   ttest = t.test(t_x,t_y)
#   runningT <- rbind(runningT, ttest$statistic)
#   pvals <- rbind(pvals, ttest$p.value)
#   rT_counter = rT_counter + 1
# }
# 
# pvals <- cbind(t_epoch_t, pvals)
# names(pvals)[2] <- "pval"
# pval.plot <- ggplot(data = pvals, aes(t_epoch_t, pval)) + geom_point()
# pval.plot

