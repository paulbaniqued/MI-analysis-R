# Motor Imagery ERP Analysis
# Initial Packages Used
library(edf)
library(ggplot2)
library(dplyr)
library(cowplot)
library(eegUtils)
library(eegkit)
library(Rmisc)

# import EDF reader package and read .EDF file in the directory
import_eeg <- read.edf("C:/Users/mnpdeb/EDF/1_s1.edf", read.annotations = TRUE, header.only = FALSE)
#import_eeg <- read.edf("C:/Users/paclab/EEG_MI/20190723140755_P04_Stream.edf", read.annotations = TRUE, header.only = FALSE)

# construct main data frame for EEG anaylsis
raw_eeg <- data.frame()

t = data.frame(import_eeg$signal$C3$t)
raw_eeg <- cbind(t)
names(raw_eeg)[1] <- "t"

C3 = data.frame(import_eeg$signal$C3$data)
raw_eeg <- cbind(raw_eeg, C3)
names(raw_eeg)[2] <- "C3"

C4 = data.frame(import_eeg$signal$C4$data)
raw_eeg <- cbind(raw_eeg, C4)
names(raw_eeg)[3] <- "C4"

CP5 = data.frame(import_eeg$signal$CP5$data)
raw_eeg <- cbind(raw_eeg, CP5)
names(raw_eeg)[4] <- "CP5"

CP6 = data.frame(import_eeg$signal$CP6$data)
raw_eeg <- cbind(raw_eeg, CP6)
names(raw_eeg)[5] <- "CP6"

C1 = data.frame(import_eeg$signal$C1$data)
raw_eeg <- cbind(raw_eeg, C1)
names(raw_eeg)[6] <- "C1"

C2 = data.frame(import_eeg$signal$C2$data)
raw_eeg <- cbind(raw_eeg, C2)
names(raw_eeg)[7] <- "C2"

FC5 = data.frame(import_eeg$signal$FC5$data)
raw_eeg <- cbind(raw_eeg, FC5)
names(raw_eeg)[8] <- "FC5"

FC6 = data.frame(import_eeg$signal$FC6$data)
raw_eeg <- cbind(raw_eeg, FC6)
names(raw_eeg)[9] <- "FC6"

raw_eeg

# Markers in a separate data frame
markers <- data.frame()

tm = data.frame(import_eeg$events$onset)
markers <- cbind(tm)
names(markers)[1] <- "tm"

trigger = data.frame(import_eeg$events$annotation)
markers <- cbind(markers, trigger)
names(markers)[2] <- "trigger"
markers$trigger <- as.numeric(markers$trigger)
within(markers, levels(trigger)[levels(trigger) == "Trigger#1"] <- "1")
markers$trigger = markers$trigger - 1
markers

markers_reduced <- dplyr::filter(markers, trigger > 2 & trigger < 5, tm >30) # tm>50 to select markers at session start 
markers_reduced
markers_left <- dplyr::filter(markers_reduced, trigger == 3)
markers_left
markers_right <- dplyr::filter(markers_reduced, trigger == 4)
markers_right
markers_reduced$trigger <- as.factor(markers_reduced$trigger)

# Visualise Markers, Raw plot of EEG for 1 channel

ggplot(raw_eeg, aes(t, C3)) +
  geom_vline(data = markers_reduced, aes(xintercept = tm, y = NULL, color = trigger), linetype = "dashed") +
  geom_path() +
  scale_color_manual(values=c('#E69F00', '#56B4E9')) +
  theme_cowplot()

# Save EEG as CSV file
#save_eeg <- cbind(raw_eeg, markers$trigger)
#write.table(raw_eeg, file = "P05_Session8.csv", sep = ",")

# Spatial Filter

# Surface Laplacian
# Right Classes: C3, FC5, C1, CP5
eeg_oc1 <- raw_eeg %>% 
  select(C3, FC5, C1, CP5) %>% # Select left hemisphere channels; right hand trials
  #mutate(C3 = C3*4 - FC5*1 - C1*1 - CP5*1) %>% # Coefficients: 4;-1;-1;-1
  mutate(C3 = C3*1 - FC5*0 - C1*0 - CP5*0) %>% # SL DISABLED!!!!
  select(C3) # Reduce dimensionality to 1 output channel
eeg_oc1

# Left Classes: C4, FC6, C2, CP6
eeg_oc2 <- raw_eeg %>% 
  select(C4, FC6, C2, CP6) %>% # Select right hemisphere channels; left hand trials
  #mutate(C4 = C4*4 - FC6*1 - C2*1 - CP6*1) %>% # Coefficients: 4;-1;-1;-1
  mutate(C4 = C4*1 - FC6*0 - C2*0 - CP6*0) %>% # SL DISABLED!!!
  select(C4) # Reduce dimensionality to 1 output channel
eeg_oc2

eeg_spatialftd <- cbind(eeg_oc1, eeg_oc2)
eeg_spatialftd

# Temporal Filter

# Bandpass filter for Mu (8-12 Hz)
eeg_tempftd <- eegfilter(eeg_spatialftd, Fs = 500, lower = 8.5, upper = 11.5, method = "butter",
          order = 4, forwardreverse = TRUE, 
          scale = FALSE, plot = FALSE)
eeg_tempftd <- data.frame(eeg_tempftd)
eeg_tempftd <- cbind(t = raw_eeg$t, eeg_tempftd)

# Plot of band-pass filtered output channels

ggplot() +
  geom_vline(data = markers_reduced, aes(xintercept = tm, y = NULL, color = trigger, alpha = 0.6, ), linetype = "dashed", show.legend = FALSE) +
  scale_color_manual(values=c("#00AFBB", "#E7B800"), name = "Trigger", labels = c("Left", "Right")) +
  geom_path(data = eeg_tempftd, aes(t,C4), colour="#E7B800", size=1.0) + # yellow, left
  geom_path(data = eeg_tempftd, aes(t,C3), colour="#00AFBB", size=1.0) + # blue, right
  ylim(-20, 20) +
  #xlim(0, 500) + #whole data
  xlim(206.5, 226.5) + #specific time
  theme_cowplot()
  
# Before proceeding with Epoching Trials, we must first define which channel (C3laplacian, C4laplacian)
# and which trials are we looking at:
# LEFT TRIALS: (C4 Activity High, C3 Baseline)
# RIGHT TRIALS: (C3 Activity Low, C4 Baseline)

# Settings: Change this
trial_name = "Right"
trial_markers = markers_right
channel1 = "C3"
channel2 = "C4"

# Time-based Epoching

#Initialisation
sampling_frequency = 500 #Hertz
epoch_length = 15 #seconds
samples = (sampling_frequency * epoch_length)
t_e <- seq(from = -6.000, to = epoch_length, length.out = 7501)
t_e <- round(t_e, digits = 3)
t_e <- data.frame(t_e)

# All the trials -- change trial in settings
epoch_counter = 1
all_trial_signals1 <- data.frame()[1:7501, ]

for (i in 1:20) #Primary Channel of Interest
{
  epoch_x <- data.frame()
  epoch_start = trial_markers$tm[epoch_counter]-6 # -- change trial in settings
  epoch_end = epoch_start + epoch_length
  epoch_start_i = which(round(eeg_tempftd$t, 3) == round(epoch_start, 3))
  epoch_end_i = epoch_start_i + samples
  epoch_x <- eeg_tempftd[epoch_start_i:epoch_end_i,]
  epoch_x <- epoch_x$C3                             # -- CHANGE MANUALLY
  epoch_x <- data.frame(epoch_x)
  all_trial_signals1 <- cbind(all_trial_signals1, epoch_x)
  names(all_trial_signals1)[epoch_counter] <- sprintf("%s_Epoch_%s", channel1, epoch_counter)
  epoch_counter = epoch_counter + 1
}

epoch_counter = 1
all_trial_signals2 <- data.frame()[1:7501, ]

for (i in 1:20) #Secondary Channel (non-interest)
{
  epoch_x <- data.frame()
  epoch_start = trial_markers$tm[epoch_counter]-6 # -- change trial in settings
  epoch_end = epoch_start + epoch_length
  epoch_start_i = which(round(eeg_tempftd$t, 3) == round(epoch_start, 3)) # 3 decimal places
  epoch_end_i = epoch_start_i + samples
  epoch_x <- eeg_tempftd[epoch_start_i:epoch_end_i,]
  epoch_x <- epoch_x$C4                            # -- CHANGE MANUALLY
  epoch_x <- data.frame(epoch_x)
  all_trial_signals2 <- cbind(all_trial_signals2, epoch_x)
  names(all_trial_signals2)[epoch_counter] <- sprintf("%s_Epoch_%s", channel2, epoch_counter)
  epoch_counter = epoch_counter + 1
}

all_trials <- data.frame()
all_trials <- cbind(t_e, all_trial_signals1, all_trial_signals2)
all_trials

# Signal Averaging
all_ave1 <- data.frame()
all_ave1 <- rowMeans(all_trial_signals1)
all_ave1 <- data.frame(all_ave1)
all_ave1 <- cbind(t_e = all_trials$t_e, Ave = all_ave1)

all_ave2 <- data.frame()
all_ave2 <- rowMeans(all_trial_signals2)
all_ave2 <- data.frame(all_ave2)
all_ave2 <- cbind(t_e = all_trials$t_e, Ave = all_ave2)

# Visualise Bandpass Filtered and Epoched Signals against Averaged Signals
ggplot() +
  geom_vline(aes(xintercept = 0, y = NULL, size = 0.25, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
  geom_vline(aes(xintercept = -3.5, y = NULL, size = 0.25, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
  geom_vline(aes(xintercept = 6.5, y = NULL, size = 0.25, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
  geom_rect(aes(xmin = -6, xmax = -3.5, ymin = -Inf, ymax = Inf), fill = 'red', alpha = 0.2) +
  geom_rect(aes(xmin = 6.5, xmax = 9, ymin = -Inf, ymax = Inf), fill = 'red', alpha = 0.2) +
  geom_rect(aes(xmin = -3.5, xmax = 0, ymin = -Inf, ymax = Inf), fill = 'orange', alpha = 0.20) +
  geom_rect(aes(xmin = 0, xmax = 6.5, ymin = -Inf, ymax = Inf), fill = 'green', alpha = 0.20) +
  geom_path(data = all_trials, aes(t_e,C3_Epoch_1), colour="#E69F00", size=0.5, alpha=0.4) +
  geom_path(data = all_trials, aes(t_e,C3_Epoch_2), colour="#56B4E9", size=0.5, alpha=0.4) +
  geom_path(data = all_trials, aes(t_e,C3_Epoch_3), colour="#009E73", size=0.5, alpha=0.4) +
  geom_path(data = all_trials, aes(t_e,C3_Epoch_4), colour="#F0E442", size=0.5, alpha=0.4) +
  geom_path(data = all_trials, aes(t_e,C3_Epoch_5), colour="#0072B2", size=0.5, alpha=0.4) +
  geom_path(data = all_trials, aes(t_e,C3_Epoch_6), colour="#D55E00", size=0.5, alpha=0.4) +
  geom_path(data = all_trials, aes(t_e,C3_Epoch_7), colour="#CC79A7", size=0.5, alpha=0.4) +
  geom_path(data = all_ave1, aes(t_e, all_ave1), colour="black", size=0.8, alpha=1) +
  #geom_path(data = all_ave2, aes(t_e, all_ave2), colour="yellow", size=0.8, alpha=1) +
  xlim(-6, 9) +
  ylim(-10, 10) +
  xlab("seconds") + ylab("uV") + 
  labs(title = sprintf("%s: All %s Trials and Trial Average", channel1, trial_name), subtitle = "Bandpass-filtered (Mu 9-11 Hz), Epoch Length = 15s") +
  theme_cowplot(12)
  
# Spectral Bandpower Estimation (Squaring)

eeg_bp1 <- data.frame()
eeg_bp1 <- (all_trial_signals1)^2
eeg_bp1 <- cbind(t_e = all_trials$t_e, eeg_bp1)
eeg_bp2 <- data.frame()
eeg_bp2 <- (all_trial_signals2)^2
eeg_bp2 <- cbind(t_e = all_trials$t_e, eeg_bp2)

# sqd_ave1 <- data.frame()
# sqd_ave1 <- rowMeans(eeg_bp1)
# sqd_ave1 <- data.frame(sqd_ave1)
# sqd_ave2 <- data.frame()
# sqd_ave2 <- rowMeans(eeg_bp2)
# sqd_ave2 <- data.frame(sqd_ave2)
# sqd_ave_t <- cbind(t_e = all_trials$t_e, sqd_ave1, sqd_ave2)
# sqd_ave_t <- data.frame(sqd_ave_t)

# # Visualise Squared Signals (Spectral Bandpower)
# ggplot() +
#   geom_vline(aes(xintercept = 0, y = NULL, size = 0.25, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
#   geom_vline(aes(xintercept = -3.5, y = NULL, size = 0.25, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
#   geom_vline(aes(xintercept = 6.5, y = NULL, size = 0.25, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
#   geom_rect(aes(xmin = -6, xmax = -3.5, ymin = -Inf, ymax = Inf), fill = 'black', alpha = 0.5) +
#   geom_rect(aes(xmin = 6.5, xmax = 9, ymin = -Inf, ymax = Inf), fill = 'black', alpha = 0.5) +
#   geom_rect(aes(xmin = -3.5, xmax = 0, ymin = -Inf, ymax = Inf), fill = 'orange', alpha = 0.10) +
#   geom_rect(aes(xmin = 0, xmax = 6.5, ymin = -Inf, ymax = Inf), fill = 'green', alpha = 0.10) +
#   geom_path(data = eeg_bp_t, aes(t_e,C3_Epoch_1), colour="#E69F00", size=0.5, alpha=0.1) +
#   geom_path(data = eeg_bp_t, aes(t_e,C3_Epoch_2), colour="#56B4E9", size=0.5, alpha=0.1) +
#   geom_path(data = eeg_bp_t, aes(t_e,C3_Epoch_3), colour="#009E73", size=0.5, alpha=0.1) +
#   geom_path(data = eeg_bp_t, aes(t_e,C3_Epoch_4), colour="#F0E442", size=0.5, alpha=0.1) +
#   geom_path(data = eeg_bp_t, aes(t_e,C3_Epoch_5), colour="#0072B2", size=0.5, alpha=0.1) +
#   geom_path(data = eeg_bp_t, aes(t_e,C3_Epoch_6), colour="#D55E00", size=0.5, alpha=0.1) +
#   geom_path(data = eeg_bp_t, aes(t_e,C3_Epoch_7), colour="#CC79A7", size=0.5, alpha=0.1) +
#   geom_path(data = sqd_ave_t, aes(t_e,sqd_ave1), colour="black", size=0.8, alpha=1) +
#   geom_path(data = sqd_ave_t, aes(t_e,sqd_ave2), colour="steelblue", size=0.8, alpha=0.6) +
#   xlim(-6, 9) +
#   ylim(0, 50) +
#   xlab("seconds") + ylab("uV^2") + 
#   labs(title = sprintf("Mu Band Power in %s Across All %s Trials", channel1, trial_name), subtitle = expression(paste("Spatially-filtered (Surface Laplacian), Bandpass-filtered (", mu, " 8-12 Hz), Epoch Length = 10s"))) +
#   theme_cowplot(12)

# Averaging Over Time

#Initialisation
t_epoch_duration = 2 #seconds
t_epoch_interval = 0.1 #seconds
sampling_frequency = 500 #samples / sec
t_epoch <- data.frame()[1:100, ]
t_epoch_t = seq(from = -6.000, to = 9.000, length.out = 100)
t_epoch_t <- round(t_epoch_t, digits = 1)
t_epoch_t <- data.frame(t_epoch_t)
t_epoch <- cbind(t_epoch_t)
sqd_ave_t$t_e <- round(all_trials$t_e, digits = 1)
epoch_counter = 1

t_epoch_sig1 <- data.frame()
t_epoch_sig2 <- data.frame()

for (i in 1:101)
{
  t_epoch_point <- t_epoch_t[epoch_counter,]
  t_epoch_start = t_epoch_point - t_epoch_duration*(0.5)
  t_epoch_end = t_epoch_point + t_epoch_duration*(0.5)
  eeg_bp_use1 <- dplyr::filter(eeg_bp1, t_e >= t_epoch_start & t_e <= t_epoch_end)
  eeg_bp_use1 <- eeg_bp_use1[2:21]
  eeg_bp_use2 <- dplyr::filter(eeg_bp2, t_e >= t_epoch_start & t_e <= t_epoch_end)
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
t_epoch_sig1_b <- cbind(t_epoch_t, t_epoch_sig1)
t_epoch_sig2 <- slice(t_epoch_sig2, 1:100)
t_epoch_sig2_b <- cbind(t_epoch_t, t_epoch_sig2)

# Baseline Correction (Absolute) by B.Herrmann
baseline1 <- dplyr::filter(t_epoch_sig1_b, between(t_epoch_t, -3.5, 0))
baseline1 <- baseline1[2:21]
baseline_ave1 <- colMeans(baseline1)
baseline_ave1 <- data.frame(t(baseline_ave1))
baseline_ave1 <- baseline_ave1[rep(1,each=100),]

baseline2 <- dplyr::filter(t_epoch_sig2_b, between(t_epoch_t, -3.5, 0))
baseline2 <- baseline2[2:21]
baseline_ave2 <- colMeans(baseline2)
baseline_ave2 <- data.frame(t(baseline_ave2))
baseline_ave2 <- baseline_ave2[rep(1,each=100),]


eeg_basecorr1 <- t_epoch_sig1 - baseline_ave1
eeg_basecorr1 <- data.frame(eeg_basecorr1)
eeg_basecorr2 <- t_epoch_sig2 - baseline_ave2
eeg_basecorr2 <- data.frame(eeg_basecorr2)

# %ERD/ERS
erd_pct1 <- ((eeg_basecorr1)/baseline_ave1)*100
erd_pct2 <- ((eeg_basecorr2)/baseline_ave2)*100
erd_pct <- data.frame()
erd_pct <- cbind(t_epoch_t, erd_pct1, erd_pct2)

sqd_ave1 <- data.frame()
sqd_ave1 <- rowMeans(erd_pct1)
sqd_ave1 <- data.frame(sqd_ave1)
sqd_ave2 <- data.frame()
sqd_ave2 <- rowMeans(erd_pct2)
sqd_ave2 <- data.frame(sqd_ave2)
sqd_ave_t <- cbind(t_e = t_epoch_t, sqd_ave1, sqd_ave2)
sqd_ave_t <- data.frame(sqd_ave_t)

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

ggplot() +
  geom_rect(aes(xmin = -6, xmax = -3.5, ymin = -Inf, ymax = Inf), fill = 'black', alpha = 0.40) +
  geom_rect(aes(xmin = 6.5, xmax = 9, ymin = -Inf, ymax = Inf), fill = 'black', alpha = 0.40) +
  geom_rect(aes(xmin = -3.5, xmax = 0, ymin = -Inf, ymax = Inf), fill = 'black', alpha = 0.20) +
  geom_vline(aes(xintercept = 0, y = NULL, size = 0.20, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
  geom_vline(aes(xintercept = -3.5, y = NULL, size = 0.15, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
  geom_vline(aes(xintercept = 6.5, y = NULL, size = 0.15, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
  geom_ribbon(data = sqd_ave2, aes(t_epoch_t, ymin=CI_2_lower, ymax=CI_2_upper), fill="#e7b800", alpha=0.4) +
  geom_ribbon(data = sqd_ave1, aes(t_epoch_t, ymin=CI_1_lower, ymax=CI_1_upper), fill="#00afbb", alpha=0.4) +
  geom_smooth(data = sqd_ave_t, aes(t_epoch_t, sqd_ave2, colour=channel2), size= 1.5, alpha=0.0, span = 0.1) +
  geom_smooth(data = sqd_ave_t, aes(t_epoch_t, sqd_ave1, colour=channel1), size= 1.5, alpha=0.0, span = 0.1) +
  scale_x_continuous(breaks=c(-3.5, 0, 6.5)) +
  ylim(-200, 200) +                                       # CHANGE THIS TO CHANGE Y AXIS
  xlab("seconds") + ylab("%ERD/ERS") + 
  annotate("text", x = -4.75, y = -190, label = "Rest") +
  annotate("text", x = -1.75, y = -190, label = "Baseline") +
  annotate("text", x = 3.25, y = -190, label = "Trial") +
  annotate("text", x = 7.75, y = -190, label = "Rest") +
  labs(title = sprintf("Percent Change in ERD/ERS of %s and %s Across All %s Trials", channel1, channel2, trial_name), subtitle = expression(paste("Relative ", mu, " Power at 9-11 Hz, (S-B)/B x100"))) +
  scale_colour_manual(name="Legend", values=c("#00afbb", "#e7b800")) +
  theme_cowplot(12)


