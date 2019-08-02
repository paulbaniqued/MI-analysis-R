# Motor Imagery ERP Analysis
# Initial Packages Used
library(edf)
library(ggplot2)
library(dplyr)
library(cowplot)
library(eegUtils)
library(eegkit)

# import EDF reader package and read .EDF file in directory
import_eeg <- read.edf("C:/Users/Paul/EEG_MI/20190726151358_P05_Stream.edf", read.annotations = TRUE, header.only = FALSE)

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

markers_reduced <- dplyr::filter(markers, trigger > 0, tm >50) # tm>50 to select markers at session start 
markers_reduced
markers_left <- dplyr::filter(markers_reduced, trigger == 1)
markers_left
markers_right <- dplyr::filter(markers_reduced, trigger == 2)
markers_right
markers_reduced$trigger <- as.factor(markers_reduced$trigger)

# Raw plot of EEG for 1 channel

ggplot(raw_eeg, aes(t, C3)) +
  geom_vline(data = markers_reduced, aes(xintercept = tm, y = NULL, color = trigger), linetype = "dashed") +
  geom_path() +
  scale_color_manual(values=c('#E69F00', '#56B4E9')) +
  theme_cowplot()

# Save EEG as CSV file
#save_eeg <- cbind(raw_eeg, markers$trigger)
#write.table(raw_eeg, file = "P05_Session8.csv", sep = ",")

# Spatial Filter

# Channel Selection 
# Right Classes: C3, FC5, C1, CP5
eeg_oc1 <- raw_eeg %>% 
  select(C3, FC5, C1, CP5) %>% # Select left hemisphere channels; right hand trials
  mutate(C3 = C3*4 - FC5*1 - C1*1 - CP5*1) %>% # Coefficients: 4;-1;-1;-1
  select(C3) # Reduce dimensionality to 1 output channel
eeg_oc1

# Left Classes: C4, FC6, C2, CP6
eeg_oc2 <- raw_eeg %>% 
  select(C4, FC6, C2, CP6) %>% # Select right hemisphere channels; left hand trials
  mutate(C4 = C4*4 - FC6*1 - C2*1 - CP6*1) %>% # Coefficients: 4;-1;-1;-1
  select(C4) # Reduce dimensionality to 1 output channel
eeg_oc2

eeg_spatialftd <- cbind(eeg_oc1, eeg_oc2)
eeg_spatialftd

# Temporal Filter

# Bandpass filter for Mu (8-12 Hz)
eeg_tempftd <- eegfilter(eeg_spatialftd, Fs = 500, lower = 7.5, upper = 12.5, method = "butter",
          order = 4L, forwardreverse = TRUE, 
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
  xlim(0, 500) + #whole data
  #xlim(206.5, 226.5) + #specific time
  theme_cowplot()
  
# Before proceeding with Epoching Trials, we must first define which channel (C3laplacian, C4laplacian)
# and which trials are we looking at:
# LEFT TRIALS: (C4 Activity High, C3 Baseline)
# RIGHT TRIALS: (C3 Activity Low, C4 Baseline)

# Settings: Change this
trial_name = "Left"
trial_markers = markers_left
channel = "C4"

# Time-based Epoching

#Initialisation
sampling_frequency = 500 #Hertz
epoch_length = 10 #seconds
samples = (sampling_frequency * epoch_length)
t_e <- seq(from = 0.000, to = epoch_length, length.out = 5001)
t_e <- round(t_e, digits = 3)
t_e <- data.frame(t_e)

# All the trials (one side) -- change trial in settings
epoch_counter = 1
all_trials <- data.frame()
all_trial_signals <- data.frame()[1:5001, ]
channel_i = toString(channel)
for (i in 1:20)
{
  epoch_x <- data.frame()
  epoch_start = trial_markers$tm[epoch_counter] # -- change trial in settings
  epoch_end = epoch_start + epoch_length
  epoch_start_i = which(eeg_tempftd$t == epoch_start) 
  epoch_end_i = epoch_start_i + samples
  epoch_x <- eeg_tempftd[epoch_start_i:epoch_end_i,]
  epoch_x <- epoch_x$C4 # -- change channel in settings                      CHANGE MANUALLY FOR NOW
  epoch_x <- data.frame(epoch_x)
  all_trial_signals <- cbind(all_trial_signals, epoch_x)
  names(all_trial_signals)[epoch_counter] <- sprintf("Epoch_%s", epoch_counter)
  epoch_counter = epoch_counter + 1
}
all_trials <- cbind(t_e, all_trial_signals)
all_trials

# Signal Averaging
all_ave <- data.frame()
all_ave <- rowMeans(all_trial_signals)
all_ave <- data.frame(all_ave)
all_ave <- cbind(t_e = all_trials$t_e, Ave = all_ave)

# Visualise Spatial and Bandpass Filtered and Epoched Signals against Averaged Signals
ggplot() +
  geom_vline(aes(xintercept = 0, y = NULL, size = 0.5, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
  geom_path(data = all_trials, aes(t_e,Epoch_1), colour="#E69F00", size=0.5, alpha=0.4) +
  geom_path(data = all_trials, aes(t_e,Epoch_2), colour="#56B4E9", size=0.5, alpha=0.4) +
  geom_path(data = all_trials, aes(t_e,Epoch_3), colour="#009E73", size=0.5, alpha=0.4) +
  geom_path(data = all_trials, aes(t_e,Epoch_4), colour="#F0E442", size=0.5, alpha=0.4) +
  geom_path(data = all_trials, aes(t_e,Epoch_5), colour="#0072B2", size=0.5, alpha=0.4) +
  geom_path(data = all_trials, aes(t_e,Epoch_6), colour="#D55E00", size=0.5, alpha=0.4) +
  geom_path(data = all_trials, aes(t_e,Epoch_7), colour="#CC79A7", size=0.5, alpha=0.4) +
  geom_path(data = all_ave, aes(t_e, all_ave), colour="black", size=0.8, alpha=1) +
  xlim(0, 10) +
  ylim(-10, 10) +
  xlab("seconds") + ylab("uV") + 
  labs(title = sprintf("%s: All %s Trials and Trial Average", channel, trial_name), subtitle = "Spatially-filtered (Surface Laplacian), Bandpass-filtered (Mu 8-12 Hz), Epoch Length = 10s") +
  theme_cowplot()
  
# Spectral Bandpower Estimation (Squaring)

eeg_bp <- data.frame()
eeg_bp <- (all_trial_signals)^2
eeg_bp <- data.frame(eeg_bp)
eeg_bp <- cbind(t_e = all_trials$t_e, eeg_bp)
eeg_bp

sqd_ave <- data.frame()
sqd_ave <- rowMeans(eeg_bp)
sqd_ave <- data.frame(sqd_ave)
sqd_ave <- cbind(t_e = all_trials$t_e, Ave = sqd_ave)
sqd_ave

# Visualise Squared Signals (Spectral Bandpower)
ggplot() +
  geom_vline(aes(xintercept = 0, y = NULL, size = 0.5, alpha = 0.6), linetype = "dashed", show.legend = FALSE) +
  geom_path(data = eeg_bp, aes(t_e,Epoch_1), colour="#E69F00", size=0.5, alpha=0.2) +
  geom_path(data = eeg_bp, aes(t_e,Epoch_2), colour="#56B4E9", size=0.5, alpha=0.2) +
  geom_path(data = eeg_bp, aes(t_e,Epoch_3), colour="#009E73", size=0.5, alpha=0.2) +
  geom_path(data = eeg_bp, aes(t_e,Epoch_4), colour="#F0E442", size=0.5, alpha=0.2) +
  geom_path(data = eeg_bp, aes(t_e,Epoch_5), colour="#0072B2", size=0.5, alpha=0.2) +
  geom_path(data = eeg_bp, aes(t_e,Epoch_6), colour="#D55E00", size=0.5, alpha=0.2) +
  geom_path(data = eeg_bp, aes(t_e,Epoch_7), colour="#CC79A7", size=0.5, alpha=0.2) +
  geom_path(data = sqd_ave, aes(t_e,sqd_ave), colour="black", size=0.8, alpha=1) +
  xlim(0, 10) +
  ylim(0, 50) +
  xlab("seconds") + ylab("uV^2") + 
  labs(title = sprintf("Mu Band Power in %s Across All %s Trials", channel, trial_name), subtitle = expression(paste("Spatially-filtered (Surface Laplacian), Bandpass-filtered (", mu, " 8-12 Hz), Epoch Length = 10s"))) +
  theme_cowplot()

