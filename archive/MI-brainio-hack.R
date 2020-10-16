# MI analysis for stroke rehab, BR41N.IO Hackathon

# Import packages
library(R.matlab)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(eegkit)
library(Rmisc)

channel_name = "C3"
class1 = "Left"
class2 = "Right"

# Import file
filename = "E:/brain-io-hack/P1_pre_training.mat"
import_eeg <- readMat(filename)
events <- import_eeg$trig
raw_eeg <- import_eeg$y
colnames(raw_eeg) <- c("FC3", "FCz", "FC4", "C5", "C3", "C1", "Cz", "C2", 
                       "C4", "C6", "CP3", "CP1", "CPz", "CP2", "CP4", "Pz")

# Account for time shift
timesteps = as.integer(import_eeg$fs*8)
timestep <- data.frame(timestep=-timesteps:(nrow(events)-timesteps-1))
raw_eeg <- cbind(timestep, events, raw_eeg)


# Epochs
trial_number <- data.frame(rep(1:40,each=2048))
colnames(trial_number) <- c("trial")
epochs_left <- dplyr::filter(raw_eeg, raw_eeg$events == 1)
time_left <- data.frame(time=rep(1:2048, length.out=length(epochs_left$FC3)))
epochs_left <- cbind(trial_number, time_left, epochs_left) %>% select(-events)

trial_number <- data.frame(trial=rep(1:39,each=2048))
excess <- data.frame(trial=rep(39,each=17))
trial_number <- rbind(trial_number, excess)
epochs_right <- dplyr::filter(raw_eeg, raw_eeg$events == -1)
time_right <- data.frame(time=rep(1:2048, length.out=length(epochs_right$FC3)))
epochs_right <- cbind(trial_number, time_right, epochs_right) %>% select(-events)



# Notch Filter
left_ftd <- eegfilter(epochs_left[4:19], Fs = 250, lower = 0.1, upper = 49.5, 
                      method = "fir1") 
right_ftd <- eegfilter(epochs_right[4:19], Fs = 250, lower = 0.1, upper = 49.5, 
                       method = "fir1")


# Butterworth Filter
left_bp <- (eegfilter(left_ftd, Fs = 250, lower = 8.5, upper = 30.5, 
                      method = "butter", order = 4, forwardreverse = TRUE, scale = FALSE, plot = FALSE))^2 
left_bp <- cbind(trial=epochs_left$trial, time_left, left_bp)

right_bp <- (eegfilter(right_ftd, Fs = 250, lower = 8.5, upper = 30.5, 
                       method = "butter", order = 4, forwardreverse = TRUE, scale = FALSE, plot = FALSE))^2
right_bp <- cbind(trial=epochs_right$trial, time_right, right_bp)

# Trial Averaging
left_bp_ave <- aggregate(left_bp[,3:18], list(left_bp$time), mean)
right_bp_ave <- aggregate(right_bp[,3:18], list(right_bp$time), mean)

# Moving average

Legend <- c("Left Trials" = "blue", "Right Trials" = "red")

ggplot() +
  geom_path(data = left_bp_ave, aes(Group.1,C4, color="Left Trials"), alpha=0.5, size = 1) +
  geom_path(data = right_bp_ave, aes(Group.1,C4, color="Right Trials"), alpha=0.5, size = 1) + 
  xlab('Samples') + ylab('u + B bandpower') + 
  labs(title = "C3 Bandpower Estimation in u+B (8-30 Hz)") +
  scale_color_manual(values = Legend)
  


