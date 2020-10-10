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
time <- data.frame(time=-timesteps:(nrow(events)-timesteps-1))
raw_eeg <- cbind(time, events, raw_eeg)


# Epochs
trial_number <- data.frame(rep(1:40,each=2048))
colnames(trial_number) <- c("trial")
epochs_left <- dplyr::filter(raw_eeg, raw_eeg$events == 1)
epochs_left <- cbind(trial_number, epochs_left)

trial_number <- data.frame(trial=rep(1:39,each=2048))
excess <- data.frame(trial=rep(39,each=17))
trial_number <- rbind(trial_number, excess)
epochs_right <- dplyr::filter(raw_eeg, raw_eeg$events == -1)
epochs_right <- cbind(trial_number, epochs_right)

# Butterworth Filter
right_ftd <- eegfilter(epochs_right[4:19], Fs = 250, lower = 8.5, upper = 30.5, 
                       method = "butter", order = 4, forwardreverse = TRUE, scale = FALSE, plot = FALSE)
right_ftd <- cbind(trial=epochs_right$trial, right_ftd)
  
left_ftd <- eegfilter(epochs_left[4:19], Fs = 250, lower = 8.5, upper = 30.5, 
                       method = "butter", order = 4, forwardreverse = TRUE, scale = FALSE, plot = FALSE)
left_ftd <- cbind(trial=epochs_left$trial, left_ftd)
  


