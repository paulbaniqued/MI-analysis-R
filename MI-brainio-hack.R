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
