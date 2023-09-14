## title: "Fit individual regressions and compare lines between groups"
## author: "Maria Ploumitsakou"
## date: "September 14, 2023"
## Email: maria.ploumitsakou@epfl.ch
## Purpose of script: fit regression lines on each subject's data and compare 
##                    fitted lines between the two subject groups

## ---------- install and load libraries -----------
if(!require(ggplot2)){install.packages("ggplot2")}
library(ggplot2)

if(!require(npmv)){install.packages("npmv")}
library(npmv)

## ---------- define data directory -----------

# IMPORTANT! You have to change the dataFolder according to where you placed 
# the my_project folder
dataFolder = "C:/Maria/PhD/courses/ORPER/ORPER_project/data/"

## ---------- read data and organize data frame -----------

# read betas from data file
betasFile <- paste(dataFolder, "sub-ALL_ses-02_task-PrePostCond_betas.tsv", 
                   sep = "")

betas <- read.csv(betasFile, header = T, sep = ",")

# the part_id variable is our subject list. Turn it into factor
betas$part_id <- as.factor(betas$part_id)

subList <- unique(betas$part_id)

# Group assignment. group ID = 1 for group 1 and 2 for group 2
groupID <- c(1,2,2,1,1,1,1,2,1,2,2,1,1,2,1,2,2,2,1,2,1,2,2,2,1,2,2,2,1,2,2,
             1,1,1,1,2,2,1,1,1,1)

# add group column to our data frame and turn it into a factor
betas$group <- rep(groupID,each=4)

betas$group <- as.factor(betas$group)

# code stimuli as continuous (first stimulus = 0 => meaningful zero, better for 
# intercept and slope calculation)
betas$stim <- 0
betas$stim[grep("GS1", betas$cond)] <- 1
betas$stim[grep("GS2", betas$cond)] <- 2
betas$stim[grep("CSminus_V", betas$cond)] <- 3

# define ROIs list and keep ROIs that are of interest to us
ROIs = unique(betas$anatROI)

ROIs = ROIs[c(1,4:9)]

#### ----------- Analysis - individually fitted regressions -----------

# define the coefficients data frame where individual slopes and 
# intercepts will be saved for each participant
columns = c("part_id", "ROI", "intercept", "slope", "group")

coefficients <- data.frame(matrix(nrow = length(subList)*(length(ROIs)), ncol = length(columns)))

colnames(coefficients) = columns

# repeat the subList as many times as the ROIs 
coefficients$part_id <- rep(subList, length(ROIs))

# add ROI column. Each ROI is repeated as many times as the number of subjects
coefficients$ROI <- rep(ROIs, each = length(subList))

# add group column. The group list is repeated as many times as the ROIs
coefficients$group <- rep(groupID, length(ROIs))

coefficients$group <- as.factor(coefficients$group)

# for every subject in the data frame, extract sub data
for (sub_i in unique(betas$part_id)) {
  currentSub_betas <- betas[betas$part_id == sub_i,]
  
  # for every ROI in this sub data, extract sub & ROI data
  for (ROI_i in unique(currentSub_betas$anatROI)){
    currentSub_currentROI_betas <- currentSub_betas[currentSub_betas$anatROI == ROI_i,]
    
    # fit LM explaining betas by stimulus
    currentLm <- lm(betas ~ stim, data = currentSub_currentROI_betas)
    
    # save LM slope and intercept into the corresponding cell of coefficients data frame
    coefficients$intercept[coefficients$part_id == sub_i &
                           coefficients$ROI == ROI_i] <- coef(currentLm)[1]
    
    coefficients$slope[coefficients$part_id == sub_i &
                       coefficients$ROI == ROI_i] <- coef(currentLm)[2]
  }
} 

# make boxplots of intercepts and slopes per group to visualize them
for (roi_i in ROIs) {
  print(roi_i)
  dataROI = coefficients[coefficients$ROI==as.character(roi_i),]
  
  print(ggplot(aes(group, slope), data = dataROI) + geom_boxplot(width = 0.1) +
          stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red") +
          stat_summary(fun.data = mean_se, geom = "errorbar", width=0.05) + theme_bw() +
          labs(title = paste(substr(as.character(roi_i),1,nchar(roi_i)-14), "Slopes")))
  
  print(ggplot(aes(group, intercept), data = dataROI) + geom_boxplot(width = 0.1) +
          stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red") +
          stat_summary(fun.data = mean_se, geom = "errorbar", width=0.05) + theme_bw() +
          labs(title = paste(substr(as.character(roi_i),1,nchar(roi_i)-14), "Intercepts")))
}

# separate each ROI's coefficients in a separate data frame and perform 
# a non parametric multivariate test to test if the lines of the two groups 
# are significantly different. 

# Important! When you run the nonpartest(), click on the console and hit enter
# to proceed to the next graph. When both graphs are plotted the statistical
# test outcome will appear on the console.

data_aINS = coefficients[coefficients$ROI=="aINS_space-EPI.nii",] # lines are different (singular system)
nonpartest(data_aINS$intercept|data_aINS$slope~data_aINS$group,data=data_aINS, permreps=1000, tests = c(1,1,1,1))

data_SMA = coefficients[coefficients$ROI=="SMA_space-EPI.nii",] # lines are different (singular system)
nonpartest(data_SMA$intercept|data_SMA$slope~data_SMA$group,data=data_SMA, permreps=1000, tests = c(1,1,1,1))

data_ACC = coefficients[coefficients$ROI=="ACC_space-EPI.nii",] # lines are different
nonpartest(data_ACC$intercept|data_ACC$slope~data_ACC$group,data=data_ACC, permreps=1000, tests = c(1,1,1,1))

data_PCC = coefficients[coefficients$ROI=="PrCu-PCC_space-EPI.nii",] # lines are not different (singular system)
nonpartest(data_PCC$intercept|data_PCC$slope~data_PCC$group,data=data_PCC, permreps=1000, tests = c(1,1,1,1))

data_THA = coefficients[coefficients$ROI=="THA_space-EPI.nii",] # lines are different (singular system)
nonpartest(data_THA$intercept|data_THA$slope~data_THA$group,data=data_THA, permreps=1000, tests = c(1,1,1,1))

data_CAU = coefficients[coefficients$ROI=="CAU_space-EPI.nii",] # lines are NOT different 
nonpartest(data_CAU$intercept|data_CAU$slope~data_CAU$group,data=data_CAU, permreps=1000, tests = c(1,1,1,1))

data_IFG = coefficients[coefficients$ROI=="IFG_space-EPI.nii",] # lines are different (singular system)
nonpartest(data_IFG$intercept|data_IFG$slope~data_IFG$group,data=data_IFG, permreps=1000, tests = c(1,1,1,1))
