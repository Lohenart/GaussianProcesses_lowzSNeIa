# GAUSSIAN PROCESS CODE
# Author: Arturo Avelino

#     DESCRIPTION
# - The hyperparameters are computed using all the LCs but making qualities cutoff  
# based on dm15, z_cmb, EBVhost, EBV_mw.
# - Once the hyperparameters are determined I fit -all- (i.e., no restrictions in the qualities cutoff) the LCs using those hyperparameters
# - The covariance matrix is computed taking into account the uncertainty in the peculiar velocity (default, by setting UsePecMatrix <- TRUE), velPecuFix, or ignoring velPecuFix by setting UsePecMatrix <- FALSE.
# - The mean GP function is computed -ignoring- velPecuFix (default, by setting ComputeMeanUsingPecMatrix <- FALSE), but can be instead used if ComputeMeanUsingPecMatrix <- TRUE and UsePecMatrix <- TRUE.


# For more information about the theory see Chapter 2 of Rasmussen and Williams's 
# book `Gaussian Processes for Machine Learning' provides a detailed explanation of the
# math for Gaussian process regression. 

# Although the nested loops are ugly, I've checked and it's about
# 30% faster than a solution using expand.grid() and apply()

# Load in the required libraries for data manipulation
# and multivariate normal distribution

# Version 3.0 (Sep 27, 2016)
# - implementation of the previous template as a PRIOR for the GP fitting 

####################################################
require(MASS)
require(plyr)
require(reshape2)
require(ggplot2)

#     USER DEFINITIONS

#-- Choosing the sample to interpolate using the GP
# sample <- 'CfA'
# sample <- 'CSP'
# sample <- 'Others'
sample <- 'AllSamples'

#-----------

# Compute the GP hyperparameters from the global likelihood PDF ('ComputeHyperpars <- TRUE')?
# or assume a fixed value set by hand ('ComputeHyperpars <- FALSE')?
# I use 'TRUE' for the paper, however, when I need to remake the GP fit for some individual LCs, then I use 'FALSE' and then set up the values of hyperparameters by hand.
ComputeHyperpars <- FALSE  # (TRUE, FALSE). 'TRUE' is the option that I use for the paper

# Peculiar velocity uncertainty (km/sec)
velPecuFix <- 0  # (0, 150, 300) to compute the template. Use 0 (zero) to compute the -normalized- template.
# velPecuFix = 0 # To compute the normalized template.


# Cutoffs on the SNe used to determine the hyperparameters of the GP kernel

zMin <- 0
zMax <- 0.09

delta15Min <- 0.78
delta15Max <- 1.6 #
EBVhostMin <- -0.4
EBVhostMax <- 0.4
EBVmwMax <- 1

MinNumDataForTempl <- 3 # Minimal number of data points in a given LC to be use for the template

#--------------------------
#-- Fitting using a 'flat' or 'template' as prior

# GP_prior = 'flat': assuming a abs mag value for GP  mean function prior, specified by 'meanPriorFix'. We will use this method for the paper.
# GP_prior = 'template': assuming an existing LC template as the GP  mean function prior, usually determined by a moving window average method.

GP_prior <- 'flat' # I use this option for the paper.
# GP_prior <- 'template'

# Location of the template (run this line even if I'm using a flat prior):
# DirTemplate <- 'Std_filters/3_Template/AllSamples/z_gr_001'
DirTemplate <- 'Std_filters/3_Template_MWA/AllSamples/z_gr_001'
TemplateName <- 'TempWeightedSmooth_Box7_Step05_Window21_Poly3_.dat'

#-----------

# Use the peculiar velocity covariance matrix in the 'k.xx' to compute the mean function?:
ComputeMeanUsingPecMatrix <- FALSE # (TRUE, FALSE). 'FALSE' is the way I proceed for the paper.
# OLD. MeanTrick <- TRUE # TRUE, FALSE

# Use the peculiar velocity covariance matrix in the k.xsx, k.xxs, kxsxs terms?:
# Note that this option does NOT have any effect on the k.xx term.
UsePecMatrix <- FALSE # (TRUE, FALSE). 'FALSE' is the CORRECT way and it is what I use for the paper.
# About this option. Initially I didn't include the peculiar velocity cov matrix on the k.xsx, k.xxs, kxsxs terms, only on k.xx, however then Kaisey suggested that the peculiar velocity cov matrix should also be included on those terms (i.e., k.xsx, k.xxs, kxsxs). So I did several tests and I realized that I obtained very narrow 'snakes' for each light curve fit. After discussing again with Kaisey we realized that the peculiar velocity cov matrix should NOT be included in k.xsx, k.xxs, kxsxs terms, only on k.xx, in this way I obtain snakes that make sense, i.e., the width of the snake is no narrower than the peculiar velocity uncertainty for each SN. So, this option was implemented to easily explore Kaisey's idea, however, the CORRECT option for this variable MUST be UsePecMatrix <- FALSE.

#--------------------------

# Step size of the test points
StepSizeTestPoints <- 0.5 # in days.

# The limits to create the filled files of data to be used in
# the hierarchical code.
phase_lowerLimit <- -35
phase_upperLimit <- 60

# Set a seed for repeatable plots
randomSeed <- 12345
set.seed(randomSeed)

#------------------------------------------

#     INFORMATION FOR EACH BAND


#     J BAND   (OK)
bandname = 'J'

#- Mean prior for the Gaussian processes
# I shift the data in y-axis to around y=0 with the help of 'meanPriorFix' in order to fit the GP as it were around zero, otherwise, the GP fit set some values around y=0 but if the data is in its actual y-axis values (~ -18 mag) then it produces a very large and bias 'snakes' in the regions with no data.
meanPriorFix =  -17 # flat GP prior  

#-- General values for fitting --
InitialGuess_length <- 10 # 16 # old 14  
InitialGuess_sigma2kern <- 1

LowerLimit_length <- 5 # 14.5 # OLD: 9.5, 11.5, 12.5 
UpperLimit_length <- 30 # 25.5

LowerLimit_sigma2kern <- 0.01
UpperLimit_sigma2kern <- 20

# If ComputeHyperpars = FALSE, then what values set by hand I use for (l, sigma)?
# Note: I can run this line no matter the decision about 'ComputeHyperpars'.
if (velPecuFix==0) {
  HyperByHand <- c(7.5352, 1.0029) # when velPecuFix <- 0 km/s
} else if (velPecuFix==150){
  HyperByHand <- c(7.4790, 0.9918) # when velPecuFix <- 150 km/s 
} else if (velPecuFix==300) {
  HyperByHand <- c(7.4343, 0.9726) # when velPecuFix <- 300 km/s
}

ymin_plot <- -14
ymax_plot <- -21

# Location of text in the plot
xText <- c(-10, -5, 0, 2, -10, -1.5, -10, -5, -10, 2, -10, 2, 8.5, 10.5, 35)
yText <- c(-15.8, -15.8, -15.8, -15.8, -15.4, -15.4, -15, -15, -14.6, -14.6, -14.2, -14.2, -14.2, -14.2, -19.5)


#------------------------------------------

"
#     Y BAND  (Ok)

bandname = 'Y'

#- Mean prior for the Gaussian processes
# I shift the data in y-axis to around y=0 with the help of 'meanPriorFix' in order to fit the GP as it were around zero, otherwise, the GP fit set some values around y=0 but if the data is in its actual y-axis values (~ -18 mag) then it produces a very large and bias 'snakes' in the regions with no data.
meanPriorFix = -17.5 # old: 18 

#-- General values for fitting --
InitialGuess_length <- 12 # # old 14  
InitialGuess_sigma2kern <- 1
LowerLimit_length <- 4 # 14 # OLD: 9.5, 13
UpperLimit_length <- 40 # 24.5 # OLD: 20, 24

LowerLimit_sigma2kern <- 0.1
UpperLimit_sigma2kern <- 20

# If ComputeHyperpars = FALSE, then what values set by hand I use for (l, sigma)?
# Note: I can run this line no matter the decision about 'ComputeHyperpars'.

if (velPecuFix==0){
  HyperByHand <- c(8.5952, 0.7428) # when velPecuFix <- 0 km/s 
} else if (velPecuFix==150){
  HyperByHand <- c(8.5175, 0.7346) # when velPecuFix <- 150 km/s 
} else if (velPecuFix==300) {
  HyperByHand <- c(8.4647, 0.7224) # when velPecuFix <- 300 km/s
}

ymin_plot <- -14.5 # tmp: -14
ymax_plot <- -19.5  # tmp: -21

# Location of text in the plot
xText <- c(-10, -5, 0, 2, -10, -1.5, -10, -5, -10, 2, -10, 2, 8.5, 10.5, 35)
yText <- c(-17, -17, -17, -17, -16.7, -16.7, -16.4, -16.4, -16.1, -16.1, -15.8, -15.8, -15.8, -15.8, -19)
"

#------------------------------------------

"
#     H BAND  (OK)

bandname = 'H'

#- Mean prior for the Gaussian processes
# I shift the data in y-axis to around y=0 with the help of 'meanPriorFix' in order to fit the GP as it were around zero, otherwise, the GP fit set some values around y=0 but if the data is in its actual y-axis values (~ -18 mag) then it produces a very large and bias 'snakes' in the regions with no data.
meanPriorFix = -18 # old: 15 

#-- General values for fitting --
InitialGuess_length <- 10
InitialGuess_sigma2kern <- 1
LowerLimit_length <- 5 # 13.5 <-- original, Ok. 
UpperLimit_length <- 30 # 24.5 <-- original, Ok.

LowerLimit_sigma2kern <- 0.3
UpperLimit_sigma2kern <- 20

# If ComputeHyperpars = FALSE, then what values set by hand I use for (l, sigma)?
# Note: I can run this line no matter the decision about 'ComputeHyperpars'.
if (velPecuFix==0){
  HyperByHand <- c(9.3354, 0.8692) # when velPecuFix <- 0 km/s 
} else if (velPecuFix==150){
  HyperByHand <- c(9.1412, 0.8227) # when velPecuFix <- 150 km/s 
} else if (velPecuFix==300) {
  HyperByHand <- c(8.9101, 0.7666) # when velPecuFix <- 300 km/s
}

ymin_plot <- -14.5
ymax_plot <- -20

# Location of text in the plot
xText <- c(-10, -5, 0, 2, -10, -1.5, -10, -5, -10, 2, -10, 2, 8.5, 10.5, 35)
yText <- c(-17, -17, -17, -17, -16.7, -16.7, -16.4, -16.4, -16.1, -16.1, -15.8, -15.8, -15.8, -15.8, -19)
"

#------------------------------------------

"
#     K BAND

bandname = 'K'

# Mean prior for the Gaussian processes
meanPriorFix = -18 

#-- General values for fitting --
InitialGuess_length <- 16
InitialGuess_sigma2kern <- 1
LowerLimit_length <- 5 # 16 <-- original, Ok. 
UpperLimit_length <- 40 # 18 <-- original, Ok.

LowerLimit_sigma2kern <- 0.1
UpperLimit_sigma2kern <- 20

# If ComputeHyperpars = FALSE, then what values set by hand I use for (l, sigma)?
# Note: I can run this line no matter the decision about 'ComputeHyperpars'.
if (velPecuFix==0){
  HyperByHand <- c(8.1163, 0.6517) # when velPecuFix <- 0 km/s 
} else if (velPecuFix==150){
  HyperByHand <- c(8.0240, 0.6533) # when velPecuFix <- 150 km/s 
} else if (velPecuFix==300) {
HyperByHand <- c(7.9214, 0.6325) # when velPecuFix <- 300 km/s
}

ymin_plot <- -14.5
ymax_plot <- -20.5

# Location of text in the plot
xText <- c(-10, -5, 0, 2, -10, -1.5, -10, -5, -10, 2, -10, 2, 8.5, 10.5, 35)
yText <- c(-17, -17, -17, -17, -16.7, -16.7, -16.4, -16.4, -16.1, -16.1, -15.8, -15.8, -15.8, -15.8, -19)
"

#------------------------------------------

#   END OF USER SUPPLIED INFORMATION. THE REST SHOULD BE AUTOMATICALLY DONE FOR ALL THE BANDS WITHOUT THE NEED OF PARTICULAR MODIFICATIONS

#############################################################
#############################################################
#     AUTOMATIC

cc = 299792.458  # Speed of light
# OLD. zerrorFix = 0.001 # Error in the estimation of redshift

# Number of lines of only information about a given LC in the text file: 
NumLinesInfoSN <- 8 

# Minimal number of data in a LC to fit it. It doesn't necessary mean that
# it will be used for the template or Hubble diagram.
MinNumDataToFitIt <- 3 

#------------------------------------------
#     Uploading the data

# The main directory containing the subfolders with data.
MainDir <- paste('/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/TheTemplates/',bandname,'_band', sep='')
MainDir

KindOfData <- paste('Std_filters/1_AllData_InitialFit/AbsMag/', sample, sep='')
KindOfData

#- Concatenate the paths
DirLCData <- file.path(MainDir, KindOfData)
DirLCData

# OLD. DirSaveOutput <- file.path(DirLCData,'Goods/')
# dir.create(DirSaveOutput)

DirSaveOutput <- DirLCData
DirSaveOutput

# Old. DirSaveAppmag <- file.path(MainDir,'Std_filters/1_AllData_InitialFit/AppMag/')
# Old. DirSaveAppmag

# Set the final directory where the data is located
setwd(DirLCData)
getwd() # Show the current directory

#--  Read the LC datafile names and create an array with them
FinalNameEachFile <- paste('*_',bandname,'.txt', sep='')
FinalNameEachFile
list_SNe <- list.files(pattern=FinalNameEachFile, all.files=TRUE) 
list_SNe
list_SNe[1]

numSNe <- length(list_SNe)
numSNe

#-----------
# Function to substract the first letters to the name of the SN file
# http://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r

substrRight <- function(x, n){
  substr(x, 1, nchar(x)-n)
}

#############################################################

#    WRITE TEXT FILES WITH INFORMATION ABOUT THE SETTINGS AND LISTING THE SNe

# Create a label for 'HyperByHand':
if (ComputeHyperpars == TRUE){ HandLabel <- FALSE
} else{ HandLabel <- TRUE}

text0 <- '# Supernova after quality cutoffs' 
text1  <- sprintf('# Band: %s', bandname)
text1_1  <- sprintf('# SN sample to interpolate using the GP: %s', sample)
text2 <- '# They have been selected based on the following criteria:'
text3 <- sprintf('# %.2f < zcmb < %.2f, %.2f < dm15 < %.2f,  %.2f < EBVhost < %.2f, EBVmw < %.2f', zMin, zMax, delta15Min, delta15Max, EBVhostMin, EBVhostMax, EBVmwMax)
text3_1 <- sprintf('# Number of LC data points >= %.f', MinNumDataForTempl)

text4 <- '# Location of the LC files used:'
text5 <- KindOfData
textLine <- '#-------------------------------------'
text7 <- '#------- Settings of the GP fitting--------'
text7_0_1 <- sprintf('# GP_prior = %s', GP_prior)
text7_0_2 <- sprintf('# If GP_prior=flat, then the assumed abs mag flat value is = %.2f', meanPriorFix)
text7_1 <- sprintf('# Compute the GP hyperparameters from the global likelihood PDF?: %s',ComputeHyperpars)
text7_2 <- sprintf('# Or, assume fixed values set by hand for the kernel parameters?: %s',HandLabel)
text9 <- sprintf('# Hyperparameters set by hand if TRUE: (l = %.4f, sigma_kern = %.4f)', HyperByHand[1], HyperByHand[2])
text11 <- sprintf('# Random number seed = %.f', randomSeed)
text11_1 <- sprintf('# Use the peculiar velocity covariance matrix to determine the -mean- GF function?:%s', ComputeMeanUsingPecMatrix)
text11_2 <- sprintf('# Note that to compute the GP hyperparameters from the global likelihood PDF it IS used the peculiar velocity covariance matrix (it is the correct way to proceed).')
text11_3 <- sprintf('# Use the peculiar velocity covariance matrix in the k.xsx, k.xxs, kxsxs terms...')
text11_4 <- sprintf('#... to compute the covariance of the individual GP fit of each light curves (the correct option MUST be FALSE)?: %s', UsePecMatrix)
text11_5 <- sprintf('# Peculiar velocity assumed = %.1f km/sec', velPecuFix)
text12<-sprintf('#---- Settings to compute the hyperparameters if TRUE ----') 
text13<-sprintf('# InitialGuess_length = %.2f', InitialGuess_length) 
text14<-sprintf('# InitialGuess_sigma2kern = %.2f', InitialGuess_sigma2kern) 
text15<-sprintf('# LowerLimit_length = %.2f', LowerLimit_length ) 
text16<-sprintf('# UpperLimit_length = %.2f', UpperLimit_length) 
text17<-sprintf('# LowerLimit_sigma2kern = %.2f', LowerLimit_sigma2kern) 
text18<-sprintf('# UpperLimit_sigma2kern = %.2f', UpperLimit_sigma2kern) 

#------------------------------

#             'List_SN_afterCutoffs_.txt'

write.table(text0, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(text1,   file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text1_1, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text2, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text3,   file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text3_1, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text4, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text5, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names='#',   col.names=FALSE, append=TRUE)
write.table('# ', file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text7,     file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text7_0_1, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
if (GP_prior == 'flat'){
write.table(text7_0_2, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE) }
if (GP_prior == 'template'){
write.table(DirTemplate, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names='# DirTemplate =', col.names=FALSE, append=TRUE)
write.table(TemplateName, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names='# TemplateName =', col.names=FALSE, append=TRUE) }
write.table(text11, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text11_1, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text11_2, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text11_3, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text11_4, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text11_5, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text7_1, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text7_2, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
if (ComputeHyperpars == FALSE) {
  write.table(text9, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)  }

if (ComputeHyperpars == TRUE){
write.table(text12, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text13, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text14, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text15, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text16, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text17, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text18, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)  }

write.table(textLine, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

#------------------------------

#            'Settings_GPFit_.txt'

write.table(text0, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(text1,   file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text1_1, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text2, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text3,   file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text3_1, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text4, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text5, file='Settings_GPFit_.txt', quote=FALSE, row.names='#',   col.names=FALSE, append=TRUE)
write.table('# ', file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text7,     file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text7_0_1, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
if (GP_prior == 'flat'){
  write.table(text7_0_2, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE) }
if (GP_prior == 'template'){
  write.table(DirTemplate, file='Settings_GPFit_.txt', quote=FALSE, row.names='# DirTemplate =', col.names=FALSE, append=TRUE)
  write.table(TemplateName, file='Settings_GPFit_.txt', quote=FALSE, row.names='# TemplateName =', col.names=FALSE, append=TRUE) }
write.table(text11, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text11_1, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text11_2, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text11_3, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text11_4, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text11_5, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text7_1, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(text7_2, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
if (ComputeHyperpars == FALSE) {
  write.table(text9, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)  }

if (ComputeHyperpars == TRUE){
  write.table(text12, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(text13, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(text14, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(text15, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(text16, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(text17, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(text18, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)  }

write.table('# ', file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)


#############################################################

#     CREATE A TEXTFILE LIST OF SNe AFTER CUTOFFS
# These will be used to determine the values of they hyperparameters (l, sigma) of the Gaussian Process kernel

count1 <- 0

for(i in 1:numSNe){
  DataLC = read.table(list_SNe[i])
  # DataLC = read.table(list_SNeAfterCutoff[1])
  
  redshiftSN <- DataLC$V1[1] # Redshift_CMB
  zcmb_error <- DataLC$V2[1] # z_CMB_error
  delta15       <- DataLC$V1[2] # Delta_15 parameter
  delta15_error <- DataLC$V2[2] # Delta_15 parameter error
  EBVmw         <- DataLC$V1[4]
  EBVmw_err     <- DataLC$V2[4]
  EBVhost       <- DataLC$V1[5]
  EBVhost_err   <- DataLC$V2[5]
  NumLCData <- dim(DataLC)[1]-NumLinesInfoSN # Number of data points in the LC

  if (redshiftSN > zMin & redshiftSN < zMax & delta15>delta15Min & delta15<delta15Max & EBVhost>EBVhostMin & EBVhost < EBVhostMax & EBVmw < EBVmwMax & NumLCData >= MinNumDataForTempl) {
    # Save the name of the SNe used in a text file
    write.table(list_SNe[i], file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
    count1 <- count1 + 1
  }
}

text6 <- sprintf('# Number of SNe in this list: %.f', count1)
write.table(text6, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

# Number of SNe after cutoffs:
count1


#------------------------------------------------------------

#   READ THE TEXT FILE OF THE SNe LIST TO BE USED TO DETERMINE THE HYPERPARAMETERS
# These will be used to determine the values of they hyperparameters (l, sigma) of the Gaussian Process kernel

# First check if I have one file with comments about the SNe; if so read that file 
# otherwise read the just created file
try1 <- try(list_SNeAfterCutoff <- read.table('List_SN_afterCutoffs_Notes_.txt'))
if(inherits(try1, "try-error")) {
  print('List_SN_afterCutoffs_Notes_.txt not found.')
  list_SNeAfterCutoff <- read.table('List_SN_afterCutoffs_.txt')
  # writeLines('Reading data from List_SN_afterCutoffs_.txt instead')
  print('Reading data from List_SN_afterCutoffs_.txt instead')
}

list_SNeAfterCutoff$V1[1] # show the first entry
# list_SNeAfterCutoff
# class(list_SNeAfterCutoff)

numSNeAfterCutoff <- dim(list_SNeAfterCutoff)[1]
numSNeAfterCutoff


#############################################################


#   INTERPOLATION OF THE PRIOR TEMPLATE
# The interpolation is used to compute the "residual" data before to fit it
# using the Gaussian Process method.
#  Interpolate the template when using it as a prior for the GP fitting

if (GP_prior == 'template'){
  #-- Reading an existing template file ---
  # TempData <- read.table(file.path(MainDir, DirTemplate, 'TempWeightedSmooth_Mean_StdError_Variance_Box3_Step05.dat'))
  TempData <- read.table(file.path(MainDir, DirTemplate, TemplateName))
  # TempData <- read.table(file.path(MainDir, DirTemplate, 'Template_phase_mu_tau_FromR_test.dat'))
  # TempData
  # TempData$V1
  # class(TempData)
  
  #- Interpolating the template
  TempInterpol <- approxfun(TempData$V1, TempData$V2)
  TempInterpol(0.3) # Test
  
  InterpolMin <- TempData$V1[1]
  InterpolMin
  InterpolMax <- TempData$V1[length(TempData$V1)]
  InterpolMax
}


#############################################################

#     KERNEL DEFINITION
# Function to compute the covariances matrices
# It uses the covariance function (i.e., kernel) of Eq. (2.31) of 'Rasmussen_William_2006_GaussianProcessesForMachineLearning_N.pdf'

kernelFunc <- function(X1,X2,l, sigma2kern) {
  # Empty matrix K with zeros as entries
  K <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(K)) {
    for (j in 1:ncol(K)) {
      # The covariance function Eq. (2.31) without the last term:
      K[i,j] <- sigma2kern*exp(-0.5*((X1[i]-X2[j])/l)^2)
    }
  }
return(K)
}


# Testing that the function works well.
# kernelFunc(9.9464906258 , 9.9464906258, 5.6, 0.4^2)
# 0.16


#----------------------------------------------------

# The peculiar velocity variance (standard error^2) function.
# Using Eq. (1) of Wood-Vasey+2008 paper.

sigma2muPecu <- function(vpec, z, zerror){
  cc = 299792.458  # Speed of light in km/s
  ( (5/(z*log(10))) * sqrt((vpec/cc)**2 + zerror**2) )^2
}

#############################################################


#   GLOBAL LOG MARGINAL LIKELIHOOOD FUNCTION
# Using ALL the light curves simultaneously to determine the GP hyperparameters

LogMargLikelFuncAll <- function(hyperpars) {
  
  # This is a function of 2 variables: (l, sigma2kern). 
  l <- hyperpars[1] # "l" in Eq. (2.31).
  sigma2kern <- hyperpars[2] # "sigma^2_f" in Eq. (2.31).
  
  logMarLikelSingle <- 0
  
  # LOOP over the supernovae (i.e, the light curves)
  for(i in 1:numSNeAfterCutoff){
    # for(i in 1:42){
    # i<-3
    
    # list_SNeAfterCutoff$V1[i]
    DataLC <- read.table(paste('', list_SNeAfterCutoff$V1[i],'', sep = ''))
    # DataLC
    
    # SNname <- list_SNeAfterCutoff$V1[i]
    # SNname
    
    LengthData <- dim(DataLC)[1]
    # LengthData
    
    redshiftSN <- DataLC$V1[1] # z_CMB
    # redshiftSN
    zcmb_error <- DataLC$V2[1] # z_CMB_error
    
    # The "training set", i.e., the observed data.
    # ASSUMING A MEAN PRIOR
    # I shift the data in y-axis (abs mag) to around y=0 with the help of 'meanPriorFix' in order to fit the GP as if the abs mag data were around zero, otherwise, the GP fit set some values around y=0 but if the data is in its actual y-axis values (~ -18 mag) then it produces a very large and bias "snakes" in the regions with no data.
    if (GP_prior == 'flat'){
      f <- data.frame(x=DataLC$V1[9:LengthData], y=(DataLC$V2[9:LengthData] - meanPriorFix), z=DataLC$V3[9:LengthData])
      counter2 <- LengthData
      # print('Flat prior is used.')
      
    } else if (GP_prior == 'template'){
      
      x_array <- numeric(LengthData) # Creating the array for x.
      # length(x_array)
      residual <- numeric(LengthData) # Creating the array.
      sigma_array <- numeric(LengthData)
      counter2 <- 0; counter3 <- 0
      
      for(j in 9:LengthData){
        if(DataLC$V1[j] > InterpolMin){
          if(DataLC$V1[j] < InterpolMax){
            x_array[j-8] <- DataLC$V1[j]
            residual[j-8] <- DataLC$V2[j] - TempInterpol(DataLC$V1[j])
            sigma_array[j-8] <- DataLC$V3[j]
            counter2 <- counter2 + 1 }
        } else { counter3 <- counter3 + 1 } # Count number of cases where DataLC$V1[j] < InterpolMin
      }
      # Create the data frame:
      f <- data.frame(x=x_array[(counter3+1):(counter3+counter2)], y=residual[(counter3+1):(counter3+counter2)], z=sigma_array[(counter3+1):(counter3+counter2)])
    }
    
    # f
    # plot(f$x, f$y, ylim=c(1.5,-1.5)) 
    # dev.off()
    
    # Assigning x as f$x
    x <- f$x
    # x
    
    # Lenght trainning set
    nn <- length(f$x)
    # nn
    
    # The standard deviation of the noise
    # sigma.n <- DataLC$V3[9:LengthData] # old
    sigma.n <- f$z
    # sigma.n
    # length(sigma.n)
    
    # PECULIAR VELOCITY COVARIANCE MATRIX
    # OLD. sigma2PecuMatrix_xx <- matrix(sigma2muPecu(velPecuFix, redshiftSN, zerrorFix), nn, nn)
    sigma2PecuMatrix_xx <- matrix(sigma2muPecu(velPecuFix, redshiftSN, zcmb_error), nn, nn)
    
    # sigma2PecuMatrix_xx
    # dim(sigma2PecuMatrix_xx)
    # sqrt(length(sigma2PecuMatrix_xx))
    
    # Covariance matrix  K(X, X). It is upper left term in matrix Eq. (2.21)
    # k_xx <- kernelFunc(x,x, 12, 0.5)
    k_xx <- kernelFunc(x,x, l, sigma2kern)
    # dim(k_xx)
    
    # Eq. (2.30)
    # The negative of the log marginal likehood.
    
    # WITH THE COVARIANCE MATRIX OF PECULIAR VELOCITY CORRELATION. ORIGINAL AND CORRECT WAY TO COMPUTE
    logMarLikelSingle <- logMarLikelSingle - (-0.5*t(f$y)%*%solve(k_xx + (sigma.n^2)*diag(1,ncol(k_xx)) + sigma2PecuMatrix_xx)%*%f$y - 0.5*log(det(k_xx + (sigma.n^2)*diag(1,ncol(k_xx)) + sigma2PecuMatrix_xx)) - 0.5*nn*log(2*pi) )[1]
    
    # WITHOUT THE COVARIANCE MATRIX OF PECULIAR VELOCITY CORRELATION
    # logMarLikelSingle <- logMarLikelSingle -(-0.5*t(f$y)%*%solve(k_xx + (sigma.n^2)*diag(1,ncol(k_xx)))%*%f$y - 0.5*log(det(k_xx + (sigma.n^2)*diag(1,ncol(k_xx)))) - 0.5*nn*log(2*pi) )[1]
  }
  logMarLikelSingle # Report the final sum 
} # End definition of the global log marginal likelihood

# test1 <- c(10,0.6)
# LogMargLikelFuncAll(test1)
# -151.3471


#---------------

#     DEFINE VALUES OF GP HYPERPARAMETERS

if (ComputeHyperpars == TRUE){ # SEARCHING THE MINIMUM SETTING LIMITS TO THE PARAMETERS
  ResultOptim <- optim(c(InitialGuess_length, InitialGuess_sigma2kern), 
                       LogMargLikelFuncAll, method='L-BFGS-B', 
                       lower = c(LowerLimit_length, LowerLimit_sigma2kern), 
                       upper=c(UpperLimit_length, UpperLimit_sigma2kern))
  l_Fix <- ResultOptim$par[1]  
  sigma2kern_fix <- ResultOptim$par[2] 
  
  text20 <- sprintf('# Best estimated value for l = %.4f', l_Fix)
  text21 <- sprintf('# Best estimated value for sigma_kern = %.4f', sqrt(sigma2kern_fix))
  write.table(text20, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(text21, file='List_SN_afterCutoffs_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(text20, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(text21, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(textLine, file='Settings_GPFit_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
} else if (ComputeHyperpars == FALSE) { # OR ASSIGN A SPECIFIC VALUE SET BY HAND
  l_Fix <- HyperByHand[1]  
  sigma2kern_fix <- (HyperByHand[2])^2 
}

# ResultOptim$par

# 6.1702528 0.1115946 # Using the correct expression (i.e., with the pec vel uncertainty)
# 6.5984559 0.1452227 # Without the pec vel matrix
# 5.8836066 0.1060577 # Using a template from GP, with pec. vel. matrix
# 6.2409010 0.1382619 # Using a template from GP, without pec. vel. matrix
# 6.0288230 0.1137106 # with pec. vel. matrix

# ResultOptim$value
# -240.2444 # Using the correct expression (i.e., with the pec vel uncertainty)
# -230.0351 # Without the pec vel matrix
# -266.2605 # Using a template from GP, with pec. vel. matrix
# -255.9029 # Using a template from GP, without pec. vel. matrix
# -243.3802 # with pec. vel. matrix

# l_Fix <- 6.170253
# l_Fix <- 5.698
# l_Fix <- ResultOptim$par[1]  
# l_Fix  # The kernel's lenght scale best estimate
# 11.31701
# 6.170253 # Applying cutoffs


# sigma2kern_fix <- 0.1115946
# sigma2kern_fix <- 0.318
# sigma2kern_fix <- 2
# sigma2kern_fix <- ResultOptim$par[2] 
# sigma2kern_fix # The kernel's variance best estimate
# 0.6404094
# 0.1115946

# sigma2kern_fix^(0.5)

# LogMargLikelFuncAll(ResultOptim$par)
# 32.91491


#############################################################
#############################################################




#     LOOP to fit the individual light curves but using the hyperparameter values computed above.
# The loop is over the good SN LCs after the quality cutoffs.

for(i in 1:numSNe){
# for(i in 1:10){
# i<-1
# for(i in 1:numSNeAfterCutoff){

  # Reset some important values
  k.xx <- 0; k.xx_mean <- 0
  
  # OLD. DataLC <- read.table(paste('', list_SNeAfterCutoff$V1[i],'', sep = ''))
  DataLC = read.table(list_SNe[i])
  # DataLC
  # list_SNe[i]
  
  LengthData <- dim(DataLC)[1]
  # LengthData
  
  # Fit the LC with at least (MinNumDataForTempl-1) data points:
  if ((LengthData-NumLinesInfoSN) >= MinNumDataToFitIt) { 
    
    # Defined the name of the SN by the name of its data file
    # Removing the extension ".txt" to the filename.
    NameSN <- substrRight(list_SNe[i], 4)
    # NameSN    
    
    redshiftSN <- DataLC$V1[1] # z_CMB
    # redshiftSN
    zcmb_error <- DataLC$V2[1] # z_CMB_error
    
    # This information is just to be printed in plots
    delta15       <- DataLC$V1[2] # Delta_15 parameter
    delta15_error <- DataLC$V2[2] # Delta_15 parameter error
    DistanceMu    <- DataLC$V1[3] # Distance modulus from Snoopy
    DistanceMu_error <- DataLC$V2[3] # Distance modulus from Snoopy

    # The "training set", i.e., the observed data.
    # ASSUMING A MEAN PRIOR
    # I shift the data in y-axis (abs mag) to around y=0 with the help of 'meanPriorFix' in order to fit the GP as if the abs mag data were around zero, otherwise, the GP fit set some values around y=0 but if the data is in its actual y-axis values (~ -18 mag) then it produces a very large and bias "snakes" in the regions with no data.
    if (GP_prior == 'flat'){
      f <- data.frame(x=DataLC$V1[9:LengthData], y=(DataLC$V2[9:LengthData]-meanPriorFix), z=DataLC$V3[9:LengthData])
      counter2 <- LengthData
      # print('Flat prior is used.')
      
    } else if (GP_prior == 'template'){
      
      x_array <- numeric(LengthData) # Creating the array for x.
      # length(x_array)
      residual <- numeric(LengthData) # Creating the array.
      sigma_array <- numeric(LengthData) # array for the data uncertainties
      counter2 <- 0; counter3 <- 0
      for(j in 9:LengthData){
        if(DataLC$V1[j] > InterpolMin){
          if(DataLC$V1[j] < InterpolMax){
            x_array[j-8] <- DataLC$V1[j]
            residual[j-8] <- DataLC$V2[j] - TempInterpol(DataLC$V1[j])
            sigma_array[j-8] <- DataLC$V3[j]
            counter2 <- counter2 + 1 }
        } else { counter3 <- counter3 + 1 } # Count number of cases where DataLC$V1[j] < InterpolMin
      }
      # Create the data frame:
      f <- data.frame(x=x_array[(counter3+1):(counter3+counter2)], y=residual[(counter3+1):(counter3+counter2)], z=sigma_array[(counter3+1):(counter3+counter2)])
    }
    
    # f
    # plot(f$x, f$y, ylim=c(1.5,-1.5))
    # dev.copy(png, 'MyFigure.png')
    # dev.off()
    
    # Assigning x as f$x
    x <- f$x
    # x
    
    # Lenght trainning set
    nn <- length(f$x)
    # nn
    
    # The standard deviation of the noise
    sigma.n <- f$z
    # sigma.n
  
    #------------------------------------------------
    
    #   COMPUTING THE MEAN FUNCTION AND ITS VARIANCE 
    
    # Define the test points at which we want to predict the functions:
    x.star <- seq(round(x[1]), round(x[nn]), by=StepSizeTestPoints)
    # x.star
    # length(x.star)
    # x[1]
    # round(x[1])
    # x[nn]
    # round(x[nn])
  
    numberTestPoints <- length(x.star)
    # numberTestPoints
    
    
    # COVARIANCE MATRICES TO CONSTRUCT EQ. (2.21)
    
    # Peculiar velocity covariance matrix for the training data, K(X, X):
    # OLD. sigma2PecuMatrix_xx <- matrix(sigma2muPecu(velPecuFix, redshiftSN, zerrorFix), nn, nn)
    sigma2PecuMatrix_xx <- matrix(sigma2muPecu(velPecuFix, redshiftSN, zcmb_error), nn, nn)
  
    # sigma2muPecu(velPecuFix, redshiftSN, zerrorFix)
    # sigma2PecuMatrix_xx
    # dim(sigma2PecuMatrix_xx)
    
    # Covariance matrix  K(X, X):
    k.xx <- kernelFunc(x,x,l_Fix, sigma2kern_fix) + (sigma.n^2)*diag(1,nn) + sigma2PecuMatrix_xx
    
    # k.xx_mean matrix to be used to determine the -mean- function based on the data and ignoring the pec vel matrix
    k.xx_mean <- kernelFunc(x,x,l_Fix, sigma2kern_fix) + (sigma.n^2)*diag(1,nn)
    # k.xx
    # dim(k.xx)
    # class(k.xx)
    # kernelFunc(x,x,l_Fix, sigma2kern_fix)
    # kernelFunc(x,x,l_Fix, 1.1)
    # solve(k.xx)
    # sigma.n^2*diag(1,nn)
    # diag(1,nn)
    # sigma.n
    # k.xx_mean
    # kernelFunc(x,x,l_Fix, sigma2kern_fix)
    
    
    # Peculiar velocity covariance matrix for  K(X, X*):
    sigma2PecuMatrix_xxs <- matrix(sigma2muPecu(velPecuFix, redshiftSN, zcmb_error), nn, numberTestPoints)
    # dim(sigma2PecuMatrix_xxs)
    # sigma2PecuMatrix_xxs
    # Covariance matrix  K(X, X*):
    if (UsePecMatrix == TRUE){
      k.xxs <- kernelFunc(x,x.star,l_Fix, sigma2kern_fix) + sigma2PecuMatrix_xxs 
      } else if(UsePecMatrix == FALSE) { 
      k.xxs <- kernelFunc(x,x.star,l_Fix, sigma2kern_fix) }
    # k.xxs
    # dim(k.xxs)
    # class(k.xxs)
    
    # Peculiar velocity covariance matrix for  K(X*, X):
    sigma2PecuMatrix_xsx <- matrix(sigma2muPecu(velPecuFix, redshiftSN, zcmb_error), numberTestPoints, nn)
    # dim(sigma2PecuMatrix_xsx)
    # sigma2PecuMatrix_xsx
    # Covariance matrix  K(X*, X):
    if (UsePecMatrix == TRUE){
      k.xsx <- kernelFunc(x.star,x,l_Fix, sigma2kern_fix) + sigma2PecuMatrix_xsx 
    } else if (UsePecMatrix == FALSE) {
      k.xsx <- kernelFunc(x.star,x,l_Fix, sigma2kern_fix) }

    # k.xsx
    # dim(k.xsx)
    # class(k.xsx)
    
    # Peculiar velocity covariance matrix for  K(X*, X*):
    sigma2PecuMatrix_xsxs <- matrix(sigma2muPecu(velPecuFix, redshiftSN, zcmb_error), numberTestPoints, numberTestPoints)
    # dim(sigma2PecuMatrix_xsxs)
    # sigma2PecuMatrix_xsxs
    # Covariance matrix  K(X*, X*):
    if (UsePecMatrix == TRUE) { 
      k.xsxs <- kernelFunc(x.star,x.star,l_Fix, sigma2kern_fix) + sigma2PecuMatrix_xsxs 
      } else if (UsePecMatrix == FALSE) { 
      k.xsxs <- kernelFunc(x.star,x.star,l_Fix, sigma2kern_fix) }
    # k.xsxs
    # dim(k.xsxs)
    # class(k.xsxs)
    
    #----------------------------
    
    # Calculate the mean and covariance functions accounting for the noise in the data:
    # THE MEAN [EQ. (2.23)].
    # In the following equation I have removed the correlated matrix of the noise in order to obtain a better determination of the mean function, otherwise, I will obtain something shifted upward or downward when the error bars are big (usually when peculiar vel uncertainty is big) and the LC is not centered with respect to the mean template.
    # f.bar.star <- k.xsx %*% solve(k.xx) %*% f$y # pec vel
    # f.bar.star <- k.xsx %*% solve(k.xx_mean) %*% f$y # no pec vel
    
    if (ComputeMeanUsingPecMatrix == FALSE) {
      f.bar.star <- k.xsx %*% solve(k.xx_mean) %*% f$y 
    } else if (ComputeMeanUsingPecMatrix == TRUE) {
      f.bar.star <- k.xsx %*% solve(k.xx) %*% f$y}
    
    # if (MeanTrick == FALSE) {f.bar.star <- k.xsx %*% solve(k.xx) %*% f$y}
    # else if (MeanTrick == TRUE) {f.bar.star <- k.xsx %*% solve(k.xx_mean) %*% f$y}
    # f.bar.star
    # dim(f.bar.star)
    
    # THE COVARIANCE [EQ. (2.24)]:
    cov.f.star <- k.xsxs - (k.xsx %*% solve(k.xx) %*% k.xxs)
    # cov.f.star
    # dim(cov.f.star)
    # k.xsxs
  
    #----------------------------
    
    #  THE STANDARD ERROR OF EACH PREDICTED POINT
    # Extracting the diagonal elements of the covariance matrix = variances of the predicted data
    varianceCovMatrix <- diag(cov.f.star)
    # length(varianceCovMatrix)
    # varianceCovMatrix
    
    # Converting the variance to standard error for the mean predicted data
    StdErrorMean <- sqrt(varianceCovMatrix)
    # length(StdErrorMean)
    # StdErrorMean
    
  
    #===============================================================
    
    # SAVE THE DATA (phase, mean, stdErrorMean)
    
    # Array (phase, mean, stdErrorMean)
    phase_mu_stdError <- 0
    
    # Adding the mean prior ('meanPriorFix' or template) to go back to the actual values of magnitude. 
    
    if (GP_prior == 'flat'){
      MagPlusTemp <-  f.bar.star + meanPriorFix
    } else if (GP_prior == 'template'){
      MagPlusTemp <-  f.bar.star + TempInterpol(x.star)}
    
    # Combing (phase, mean, stdErrorMean) arrays in just one array.
    # phase_mu_stdError <- cbind(x.star, f.bar.star+meanPriorFix, StdErrorMean) # old
    phase_mu_stdError <- cbind(x.star, MagPlusTemp, StdErrorMean)
    
    # Adding a name to each column:
    colnames(phase_mu_stdError) <- c("phase", "mean", "std error")
    # phase_mu_stdError
    
    MyPathAndName <- file.path(DirSaveOutput, paste(NameSN, '_GP_mean_sigma.dat', sep = ''))
    # OLD. MyPathAndName <- file.path(MainDir, paste(KindOfData, '/', NameSN, "_GP_mean_sigma", ".dat", sep = ""))
    
    # Writting the data to a file
    write.table(phase_mu_stdError, file=MyPathAndName, sep="   ", row.names = FALSE, col.names = FALSE)
    
    #----------------------------
    
    #     CREATING AND SAVING THE 'FILLED' FILE TO BE USED IN THE HIERARCHICAL CODE
    
    # The actual predicted values from GP
    # df_GPfit <- data.frame(phase=x.star, mean=f.bar.star+meanPriorFix, stdErr=StdErrorMean) # old
    df_GPfit <- data.frame(phase=x.star, mean=MagPlusTemp, stdErr=StdErrorMean)
    # df_GPfit 
    
    #   The filling values
    
    xFillDown <- seq(phase_lowerLimit, round(x[1])-StepSizeTestPoints, by=StepSizeTestPoints)
    # xFillDown
    yFillDown <- rep(20,length(xFillDown))
    # yFillDown
    yErrorFillDown <- rep(21,length(xFillDown))
    # yErrorFillDown
    dfFill.down <- data.frame(phase = xFillDown, mean = yFillDown, stdErr = yErrorFillDown)
    
    # In case of there are GP predicted values with high phase:
    if(round(x[nn]) < phase_upperLimit ){
    xFillUp <- seq(round(x[nn])+StepSizeTestPoints, phase_upperLimit, by=StepSizeTestPoints)
    # xFillUp
    yFillUp <- rep(20,length(xFillUp))
    # yFillUp
    yErrorFillUp <- rep(21,length(xFillUp))
    # yErrorFillUp
    dfFill.up <- data.frame(phase = xFillUp, mean = yFillUp, stdErr = yErrorFillUp)
    
    # Putting together the filled data.
    df_GPfit_Filled <- rbind(dfFill.down, df_GPfit, dfFill.up)
    } else { # Putting together the filled data.
      df_GPfit_Filled <- rbind(dfFill.down, df_GPfit)  }
    
    # Writting the file
    MyPathAndName2 <- file.path(DirSaveOutput, paste(NameSN, '_GP_mean_sigma_Filled.dat', sep = ''))
    # OLD. MyPathAndName2 <- file.path(MainDir, paste(KindOfData, '/', NameSN, '_GP_mean_sigma_Filled', '.dat', sep = ''))
    write.table(df_GPfit_Filled, file=MyPathAndName2, sep='  ', row.names = FALSE, col.names = FALSE)
    
    #----------------------------
    # COVARIANCE MATRIX OF THE -DATA- (i.e., training data) only.
    "
    CovMatrix_Data <- k.xx
    # dim(CovMatrix_Data)
    
    # Writting the matrix to a file
    # MyPathAndName1 <- file.path(MainDir, paste(KindOfData, '/', NameSN, '_CovMat.dat', sep = ''))
    # write.table(CovMatrix_Data, file=MyPathAndName1, sep='   ', row.names = FALSE, col.names = FALSE)
    
    # Combining data and the their covariance matrix in a single array
    CovMatrixData_Data <- cbind(DataLC$V1[9:LengthData], DataLC$V2[9:LengthData], DataLC$V3[9:LengthData], CovMatrix_Data)
    
    MyPathAndName2 <- file.path(DirSaveOutput, paste(NameSN, '_CovMatData.dat', sep = ''))
    text10 <- '# The first 3 column are the data (phase, Abs Mag, error_abs_mag), and the rest of the columns is the covariance matrix of the data'
    write.table(text10, file=MyPathAndName2, quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(CovMatrixData_Data, file=MyPathAndName2, sep='   ', row.names = FALSE, col.names = FALSE, append=TRUE)
    "
    
    #-------------------
    # SAVING THE OTHER COVARIANCE MATRICES.
    
    "
    # Ok but I don't need this part anymore.

    MyPathAndName2 <- file.path(DirSaveAppmag, paste(NameSN, '_CovMatAppmag.dat', sep = ''))
    text10 <- '# Covariance matrix of the apparent magnitude data. Cov = K(x,x) + sigma^2_i, i.e., without the peculiar velocity uncertainty.'
    write.table(text10, file=MyPathAndName2, quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(k.xx_mean, file=MyPathAndName2, sep='   ', row.names = FALSE, col.names = FALSE, append=TRUE)
  
    
    # Writting the matrix to a file
    MyPathAndName1 <- file.path(MainDir, paste(KindOfData, '/', NameSN, '_CovMat_xxs.dat', sep = ''))
    write.table(k.xxs, file=MyPathAndName1, sep='   ', row.names = FALSE, col.names = FALSE)
    
    # Writting the matrix to a file
    MyPathAndName1 <- file.path(MainDir, paste(KindOfData, '/', NameSN, '_CovMat_xsxs.dat', sep = ''))
    write.table(k.xsxs, file=MyPathAndName1, sep='   ', row.names = FALSE, col.names = FALSE)  
    
    # Writting the matrix to a file
    MyPathAndName1 <- file.path(MainDir, paste(KindOfData, '/', NameSN, '_CovMat_xsx.dat', sep = ''))
    write.table(k.xsx, file=MyPathAndName1, sep='   ', row.names = FALSE, col.names = FALSE)  
    "
  
    #=================================================================
    

    #     PLOTTING 1
    
    # The standard deviation of the noise including pecular velocities.
    # THIS IS USED FOR PLOTTING PURPOSES ONLY
    sigmaData_PecVel.n <- sqrt((sigma.n)^2 + sigma2muPecu(velPecuFix, redshiftSN, zcmb_error))
    # sigmaData_PecVel.n <- sqrt((sigma.n)^2) #TMP
    # sigmaData_PecVel.n
    # length(sigmaData_PecVel.n)
    
    # Putting the predicted values of the mean function in a data frame:
    d <- 0
    
    # d <- data.frame(phase=phase_mu_stdError[c(1:numberTestPoints),1], 
    #                 mean=phase_mu_stdError[c(1:numberTestPoints),2], 
    #                 lower=phase_mu_stdError[c(1:numberTestPoints),2]-phase_mu_stdError[c(1:numberTestPoints),3], 
    #                 upper=phase_mu_stdError[c(1:numberTestPoints),2]+phase_mu_stdError[c(1:numberTestPoints),3])
    
    d <- data.frame(phase=x.star,   mean=(MagPlusTemp), 
                    lower = (MagPlusTemp)-StdErrorMean,
                    upper = (MagPlusTemp)+StdErrorMean)
    
    # Path and name of the figure to be saved
    # MyPathAndNamePlot <- file.path('/Users/arturo/Dropbox/Research/Articulos/10-AndyKaisey/10Compute/NIR_Template/GaussianProcess/TheTemplates/J_band/2_Fit_GP_IndividualSNe_R', paste(KindOfData, '/GP_fit/', NameSN, '_GP_plot', '.png', sep = ''))
    
    MyPathAndNamePlot <- file.path(DirSaveOutput, paste(NameSN, '_GP_plot.png', sep = ''))
    # OLD. MyPathAndNamePlot <- file.path(MainDir, paste(KindOfData, '/', NameSN, '_GP_plot', '.png', sep = ''))
    # MyPathAndNamePlot
    
    
    # Plotting the data, the mean function and its standard error
    gg2 <- ggplot() +
       # PLOT THE GP VARIANCE BAND
       geom_errorbar(data=d, mapping=aes(x=phase, ymin=upper, ymax=lower), width=1.1, size=1, color='gray', alpha=1) + 
       # PLOT THE GP MEAN FUNCTION
       geom_line(data=d, aes(x=phase, y=mean), colour='red', size=0.8) + 
       theme_bw() + # Making the plot with white background
       ggtitle(NameSN) +
       # labs(x='Phase = (MJD - T_Bmax)/(1 + z_hel)', y='Absolute Magnitude') +
       labs(x=expression(Phase = (MJD - T[Bmax])/(1 + z[hel])), y='Absolute Magnitude') +
       # Characteristics of the text of the title and axis labels
       theme(plot.title = element_text(family = 'Trebuchet MS', color='#666666', face='bold', size=10)) +
       theme(axis.title = element_text(family = 'Trebuchet MS', color='#666666', face='bold', size=9)) +
       scale_y_reverse(lim=c(ymin_plot, ymax_plot)) + 
       # scale_y_reverse(lim=c(-6,-18)) +  # FOR LOW LUMINOSITY SNe
       xlim(-10,60) +
       # PLOT THE DATA POINTS AND DATA BARS
       # geom_point(data=f,aes(x=x,y=(y+TempInterpol(x))), size=0.5) +
       if (GP_prior == 'flat') {
       geom_errorbar(data=f,aes(x=x,y=(y+meanPriorFix),ymin=(y+meanPriorFix)-sigmaData_PecVel.n, ymax=(y+meanPriorFix)+sigmaData_PecVel.n), width=0.5, color='black')
       } else if (GP_prior == 'template') { 
         geom_errorbar(data=f,aes(x=x,y=(y+TempInterpol(x)),ymin=(y+TempInterpol(x))-sigmaData_PecVel.n, ymax=(y+TempInterpol(x))+sigmaData_PecVel.n), width=0.5, color='black') }
    
    gg2 + if (GP_prior == 'flat') {
      geom_point(data=f,aes(x=x,y=(y+meanPriorFix)), size=0.5)
    } else if (GP_prior == 'template') { geom_point(data=f,aes(x=x,y=(y+TempInterpol(x))), size=0.5) }
      
    # Adding text to the plot
    gg2 + annotate('text', x = xText, y = yText, label = c('d15 = ', round(delta15,3), '+-', round(delta15_error,3), 'z_CMB =', round(redshiftSN, 4), 'l =', round(l_Fix,3), 'sigma_kern =', round(sqrt(sigma2kern_fix),3), 'distance mu=', round(DistanceMu,3), '+-', round(DistanceMu_error,3), '68.3% confidence interval' ), size = 2, hjust = 0)
    
    # ggsave(MyPathAndNamePlot, width = 12, height = 9, units = 'cm')
    ggsave(MyPathAndNamePlot, width = 10, height = 7.5, units = 'cm')
    # ggsave(MyPathAndNamePlot, width = 8, height = 6, units = 'cm'')
    # ggsave(MyPathAndNamePlot)
    
    #----------------------------
    
    #     PLOTTING 2
    # This is basically a copy of the section PLOTTING 1
    # In this section I plot the actual residual data (data minus moving window template) that
    # it is what I'm actually fitting with Gaussian Processes.
    
    # Putting the predicted values of the mean function in a data frame:
    d <- 0
    
    # d <- data.frame(phase=phase_mu_stdError[c(1:numberTestPoints),1], 
    #                 mean=phase_mu_stdError[c(1:numberTestPoints),2], 
    #                 lower=phase_mu_stdError[c(1:numberTestPoints),2]-phase_mu_stdError[c(1:numberTestPoints),3], 
    #                 upper=phase_mu_stdError[c(1:numberTestPoints),2]+phase_mu_stdError[c(1:numberTestPoints),3])
    
    d <- data.frame(phase=x.star, mean=(f.bar.star), 
                    lower = (f.bar.star)-StdErrorMean,
                    upper = (f.bar.star)+StdErrorMean)
    
    # Path and name of the figure to be saved
    # MyPathAndNamePlot <- file.path('/Users/arturo/Dropbox/Research/Articulos/10-AndyKaisey/10Compute/NIR_Template/GaussianProcess/TheTemplates/J_band/2_Fit_GP_IndividualSNe_R', paste(KindOfData, '/GP_fit/', NameSN, '_GP_plot', '.png', sep = ''))
    
    MyPathAndNamePlot <- file.path(DirSaveOutput, paste(NameSN, '_GPResidual.png', sep = ''))
    # OLD. MyPathAndNamePlot <- file.path(MainDir, paste(KindOfData, '/', NameSN, '_GPResidual', '.png', sep = ''))
    # MyPathAndNamePlot
    
    # Plotting the data, the mean function and its standard error
    gg3 <- ggplot() +
       # PLOT THE GP VARIANCE BAND
       geom_errorbar(data=d, mapping=aes(x=phase, ymin=upper, ymax=lower), width=1.1, size=1, color='gray', alpha=1) + 
       # PLOT THE GP MEAN FUNCTION
       geom_line(data=d, aes(x=phase, y=mean), colour='red', size=0.8) + 
       theme_bw() + # Making the plot with white background
       ggtitle(NameSN) +
       # labs(x='Phase = (MJD - T_Bmax)/(1 + z_hel)', y='Absolute Magnitude') +
       labs(x=expression(Phase = (MJD - T[Bmax])/(1 + z[hel])), y='Absolute Magnitude') +
       # Characteristics of the text of the title and axis labels
       theme(plot.title = element_text(family = 'Trebuchet MS', color='#666666', face='bold', size=10)) +
       theme(axis.title = element_text(family = 'Trebuchet MS', color='#666666', face='bold', size=9)) +
       scale_y_reverse(lim=c(1.5, -1.5)) + 
       # ylim(-1.5,1.5) +
       xlim(-10,60) +
       # PLOT THE DATA POINTS AND DATA BARS
       geom_point(data=f,aes(x=x,y=y), size=0.5) +
       geom_errorbar(data=f,aes(x=x,y=NULL,ymin=y-sigmaData_PecVel.n, ymax=y+sigmaData_PecVel.n), width=0.5, color='black')
  
    # Adding text to the plot
    # gg3 + annotate('text', x = xText, y = yText, label = c('d15 = ', round(delta15,3), '+-', round(delta15_error,3), 'z_CMB =', round(redshiftSN, 4), 'l =', round(ResultOptim$par[1],3), 'sigma_kern =', round(sqrt(ResultOptim$par[2]),3), 'distance mu=', round(DistanceMu,3), '+-', round(DistanceMu_error,3), '68.3% confidence interval' ), size = 2, hjust = 0)
    
    # ggsave(MyPathAndNamePlot, width = 12, height = 9, units = 'cm')
    ggsave(MyPathAndNamePlot, width = 10, height = 7.5, units = 'cm')
    # ggsave(MyPathAndNamePlot, width = 8, height = 6, units = 'cm'')
    # ggsave(MyPathAndNamePlot)
    
    #---------------------------- 
    
    #     PLOTTING 3
    # Plot the f data frame. The dots in this plot should be the same to the of PLOTTING 2.
    
    "
    MyPathAndNamePlot <- file.path(MainDir, paste(KindOfData, '/', NameSN, '_GPResidual2_plot', '.png', sep = ''))
    plot(f$x, f$y, ylim=c(1, -1))
    grid()
    dev.copy(png, MyPathAndNamePlot)
    dev.off()  
    "
    
    cat(NameSN)
    cat(' -  ')
  }
}

cat('--- All', i, 'SNe fitted smoothly ---')

################### end main loop #########################

# Copy/pasting and removing the photometric LC files that pass the quality cutoffs to the the 'Goods' folder. 

DirSaveOutputGood <- file.path(DirLCData,'Goods/')
dir.create(DirSaveOutputGood)

for(i in 1:numSNeAfterCutoff){
  
  # Defined the name of the SN by the name of its data file
  # Removing the extension ".txt" to the filename.
  NameSN_copy <- substrRight(paste('', list_SNeAfterCutoff$V1[i],'', sep = ''), 4)
  
  #      Copy the files
  # file.copy(paste(NameSN_copy, '_CovMatData.dat', sep = ''), to=DirSaveOutputGood)
  file.copy(paste(NameSN_copy, '_GP_mean_sigma_Filled.dat', sep = ''), to=DirSaveOutputGood)
  file.copy(paste(NameSN_copy, '_GP_mean_sigma.dat', sep = ''), to=DirSaveOutputGood)
  file.copy(paste(NameSN_copy, '_GP_plot.png', sep = ''), to=DirSaveOutputGood)
  file.copy(paste(NameSN_copy, '_GPResidual.png', sep = ''), to=DirSaveOutputGood)
  file.copy(paste(NameSN_copy, '.txt', sep = ''), to=DirSaveOutputGood)
  
  #      Now remove the originals to avoid to have duplicates
  # file.remove(paste(NameSN_copy, '_CovMatData.dat', sep = ''))
  file.remove(paste(NameSN_copy, '_GP_mean_sigma_Filled.dat', sep = ''))
  file.remove(paste(NameSN_copy, '_GP_mean_sigma.dat', sep = ''))
  file.remove(paste(NameSN_copy, '_GP_plot.png', sep = ''))
  file.remove(paste(NameSN_copy, '_GPResidual.png', sep = ''))
  file.remove(paste(NameSN_copy, '.txt', sep = ''))
}

cat('--- All done:', i, 'fitted and good (passed the cutoffs) SNe were copied/pasted/removed to the good folder---')

