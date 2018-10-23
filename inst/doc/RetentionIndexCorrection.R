#########################################################################

## RI correction when standards are not co-injected with biological Samples.
## ----------------------------------------------------------------------
##
## TargetSearch assumes that the retention index markers (RIM) are injected
## together with the biological samples. A different approach is to inject
## the n-alkanes or FAMEs separately, between the biological samples.
## For example, your GC-MS run may look like this

## Measurement Order  Sample Type
## -----------------  -----------
##            1       Alkanes
##            2       Biological
##            3       Biological
##            4       Biological
##            5       Biological
##            6       Alkanes
##            7       Biological
##            8       Biological
##            9       Biological
##           10       Biological
##           11       Alkanes
##           12       Biological
##           13       Biological
##           14       Biological
##           15       Biological

## In the example, samples 1, 6, 11 are RI markers and the rest are biological
## samples. The assumption is that the retention time shifts between consecutive
## runs are not significant, so sample #1 is used to correct samples #2-5,
## sample #6 corrects samples #7-10, and so on. This script does exactly that
## with the available chromatograms of package TargetSearchData.

#########################################################################

library(TargetSearch)
library(TargetSearchData)

# get the directory where TargetSearchData's CDF files are located
cdfPath <- file.path(find.package("TargetSearchData"), "gc-ms-data")
# show the cdf files
dir(cdfPath, pattern="cdf$")

# create a sample object using all chromatograms 
( samples.all <- ImportSamplesFromDir(cdfPath) )

# set the working directory where the RI files will be saved.
# ("." is the current directory. Run "getwd()" to find out which it is.) 
RIpath(samples.all) <- "."

#import retention index marker times limits
rimLimits <- ImportFameSettings(file.path(cdfPath,"rimLimits.txt"))
rimLimits

# run Retention index correction (Intensity Threshold = 50;
# peak detection method = "ppc", window = 15)
RImatrix <- RIcorrect(samples.all, rimLimits,
            Window = 15, pp.method = "ppc", IntThreshold = 50)

# Create a logical vector indicating the RI standards:
# set TRUE to the samples that contain the RI standards (samples 1,6,11 in this example)
# set FALSE otherwise. Note that there are 15 chromatograms.
isRIMarker <- c(T,F,F,F,F,T,F,F,F,F,T,F,F,F,F)

# make a copy of RImatrix
RImatrix2 <- RImatrix

# Copy the retention times of the standard with the corresponding biological samples,
RImatrix2[, 2:5]    <- RImatrix[,1]
RImatrix2[, 7:10]   <- RImatrix[,6]
RImatrix2[, 12:15]  <- RImatrix[,11]

# update the RIs of the biological samples.
fixRI(samples.all, rimLimits, RImatrix2, which(!isRIMarker))

# remove the standards, since we don't need them anymore.
samples <- samples.all[!isRIMarker]
RImatrix <- RImatrix2[, !isRIMarker]

# If you prefere to use the GUI, run TargetSearchGUI() from this point
# and import the RI files by selecting the option "Apex Data" (don't
# import the standard files)

# continue with the normal work flow (see TargetSearch vignette)
