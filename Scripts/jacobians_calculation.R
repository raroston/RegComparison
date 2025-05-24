# Calculate Jacobians
# Rachel Roston, Ph.D.
library(ANTsR)

setwd("/home2/rachel/P01/Github/RegComparison/")

source("./Scripts/saveJacobian.R")

experiment = "Exp01"

dir.input = paste0("./Results/", experiment, "/Transforms/")
dir.output = paste0("./Results/", experiment, "/Jacs-1Warp")
transformPattern = "1Warp"


if(!dir.exists(dir.output)) dir.create(dir.output)

ref.img = antsImageRead("/home2/rachel/P01/MagaLab/Data/RigidAligned/maskedtemplate0__lowRes.nii.gz")

transforms = dir(dir.input, pattern = transformPattern)
jacs.filenames = gsub(".nii.gz", "-Jac.nrrd", transforms)

for(i in 1:length(transforms)){
  saveJacobian(transforms[i], 
               dir.transforms = dir.input, 
               dir.out = dir.output,
               ref = ref.img)
  print(paste("Completed", i, "of", length(transforms)))
}
