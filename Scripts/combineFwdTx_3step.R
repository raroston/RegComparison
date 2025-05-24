library(ANTsR)

setwd("/home2/rachel/P01/Github/RegComparison/")

experiment = "Exp04"

out.dir = paste0("./Results/", experiment, "/Transforms/")
img.dir = "/home2/rachel/P01/MagaLab/Data/RigidAligned/lowRes_Volumes/"

ref.img = antsImageRead("/home2/rachel/P01/MagaLab/Data/RigidAligned/maskedtemplate0__lowRes.nii.gz")

scanIDs = unlist(lapply(strsplit(dir(out.dir, pattern = "inv"), "_Exp"), "[", 1))

for(i in 1:length(scanIDs)){
  mov.img = antsImageRead(dir(img.dir, pattern = scanIDs[i], full.names = T))
  Warp2 = paste0(out.dir, scanIDs[i], "_", experiment, "_Step3_syn2_2Warp.nii.gz")
  Warp1 = paste0(out.dir, scanIDs[i], "_", experiment, "_Step3_syn2_1Warp.nii.gz")
  GenericAffine = paste0(out.dir, scanIDs[i], "_", experiment, "_Step3_syn2_0GenericAffine.mat")
  
  antsApplyTransforms(fixed = ref.img, 
                      moving = mov.img, 
                      transformlist = c(Warp2, Warp1, GenericAffine),
                      interpolator = "linear", 
                      compose = paste0(out.dir, scanIDs[i], "_", experiment, "_allFwdTx-"))
  remove(mov.img)
  gc()
}