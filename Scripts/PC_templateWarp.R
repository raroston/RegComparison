setwd("/home2/rachel/P01/Github/RegComparison/")

library(ANTsR)

experiment = "Exp04"

pca.dir = paste0("./Results/", experiment, "/PCA-allFwdTx")

load(paste0(pca.dir, "/PCA"))

ref.img.path = "/home2/rachel/P01/MagaLab/Data/RigidAligned/maskedtemplate0__lowRes.nii.gz"
ref.mask.path = "/home2/rachel/P01/MagaLab/Data/RigidAligned/maskedtemplate0__lowRes-wholebody-label.nii.gz"
ref.label.path = "/home2/rachel/P01/MagaLab/Data/RigidAligned/maskedtemplate0__lowRes-15-label.nii.gz"
ref.img = antsImageRead(ref.img.path)
ref.mask = antsImageRead(ref.mask.path)
ref.label = antsImageRead(ref.label.path)

pca_mask = iMath(ref.mask, "MD", 15)

PCs = 3:4
PC_warpedImg.dir = paste0(pca.dir, "/plots/PC_warpedTemplate/")
if(!dir.exists(PC_warpedImg.dir))dir.create(PC_warpedImg.dir)
scales = 100*(c(-8:-1,1:8))

for(i in PCs){
  warped.img = list()
  warped.label = list()
  for(j in 1:length(scales)){
    PC_scale = scales[j]
    upscaled = vectorToMultichannel(PC_scale*pca$pca$v[,i], pca_mask)   #this should be the mask you used for your PCA
    upscaledTX = antsrTransformFromDisplacementField(upscaled)
    warped.img[[j]] = applyAntsrTransform(upscaledTX, 
                                          data=ref.img,
                                          reference= ref.img, 
                                          interpolation = "linear")
    warped.label[[j]] = applyAntsrTransform(upscaledTX, 
                                            data = ref.label,
                                            reference= ref.label, 
                                            interpolation = "nearestNeighbor")
    antsImageWrite(warped.img[[j]], paste0(PC_warpedImg.dir, "PC", i, "_x", PC_scale, "_", experiment, ".nrrd"))
    antsImageWrite(warped.label[[j]], paste0(PC_warpedImg.dir, "PC", i, "_x", PC_scale, "_", experiment, "-label.nii.gz"))
  }
  
}

rm(list = ls())
gc()
