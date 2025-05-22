# Registration Exp03
# Rachel Roston

# Exp01 = antsRegistration only
# Exp02 = antsRegistration affine -> labelImageRegistration SyN only -> antsRegistration SyN only
# Exp03 = labelImageRegistration affine & SyN -> antsRegistration antsRegistration SyN only
# Exp04 = labelImageRegistration similarity & SyN -> antsRegistration SyN only

experiment = "Exp03"

wd = "/home2/rachel/P01/Github/RegComparison/"
setwd(wd)

library(ANTsR)

  exp.dir = paste0("./Results/", experiment)
  if(!dir.exists(exp.dir)) dir.create(exp.dir)
  
  out.dir = paste0("./Results/", experiment, "/Transforms")
  if(!dir.exists(out.dir)) dir.create(out.dir)

  img.dir = "/home2/rachel/P01/MagaLab/Data/RigidAligned/lowRes_Volumes/"
  label.dir = "/home2/rachel/P01/MagaLab/Data/RigidAligned/lowRes_15Labels/"
  mask.dir = "/home2/rachel/P01/MagaLab/Data/RigidAligned/lowRes_Masks-wholebody/"
  
  ref.image = antsImageRead("/home2/rachel/P01/MagaLab/Data/RigidAligned/maskedtemplate0__lowRes.nii.gz")
  ref.label = antsImageRead("/home2/rachel/P01/MagaLab/Data/RigidAligned/maskedtemplate0__lowRes-15-label.nii.gz")
  ref.mask = antsImageRead("/home2/rachel/P01/MagaLab/Data/RigidAligned/maskedtemplate0__lowRes-wholebody-label.nii.gz")
  
  subjects = c("Scan_0082",
               "Scan_0105",
               "Scan_0110",
               "Scan_0125",
               "Scan_0162",
               "Scan_0166",
               "Scan_0176",
               "Scan_0011",
               "Scan_0035",
               "Scan_0038",
               "Scan_0044",
               "Scan_0086",
               "Scan_0093",
               "Scan_0095")
  
  for(i in 1:length(subjects)){
    ScanID = subjects[i]
    out.prefix = paste0(out.dir, "/", ScanID, "_", experiment, "_")
    
    mov.image = antsImageRead(paste0(img.dir, ScanID, "__rec-18um-tx-low.nrrd"))
    mov.label = antsImageRead(paste0(label.dir, ScanID, "__tx-low-15-label.nii.gz"))
    mov.mask = antsImageRead(paste0(mask.dir, ScanID, "__tx-low-wholebody-label.nii.gz"))
    
    # affine.reg = antsRegistration(fixed = ref.image, 
    #                               moving = mov.image,  
    #                               typeofTransform = "Affine")
    # 
    label.reg = labelImageRegistration(fixedLabelImages = ref.label, 
                                       movingLabelImages = mov.label, 
                                       fixedIntensityImages = ref.image, 
                                       movingIntensityImages = mov.image, 
                                       fixedMask = ref.mask, 
                                       movingMask = mov.mask, 
                                       initialTransforms = "affine",
                                       typeOfDeformableTransform = "antsRegistrationSyN[so]", 
                                       outputPrefix = paste0(out.prefix, "Step2_labelReg_"))
    
    syn2 = antsRegistration(fixed=ref.image, 
                            moving=mov.image, 
                            mask = ref.mask, 
                            movingMask = mov.mask,
                            typeofTransform = "antsRegistrationSyN[so]", 
                            initialTransform = label.reg$fwdtransforms, 
                            outprefix = paste0(out.prefix, "Step3_syn2_"))
    
    antsApplyTransforms(fixed= ref.image, 
                        moving = mov.image, 
                        transformlist = syn2$fwdtransforms[1:2], 
                        compose = paste0(out.prefix, "combinedWarps-"))
    
    # Linear transformed image
    linear.txImg = antsApplyTransforms(fixed = ref.image, 
                                       moving = mov.image, 
                                       transformlist = syn2$fwdtransforms[3], 
                                       interpolator = "linear")
    antsImageWrite(linear.txImg, paste0(out.prefix, "Step1_linearTx_fwdWarpedImage.nrrd"))
    
    # Linear+LabelRegistration transformed image
    label.reg.txImg = antsApplyTransforms(fixed = ref.image, 
                                          moving = mov.image, 
                                          transformlist = label.reg$fwdtransforms, 
                                          interpolator = "linear")
    antsImageWrite(label.reg.txImg, paste0(out.prefix, "Step2_Warp1Tx_fwdWarpedImage.nrrd"))
    
    # Linear+Label+syn2 transformed image
    syn2.txImg = antsApplyTransforms(fixed = ref.image, 
                                     moving = mov.image, 
                                     transformlist = syn2$fwdtransforms, 
                                     interpolator = "linear")
    antsImageWrite(syn2.txImg, paste0(out.prefix, "Step3_Warp2Tx_fwdWarpedImage.nrrd"))
    
    # Inverse Labels
    inverse.labels = antsApplyTransforms(fixed = mov.image, 
                                         moving = ref.label,
                                         transformlist = c(label.reg$invtransforms[1],
                                                           label.reg$invtransforms[2],
                                                           syn2$invtransforms[3]),
                                         whichtoinvert = c(T,F,F), 
                                         interpolator = "genericLabel")
    antsImageWrite(inverse.labels, paste0(out.prefix, "inverse-label.nii.gz"))
    
    remove(mov.image, mov.label, mov.mask, 
           affine.reg, label.reg, syn2, 
           linear.txImg, label.reg.txImg, syn2.txImg)
    gc()
  }