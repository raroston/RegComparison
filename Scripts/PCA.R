
library(ANTsR)
library(stringr)
library(viridis)

setwd("/home2/rachel/P01/Github/RegComparison/")

# Input Variables
experiment = "Exp04"
transformPattern = "combinedWarps"
resultSubdir = "PCA-combinedWarps"
makePlots = TRUE

# Output directory
results.dir = paste0("./Results/", experiment, "/", resultSubdir, "/")
if(! dir.exists(results.dir)) dir.create(results.dir)

# Input: transforms
transform.dir1 = paste0("./Results/", experiment, "/Transforms")
txName1 = transformPattern

# Input files
specimen.database = read.csv("./Data/scan_metadata.csv")
ref.img.path = "/home2/rachel/P01/MagaLab/Data/RigidAligned/maskedtemplate0__lowRes.nii.gz"
ref.mask.path = "/home2/rachel/P01/MagaLab/Data/RigidAligned/maskedtemplate0__lowRes-wholebody-label.nii.gz"
ref.img = antsImageRead(ref.img.path)
ref.mask = antsImageRead(ref.mask.path)

# Organize metadata
transformpaths = dir(transform.dir1, pattern = txName1, full.names = T)
tmp = unlist(strsplit(transformpaths, "/"))
tmp = tmp[grep("Scan", tmp)]
subjects = paste0("Scan_", unlist(lapply(strsplit(tmp, "_"), "[", 2)))
genotypes = vector()
for(i in 1:length(subjects)){
  genotypes[i] = specimen.database[which(specimen.database$ScanID == subjects[i]), "Genotype"]
}
metadata = data.frame(ScanID = subjects, Genotype = genotypes, Registration = experiment)

plotFactor = as.factor(metadata$Genotype)
plotColors = c("black", "blue")
plotShapes =c(16,17)
  
pca_mask = iMath(ref.mask, "MD", 15)

transforms = lapply(transformpaths, antsImageRead)

pca <- multichannelPCA(x = transforms, 
                       mask = pca_mask, 
                       pcaOption = "randPCA")

#Save PCA output
## PCA
out.dir = results.dir
save(objects = pca, file = paste0(out.dir, "/PCA"))
## PCA warps (They do not automatically save with save() )
warp.out.dir = paste0(out.dir, "/pcaWarps")
if(!dir.exists(warp.out.dir)) dir.create(warp.out.dir)
for(j in 1:length(pca$pcaWarps)){
  antsImageWrite(image = pca$pcaWarps[[j]], filename = paste0(warp.out.dir, "/pcaWarp_", j, ".nrrd"))
}
## transformfiles included in the analysis
write.csv(data.frame(transformfiles = transformpaths),
          paste0(out.dir, "/transformlist.csv"))

# PLOTS
if(makePlots == TRUE
   & !is.null(plotFactor)
   & !is.null(plotColors)
   & !is.null(plotShapes)){
  plot.dir = paste0(out.dir, "/plots")
  if(! dir.exists(plot.dir)) dir.create(plot.dir)
  
  ## Scree plot
  jpeg(filename = paste0(plot.dir, "/scree.jpeg"), 
       quality = 100, 
       res = 300,
       units = "in",
       height = 5,
       width = 7)
  eigen.sum=sum(pca$pca$d)
  PC.percent = pca$pca$d/eigen.sum*100
  plot(x = 1:length(pca$pca$d),
       y = PC.percent,
       xlab = "PC",
       ylab = "% Variance Explained",
       main = paste0(experiment, ", ", transformPattern))
  dev.off()
  
  ## PC plots
  if(length(transformpaths) < 17){
    maxPC = floor((length(transformpaths)-1)/2)
  } else if(length(transformfiles) >= 17 ){
    maxPC = 8
  }
  
  for(p in 1:maxPC){
    PCs = c(2*p-1, 2*p)
    jpeg(filename = paste0(plot.dir, "/PC_plot", PCs[1], "-", PCs[2], ".jpeg"), 
         quality = 100,
         res = 300, 
         units = "in",
         height = 5,
         width = 7)
    par(mar= c(5,5,5,8), xpd = TRUE)
    
    #wt = pca$pca$u[which(myfactor == "WT/WT"),PCs]
    #vertices = chull(wt)
    plot(x = pca$pca$u[,PCs[1]], 
         y = pca$pca$u[,PCs[2]], 
         xlab = paste0("PC", PCs[1], " (", signif(PC.percent, digits = 3)[PCs[1]], "%)"), 
         ylab = paste0("PC", PCs[2], " (", signif(PC.percent, digits = 3)[PCs[2]], "%)"), 
         col = plotColors[plotFactor],
         pch = plotShapes[plotFactor],
         asp = 1, 
         main = paste0(experiment, ", ", transformPattern))
    #text(x = pca$pca$u[,PCs[1]], 
    #     y = pca$pca$u[,PCs[2]], 
    #     labels = paste(mystrain, genotypes))
    #polygon(wt[vertices,])
    legend('topright',
           inset = c(-0.22,0),
           legend = levels(plotFactor),
           col = plotColors,
           pch = plotShapes,
           xpd = TRUE)
    dev.off()
  }
  
  # Heatmaps
  for(h in 1:(maxPC*2)){
    tiff(filename = paste0(plot.dir, "/heatmap_PC", h, ".tiff"), 
         res = 300, 
         units = "in",
         height = 5,
         width = 7)
    PC.warp = Reduce('+', x = pca$pcaWarps[h])
    temp =  splitChannels(PC.warp)
    PC.warp.magnitude = iMath( abs( temp[[1]] ) + abs( temp[[2]]) + abs( temp[[3]]  ), "Normalize" )
    
    txt <- list(x = 4, 
                y = 12,
                label = paste("PC", h),
                cex = 1,
                col = 'white')
    plot(x = ref.img, 
         y = PC.warp.magnitude, 
         color.overlay = "viridis",
         axis = 3,
         alpha = 0.7, 
         nslices = 10, 
         ncolumns=10, 
         text = txt, 
         colorbar=TRUE, 
         window.overlay = c(0,1), 
         useAbsoluteScale=TRUE,
         title.img = paste0(experiment, ", ", transformPattern), 
         title.line = -8)
    dev.off()
  }
  
  # PC1 Deformations
  PC1_warpedImg.dir = paste0(plot.dir, "/PC1_warpedTemplate/")
  if(!dir.exists(PC1_warpedImg.dir))dir.create(PC1_warpedImg.dir)
  
  scales = 100*(1:8)
  warped = list()
  for(j in 1:length(scales)){
    PC1_scale = scales[j]
    upscaled = vectorToMultichannel(PC1_scale*pca$pca$v[,1], pca_mask)   #this should be the mask you used for your PCA
    upscaledTX = antsrTransformFromDisplacementField(upscaled)
    warped[[j]] = applyAntsrTransform(upscaledTX, data=ref.img,reference= ref.img)
    antsImageWrite(warped[[j]], paste0(PC1_warpedImg.dir, "PC1_x", PC1_scale, "_", experiment, ".nrrd"))
  }
}

rm(list = ls())
gc()