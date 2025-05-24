
library(ANTsR)
setwd("/home2/rachel/P01/Github/RegComparison")
source("./Scripts/jacobianStats.R")

experiment = "Exp02"
transformPattern = "allFwdTx"
jacs.dir = paste0("./Results/", experiment, "/Jacs-allFwdTx")
pAdj = "fdr"


results.dir = paste0(jacs.dir, "/TBM/")
if(!dir.exists(results.dir)) dir.create(results.dir)
filenames = paste0(experiment, "_", transformPattern, "_", pAdj, "_")

# Input files
specimen.database = read.csv("./Data/scan_metadata.csv")
ref.img.path = "/home2/rachel/P01/MagaLab/Data/RigidAligned/maskedtemplate0__lowRes.nii.gz"
ref.mask.path = "/home2/rachel/P01/MagaLab/Data/RigidAligned/maskedtemplate0__lowRes-wholebody-label.nii.gz"
ref.img = antsImageRead(ref.img.path)
ref.mask = antsImageRead(ref.mask.path)


# Wrangle data
jacs = lapply(dir(jacs.dir, pattern = "_Jac.nrrd", full.names = T), antsImageRead)
subjects = lapply(strsplit(dir(jacs.dir, pattern = "_Jac.nrrd"), "_"), "[", 1:2)
subjects = unlist(lapply(subjects, function(X){paste(X, collapse = "_")}))
genotypes = vector()
for(i in 1:length(subjects)){
  genotypes[i] = specimen.database[which(specimen.database$ScanID == subjects[i]), "Genotype"]
}
genotypes = relevel(as.factor(genotypes), ref = "WT/WT")

jacobians = jacs
groups = genotypes
binaryLabel = ref.mask
P = 0.05
PadjustMethod = pAdj
outputFolder = results.dir
filePrefix = filenames
save.pImg = T
save.betaImg = T


#set the regression
print("Starting statistical analysis")
j_mat = imageListToMatrix( jacobians, binaryLabel )
j_mdl = lm( j_mat ~ groups)    # raw volume data
j_bmdl = bigLMStats( j_mdl , 1.e-5 )
print( min( p.adjust( j_bmdl$beta.pval, 'none' ) , na.rm=T ) )

#do multiple corrections
#P defines the significance level we want to evaluate things after FDR correction
volstats = j_bmdl$beta.pval

correct.ps= p.adjust(volstats, PadjustMethod)
pImg = makeImage(binaryLabel, correct.ps)

pMask = thresholdImage(pImg, 10^-16, P, 1, 0)
betaImg = makeImage( binaryLabel, j_bmdl$beta ) * pMask

print("Writing images")
if(save.pImg == TRUE){antsImageWrite(pImg, paste0(outputFolder, filePrefix, "_pImg.nii.gz"))}
if(save.betaImg == TRUE){antsImageWrite(betaImg, paste0(outputFolder, filePrefix, "_betaImg.nii.gz"))}

              