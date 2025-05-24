# Save Jacobian determinants
# Author: R.A. Roston, Ph.D., A.M. Maga, Ph.D.

saveJacobian <- function(transform,
                         dir.transforms,
                         dir.out,
                         ref){
  
  library(ANTsR)
  
  # read transform
  warp = antsImageRead(paste0(dir.transforms, transform))
  
  # get jac filename from transform
  jacs.filename = gsub(".nii.gz", "_Jac.nrrd", transform)
  
  # calculate & save jacobian determinant
  if(identical(antsGetSpacing(ref), antsGetSpacing(warp)) == FALSE){
    print("ref image spacing doesn't match transform")
  } else {
    tmp = createJacobianDeterminantImage(domainImg = ref, 
                                         tx = warp, 
                                         doLog = TRUE, 
                                         geom = TRUE )
    jacs.tmp = smoothImage(inimg = tmp, 
                           2*antsGetSpacing(ref)[1], 
                           FALSE ) #smooth the edges.
    
    antsImageWrite(image = jacs.tmp, 
                   filename = paste0(dir.out, "/", jacs.filename))
    
    remove(warp, jacs, jacs.tmp)
    gc()
  }
}