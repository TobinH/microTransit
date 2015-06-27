# scrap of script for stress-strain analytical modeling of collagen fiber orientation

# code plan:
# 1. per-pixel stress, from cross-sectional properties and mass, matrix notation
#   1a. bending, using "best fit" for neutral axis
#   1b. torsion
# 2. rotate stress matrices opposite CFO 
# 3. calculate per-pixel strain matrix for each state w/:
#    2a. empirical isotropic compliance matrix
#    2b. empirical linear orthotropic compliance matrix
# 4. rotate strain matrices back to global reference frame
# 5. grayscale representation:
#    a. predicted bending strain intensity (iso, ortho, diff)
#    b. predicted torsion strain intensity (iso, ortho, diff)

elastic.model<-function(r.qPLMtab,
                        mass,
                        J){
  norm<-(1/20)
  pois<-(-0.4/20)
  shear<-(1.4/20)
  compliance.iso<-matrix(c(norm, pois, pois,0,0,0,
                         pois,norm,pois,0,0,0,
                         pois,pois,norm,0,0,0,
                         0,0,0,shear,0,0,
                         0,0,0,0,shear,0,
                         0,0,0,0,0,shear), nrow=6)
  # isotropic compliance matrix, estimates of E and v from Carter
  stiffness.lo<-matrix(c(7.02, 3.30, 5.87,0,0,0,
                         3.30,7.02,5.87,0,0,0,
                         5.87,5.87,22.17,0,0,0,
                         0,0,0,6.6,0,0,
                         0,0,0,0,6.6,0,
                         0,0,0,0,0,3.72), nrow=6)
  # linear orthotropic stiffness matrix from Speisz and Zysset 2015
  compliance.lo<-solve(stiffness.lo)
  # linear orthotropic compliance matrix
  sig.zx<-mass*sin(atan2(r.qPLMtab[,7],r.qPLMtab[,6]))*sqrt((r.qPLMtab[,6]^2)+(r.qPLMtab[,7]^2))/J
  sig.zy<-mass*cos(atan2(r.qPLMtab[,7],r.qPLMtab[,6]))*sqrt((r.qPLMtab[,6]^2)+(r.qPLMtab[,7]^2))/J
  # zx and zy torsional stress components modeled proportional to mass
  result<-matrix(0, nrow=nrow(r.qPLMtab), ncol=8)
  result[,1:2]<-r.qPLMtab[,6:7]
  
  for (i in 1:nrow(r.qPLMtab)){
    stress<-matrix(c(0, 0, -sig.zx[i],
                     0, 0, -sig.zy[i],
                     sig.zx[i],sig.zy[i],0), nrow=3)
    # estimated stress in matrix notation
    voigt.stress.iso<-matrix(c(stress[1,1],stress[2,2],stress[3,3],stress[2,3],stress[1,3],stress[1,2]), nrow=6)
    # voigt notation of stress
    Rz<-matrix(c(-cos(r.qPLMtab[i,2]),sin(r.qPLMtab[i,2]),0,
                 -sin(r.qPLMtab[i,2]),-cos(r.qPLMtab[i,2]),0,
                 0,0,1), nrow=3)
    # rotation matrix about z
    Rx<-matrix(c(1,0,0,
                 0,cos(r.qPLMtab[i,1]),sin(r.qPLMtab[i,1]),
                 0,-sin(r.qPLMtab[i,1]),cos(r.qPLMtab[i,1])), nrow=3)
    # rotation matrix about x
    stress.lo<-Rz%*%stress%*%t(Rz)
    # stress rotated in z
    stress.lo<-Rx%*%stress.lo%*%t(Rx)
    # stress rotated in x
    voigt.stress.lo<-matrix(c(stress.lo[1,1],stress.lo[2,2],stress.lo[3,3],stress.lo[2,3],stress.lo[1,3],stress.lo[1,2]), nrow=6)
    # voigt notation of linear isotropic stress
    voigt.strain.iso<-compliance.iso%*%voigt.stress.iso
    # strain given isotropic matrix
    voigt.strain.lo<-compliance.lo%*%voigt.stress.lo
    # strain given linear orthotropic matrix
    strain.iso<-matrix(c(voigt.strain.iso[1,1], voigt.strain.iso[6,1], voigt.strain.iso[5,1],
                       voigt.strain.iso[6,1], voigt.strain.iso[2,1], voigt.strain.iso[4,1],
                       voigt.strain.iso[5,1], voigt.strain.iso[4,1], voigt.strain.iso[3,1]), nrow=3)
    # isotropic material strain into 3x3 matrix form
    strain.lo<-matrix(c(voigt.strain.lo[1,1], voigt.strain.lo[6,1], voigt.strain.lo[5,1],
                        voigt.strain.lo[6,1], voigt.strain.lo[2,1], voigt.strain.lo[4,1],
                        voigt.strain.lo[5,1], voigt.strain.lo[4,1], voigt.strain.lo[3,1]), nrow=3)
    # linear orthotropic material strain into 3x3 matrix form
    strain.lo<-t(Rx)%*%strain.lo%*%Rx
    # rotation about x back to global reference frame
    strain.lo<-t(Rz)%*%strain.lo%*%Rz
    # rotation about z back to global reference frame
    torsion.mag.iso<-sqrt((strain.iso[3,1]^2)+(strain.iso[3,2]^2))
    # magnitude of torsional strain with isotropic material
    torsion.mag.lo<-sqrt((strain.lo[3,1]^2)+(strain.lo[3,2]^2))
    # magnitude of torsional strain with linear orthotropic material
    torsion.diff<-torsion.mag.iso-torsion.mag.lo
    result[i,3:5]<-c(torsion.mag.iso, torsion.mag.lo, torsion.diff)
  }
  result[,6:8]<-result[,3:5]/sqrt((result[,1]^2)+(result[,2]^2))
  invisible(result)
  return(result)
}
