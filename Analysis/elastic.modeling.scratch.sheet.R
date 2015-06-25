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


SZ15.lvl4.compliance<-matrix(c(7.02, 3.30, 5.87,0,0,0,3.30,7.02,5.87,0,0,0,5.87,5.87,22.17,0,0,0,0,0,0,6.6,0,0,0,0,0,0,6.6,0,0,0,0,0,0,3.72), nrow=6)

SZ15.lvl4.stiffness<-solve(SZ15.lvl4.compliance)

