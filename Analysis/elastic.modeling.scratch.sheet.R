# scrap of script for stress-strain analytical modeling of collagen fiber orientation

# code plan:
# 1. per-pixel stress tensor from cross-sectional properties and length, Voigt notation (6x1)
# 2. rotate per-pixel compliance matrix to align with CFO (dip direction ambiguity?)
#    2a. isotropic reduction of empirical compliance matrix
#    2b. empirical linear orthotropic compliance matrix
# 3. calculate per-pixel strain tensor for each state
# 4. grayscale:
#    a. predicted normal strain intensity (iso, ortho, diff)
#    b. predicted radial strain intensity (iso, ortho, diff)
#    c. predicted circumferential strain intensity (iso, ortho, diff)
#    d. predicted principal strain intensity (iso, ortho, diff)
# 5. false-color predicted principal strain orientation (iso, ortho, diff)
# 6. line probes to help out the dichromats

SZ15.lvl4.compliance<-matrix(c(7.02, 3.30, 5.87,0,0,0,3.30,7.02,5.87,0,0,0,5.87,5.87,22.17,0,0,0,0,0,0,6.6,0,0,0,0,0,0,6.6,0,0,0,0,0,0,3.72), nrow=6)

SZ15.lvl4.stiffness<-solve(SZ15.lvl4.compliance)

# for a uniform modeled distribution of strain