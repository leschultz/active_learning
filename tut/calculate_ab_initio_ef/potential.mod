# NOTE: This script can be modified for different pair styles 
# See in.elastic for more info.


# ---------- Define Interatomic Potential --------------------- 
pair_style eam/alloy
pair_coeff * * Farkas.eam.alloy Nb

 
# Set neighbor style
neighbor   2.0 bin 
#neigh_modify every 1 delete 0 check yes 

# Setup minimization style
min_style	     cg
min_modify	     dmax ${dmax} line quadratic



