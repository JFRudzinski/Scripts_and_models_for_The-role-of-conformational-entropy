#!/bin/bash

#####################################################################

#System Specific
sys="aaqaa"

blen=0.6
tol=5.0

######################################################################

#Directories
bin="/gpfs/work/jfr148/pkg/gromacs-4.5.3-dp/bin"
suff="-xj_sn-4"

######################################################################

#Executables
g_rdf="${bin}/g_rdf${suff}"
g_bond="${bin}/g_bond${suff}"
g_angle="${bin}/g_angle${suff}"
g_rmsd="${bin}/g_rmsdist${suff}"
g_rmsf="${bin}/g_rmsf${suff}"

#g_rdf="g_rdf4.5.3"
#g_bond="g_bond4.5.3"
#g_angle="g_angle4.5.3"
#g_rmsd="g_rmsdist4.5.3"
#g_rmsf="g_rmsf4.5.3"

######################################################################

#File Names
tpr="${sys}.tpr"
trr="${sys}.xtc"
ndx="index_dist_all1-4bonds.ndx"
native="${sys}.gro"

######################################################################

#Calculate Distributions

#bonds
frames=50001

N14=12
for i in `seq 1 ${N14}`;
do
    g_bond -f ${trr} -s ${tpr} -blen ${blen} -tol ${tol} -n ${ndx} -d distance_bonds_$i-$(($i + 3)).xvg -o bonds_$i-$(($i + 3)).xvg <<-EOF
           $(($i - 1))
EOF

    tail -n ${frames} distance_bonds_$i-$(($i + 3)).xvg &> tmp.xvg
    mv tmp.xvg distance_bonds_$i-$(($i + 3)).xvg
done

g_gyrate -f ${trr} -s ${tpr} -n index_types.ndx -o Rg.xvg <<-EOF
       0
EOF

tail -n ${frames} Rg.xvg &> tmp.xvg
mv tmp.xvg Rg.xvg

g_rama -f ${trr} -s ../../aaqaa.for-rama.tpr -o rama.xvg

grep "ALA-2" rama.xvg &> rama_ALA2.xvg
grep "ALA-3" rama.xvg &> rama_ALA3.xvg
grep "GLN-4" rama.xvg &> rama_GLN4.xvg
grep "ALA-5" rama.xvg &> rama_ALA5.xvg
grep "ALA-6" rama.xvg &> rama_ALA6.xvg
grep "ALA-7" rama.xvg &> rama_ALA7.xvg
grep "ALA-8" rama.xvg &> rama_ALA8.xvg
grep "GLN-9" rama.xvg &> rama_GLN9.xvg
grep "ALA-10" rama.xvg &> rama_ALA10.xvg
grep "ALA-11" rama.xvg &> rama_ALA11.xvg
grep "ALA-12" rama.xvg &> rama_ALA12.xvg
grep "ALA-13" rama.xvg &> rama_ALA13.xvg
grep "GLN-14" rama.xvg &> rama_GLN14.xvg
grep "ALA-15" rama.xvg &> rama_ALA15.xvg
grep "ALA-16" rama.xvg &> rama_ALA16.xvg

