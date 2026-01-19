#!/bin/bash


ln -s ../../../cpp/sketches_SDL/Molecular/data
# ====== Paths

#name_def=test_BondAdaptedMesh    # optimal sampling grid around molecule

# ===  EFF

#name_def=test_eFF                 # Electron Force Field
#name_def=test_eFFMC               # compares EFF (eff) and CLCFGO (ff)
#name_def=test_eFF_old

# ===  Multi-Center EFF

#name_def=test_CLCFGO             # Compact linear combination of floating gaussian orbitals ( multi-center electron forcefield )
#name_def=test_CLCFSF             # Compact linear combination of finite spherical functions  ( multi-center electron forcefield )

# ===  Molecular-dynamics like

#name_def=test_ConfDynamics       # Sampling empty free space in collision potential
#name_def=test_SoftMolecularDynamics # Somehow not working now, large balls, large forces, why?

# ===  Reactive Force-fields

#name_def=test_FARFF              # GOOD, Flexible Atom sp-hybridization forcefield from cpp/common/molecular/FlexibleAtomReactiveFF.h and cpp/common/molecular/FlexibleAtomReactiveFF_dyn.h     
#name_def=test_PBD_LJ_cluster     # GOOD, Position-Based Dynamics with short-range approximate Lennard-Jones potential and cluster constraints

#name_def=test_RARFF              # Rigid Atom Reactive Force-field   "RARFF.h"
#name_def=test_RARFF2             # Rigid Atom Reactive Force-field 2 "RARFF2.h"
#name_def=test_RARFFarr           # Rigid Atom Reactive Force-field using array "RARFFarr.h"
name_def=test_RARFF_SR           # GOOD, large scale simulation of grid-based Bucket accelerated "RARFF_SR.h"

# ===  sp3 semi-quantum forcefield

#name_def=test_RspFF               # Interesting : Flexible Atom sp-hybridization forcefield, behaves strangely, defined in FspFFclean.h
#name_def=test_sp3space            # relaxation in sp3 space does not seem to converge
#name_def=test_spRotations          # Testing rotations of sp3 orbitals like Slater-Koster tables, user Approx::AutoApprox from cpp/common/Approx.h

# ======= NOT DOING SHIT

#name_def=test_MMFFmini             # EXIT with ERROR: nff.pairMask is not sorted => exit 
#name_def=test_EOFF                 # Demonstration of EOFF (Electron Orbital Force Field) from cpp/common/molecular/EOFF.h , does not render anything
#name_def=test_ESFF                 # Demonstration of ESFF (Electron Spin Foce Field)     from cpp/common/molecular/ESFF.h , does not render anything
#name_def=test_FTRFF                # Demonstration of FTRFF (Fixed Type Reactive Force-field) from cpp/common/molecular/FTRFF.h , does not render anything
#name_def=test_RRFF                 # Demonstration of RRFF (Reactive Reactive Force-field) from cpp/common/molecular/RRFF.h , does not render anything

# ======= NOT WORKING 

#name_def=test_Multipoles     


# Set $name to $1 if provided, otherwise to default_name
name="${1:-$name_def}"

echo $name

#dir=../../../cpp/Build/sketches_SDL/Molecular
dir=../../../cpp/Build/sketches_SDL/Molecular
#ln -s ../../../cpp/apps/OrbitalWar/data
ln -s ../../../cpp/common_resources

# ====== Multiprocesing
ncpu=`nproc`
ncpu=$(($ncpu - 1))     # let one CPU free for user interaction
echo "compile using ncpu="$ncpu
OMP_NUM_THREADS=$ncpu
export OMP_NUM_THREADS

# ====== Compilation
wd=`pwd`
cd $dir
pwd
rm $name
make -j$ncpu $name   # 2>$wd/compile_err.log
cd $wd
ln -s $dir/$name .

# ====== ASan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
echo   $LD_PRELOAD
export LD_PRELOAD

# ====== RUN

./$name
#./$name -s data/ship_ICF_interceptor_1.lua
#./$name -s data/ship_ICF_marksman_1.lua
#./$name -s data/ship_ICF_marksman_2.lua
#./$name -s data/ship_NFPP_pendulum_1.lua
#./$name -s data/ship_NTR_marksman_1.lua
#./$name -s data/colony_1.lua




