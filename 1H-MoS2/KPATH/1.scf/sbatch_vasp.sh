#!/bin/sh
#SBATCH --job-name=x2_runv  --output=%x.o%j --error=%x.e%j
#SBATCH -p th1-2020-64 --time=72:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT --mail-user=h.kim@fz-juelich.de
#SBATCH --nodes=1  --tasks-per-node=64

  CURDIR=$SLURM_SUBMIT_DIR
  CURDIR_=`pwd | cut -d '/' -f 4-`
  SCRDIR=/SCRATCH/$CURDIR_
  SCRDIR_=/SCRATCH/$HOSTNAME/$CURDIR_
  CURDATE=`date '+%m-%d-%y %H:%M' `
  echo "Start VASP at $CURDATE jobid_${SLURM_JOBID}.running SUBMIT_DIR= $CURDIR SCR_DIR= $SCRDIR_" >> $HOME/job_status.dat
  echo "##    SCR_DIR  : $SCRDIR_"
  echo "##    WORK_DIR : $CURDIR "
  source compiler-select intel
 
 #VASP=$HOME/bin/vasp.5.4.4_noSOC 
  VASP=$HOME/bin/vasp.5.4.4_SOC_th64
  TBFIT=$HOME/code/bin/tbfit
  irun=0
  irun_scratch=0

  if [ $irun_scratch -eq 1 ] ; then
    srun  mkdir -p $SCRDIR
    srun  rm       $SCRDIR/*
    srun  cp *     $SCRDIR
          cd       $SCRDIR
          rm       ${SLURM_JOBNAME}.*${SLURM_JOBID}
    rm SCRATCH_*
    ln -s $SCRDIR_ $CURDIR/SCRATCH-$HOSTNAME
  fi

# 0: static ;1:opt ;2:opt + scf; 3: opt + scf + band; 4:scf+band ; 5: phonon ; 6: lattice ; 
# 8: scf + band (HSE, using WAVECAR) 9: opt + scf + band (HSE, using WAVECAR)
## every calculations should start in SCF directory, except for 0(staic), 5 calculations)

## 0. Static calculation
 if [ $irun -eq 0 ];then
 echo "## 0. Static calculation start" "irun = $irun, STATIC : ` date '+%m-%d-%y %H:%M'`"
  srun $VASP > log.out
# python unfolding.py

# rm CHG CHGCAR DOSCAR vasprun* XDAT* vasprun*
# gzip PARCHG # CHGCAR
# gzip PROCAR # CHGCAR
# cp SW* sw* $CURDIR
# cp WAVECAR $CURDIR
# cp OUTCAR PROCAR* PARCHG* log.out OSZICAR EIGENVAL CONTCAR IBZKPT CHGCAR* $CURDIR
# rm *
 echo "## 0. Static calculation end : ` date '+%m-%d-%y %H:%M'`"
 fi

## 1. OPTIMIZATION
 if [ $irun -eq 1 ];then
 echo "## 1. OPTIMIZATION start" "irun = $irun, OPTIMIZATION : ` date '+%m-%d-%y %H:%M'`"
 srun  $VASP > log.out
 mkdir OPTIMIZATION
 cp  vdw_kernel.bindat OUTCAR EI* POT* POS* INC* OS* CONT* IBZ* at* log* mpi* KPOINTS XDATCAR OPTIMIZATION
 echo "## 1. OPTIMIZATION end : ` date '+%m-%d-%y %H:%M'`"
 fi

## 2. OPT+SCF 
 if [ $irun -eq 2 ];then
 echo "## 1. OPTIMIZATION start" "irun = $irun, OPT+SCF : ` date '+%m-%d-%y %H:%M'`"
 srun $VASP > log.out
 mkdir OPTIMIZATION
 cp  vdw_kernel.bindat OUTCAR EI* POT* POS* INC* OS* CONT* IBZ* at* log* mpi* KPOINTS XDATCAR OPTIMIZATION
 echo "## 1. OPTIMIZATION end : ` date '+%m-%d-%y %H:%M'`"
 mv CONTCAR POSCAR
 #sed 's/ICHARG  =         2/ICHARG  =         1/g' INCAR > temp;mv temp INCAR
 sed 's/LCHARG  =        .FALSE./LCHARG  =        .TRUE. /g' INCAR > temp;mv temp INCAR
 sed 's/   EDIFFG/#  EDIFFG/g'                     INCAR > temp;mv temp INCAR
 sed 's/   NSW /#  NSW /g' INCAR > temp;mv temp INCAR
 sed 's/   IBRION /#  IBRION /g' INCAR > temp;mv temp INCAR
#sed 's/   ISYM    =         2/   ISYM    =        -1/g' INCAR > temp;mv temp INCAR
 sed 's/# EMIN /  EMIN /g'                 INCAR > temp;mv temp INCAR
 sed 's/# EMAX /  EMAX /g'                 INCAR > temp;mv temp INCAR
 sed 's/# NEDOS /  NEDOS /g'             INCAR > temp;mv temp INCAR
 sed 's/# LORBIT /  LORBIT /g'                   INCAR > temp;mv temp INCAR
 echo "## 2. SCF start : ` date '+%m-%d-%y %H:%M'`"
 srun $VASP > log.out
 rm vasprun* CHG
 gzip PROCAR DOSCAR CHGCAR
 echo "## 2. SCF end : ` date '+%m-%d-%y %H:%M'`"
 fi

## 3. OPT+SCF+BAND
 if [ $irun -eq 3 ];then
 echo "## 1. OPTIMIZATION start" "irun = $irun, OPT+SCF+BAND : ` date '+%m-%d-%y %H:%M'`"
 srun $VASP > log.out
 mkdir OPTIMIZATION
 cp vdw_kernel.bindat OUTCAR EI* POT* POS* INC* OS* CONT* IBZ* at* log* mpi* KPOINTS XDATCAR OPTIMIZATION
 echo "## 1. OPTIMIZATION end : ` date '+%m-%d-%y %H:%M'`"
 mv CONTCAR POSCAR
#sed 's/ICHARG  =         2/ICHARG  =         1/g' INCAR > temp;mv temp INCAR
 sed 's/LCHARG  =        .FALSE./LCHARG  =        .TRUE. /g' INCAR > temp;mv temp INCAR
 sed 's/   NSW /#  NSW /g' INCAR > temp;mv temp INCAR
 sed 's/   EDIFFG/#  EDIFFG/g'                     INCAR > temp;mv temp INCAR
 sed 's/   IBRION /#  IBRION /g' INCAR > temp;mv temp INCAR
#sed 's/   ISYM    =         2/   ISYM    =        -1/g' INCAR > temp;mv temp INCAR
#sed 's/# EMIN /  EMIN /g'                 INCAR > temp;mv temp INCAR
#sed 's/# EMAX /  EMAX /g'                 INCAR > temp;mv temp INCAR
#sed 's/# NEDOS /  NEDOS /g'             INCAR > temp;mv temp INCAR
 sed 's/# LORBIT /  LORBIT /g'                   INCAR > temp;mv temp INCAR
 echo "## 2. SCF start : ` date '+%m-%d-%y %H:%M'`"
 srun $VASP > log.out
 rm vasprun* CHG DOSCAR WAVECAR PCDAT XDATCAR
 echo "## 2. SCF end : ` date '+%m-%d-%y %H:%M'`"
 mkdir BAND
 cp vdw_kernel.bindat CHGCAR INCAR POSCAR POTCAR KPOINTS_BAND $CURDIR/BAND/
 cd $CURDIR/BAND/
 mv KPOINTS_BAND KPOINTS
 sed 's/ISTART  =         0/ISTART  =         1/g' INCAR > temp;mv temp INCAR
 sed 's/ICHARG  =         2/ICHARG  =        11/g' INCAR > temp;mv temp INCAR
#sed 's/  EMIN /# EMIN /g'                 INCAR > temp;mv temp INCAR
#sed 's/  EMAX /# EMAX /g'                 INCAR > temp;mv temp INCAR
#sed 's/  NEDOS /# NEDOS /g'             INCAR > temp;mv temp INCAR
#sed 's/LWAVE   =        .FALSE./LWAVE   =        .TRUE. /g' INCAR > temp;mv temp INCAR
 echo "## 3. BAND start : ` date '+%m-%d-%y %H:%M'`"
 srun $VASP > log.out
 rm CHG* vaspr* DOSCAR WAVECAR PCDAT XDATCAR ; gzip PRO*
 echo "## 3. BAND end : ` date '+%m-%d-%y %H:%M'`"
 cd $CURDIR
 gzip PROCAR DOSCAR CHGCAR
 fi

## 4. SCF+BAND
 if [ $irun -eq 4 ];then
 echo "## 1. SCF start" "irun = $irun, SCF+BAND : ` date '+%m-%d-%y %H:%M'`"
 srun $VASP > log.out
 rm vasprun* CHG
 echo "## 1. SCF end : ` date '+%m-%d-%y %H:%M'`"
 mkdir BAND ;  cp vdw_kernel.bindat CHGCAR INCAR POSCAR POTCAR KPOINTS_BAND $CURDIR/BAND/
 cd $CURDIR/BAND/  ;  mv KPOINTS_BAND KPOINTS
 sed 's/ISTART  =         0/ISTART  =         1/g' INCAR > temp;mv temp INCAR
 sed 's/ICHARG  =         2/ICHARG  =        11/g' INCAR > temp;mv temp INCAR
 sed 's/# LORBIT /  LORBIT /g'                   INCAR > temp;mv temp INCAR
#sed 's/  EMIN /# EMIN /g'                 INCAR > temp;mv temp INCAR
#sed 's/  EMAX /# EMAX /g'                 INCAR > temp;mv temp INCAR
#sed 's/  NEDOS /# NEDOS /g'             INCAR > temp;mv temp INCAR
#sed 's/LWAVE   =        .FALSE./LWAVE   =        .TRUE. /g' INCAR > temp;mv temp INCAR
 echo "## 2. BAND start : ` date '+%m-%d-%y %H:%M'`"
 srun $VASP > log.out
 rm CHG* vaspr* DOSCAR  ;  gzip PRO*
 echo "## 2. BAND end : ` date '+%m-%d-%y %H:%M'`"
 cd $CURDIR  ;  gzip PROCAR DOSCAR CHGCAR LOCPOT
 fi

## 5. PHONON (static)
 if [ $irun -eq 5 ];then
 echo "## 5. PHONON start, supercell with displacements:" "irun = $irun, PHONON  : ` date '+%m-%d-%y %H:%M'`"
 for i in 00{1..6}
 do
 echo "## $i th structure with the displacements.. : ` date '+%m-%d-%y %H:%M'`"
 mkdir disp-$i
 cp INCAR POSCAR-$i POTCAR KPOINTS atoms.d disp-$i/
 cd disp-$i
 mv POSCAR-$i POSCAR
 srun $VASP > log.out
 rm PROCAR DOSCAR WAVECAR CHGCAR CHG
 cd $CURDIR
 done
 echo "## 5. PHONON end : ` date '+%m-%d-%y %H:%M'`"
 gzip disp-*/vasprun.xml
 fi

## 6. LATT (static)
 if [ $irun -eq 6 ];then
 OUT_ENERGY="$CURDIR/ENERGY.dat"
 echo "## 6. LATT start : irun = $irun, LATT : ` date '+%m-%d-%y %H:%M'`"
 a0=4.0312
#c0=15.0
 printf " %s \n"   "## a0(initial guess)= $a0 Ang"         > $OUT_ENERGY
 printf " %s \n"   "## a/a0     a(Ang)      Energy(eV)"   >> $OUT_ENERGY
 for i in 0.94 0.96 0.98 1.00 1.02 1.04 1.06
 do
 echo "## a= a0 * $i ..optimization... : ` date '+%m-%d-%y %H:%M'`"
####### You shuould modify lattice vector elements accordingly ##
#    A11  A12   A13                                             #
#    A21  A22   A2                                              #
#    A31  A32   A33                                             #
 A11=`echo "$i * $a0" | bc -l`                                  #
#A12=`echo " 0      " | bc -l`                                  #
#A13=`echo " 0      " | bc -l`                                  #
 A21=`echo " 1 * $A11 / 2" | bc -l`                             #
 A22=`echo " 1 * $A21 * sqrt(3)" | bc -l`                       #
#A23=`echo " 0                   " | bc -l`                     #
#A31=`echo " 0                   " | bc -l`                     #
#A32=`echo " 0                   " | bc -l`                     #
#A33=`echo "     $c0             " | bc -l`                     #
 cp POSCAR_latt POSCAR_$i                                       #
 sed "s/A11/$A11/g" POSCAR_$i > temp;mv temp POSCAR_$i          #
 sed "s/A21/$A21/g" POSCAR_$i > temp;mv temp POSCAR_$i          #
 sed "s/A22/$A22/g" POSCAR_$i > temp;mv temp POSCAR_$i          #
#################################################################
 mkdir latt-$i
 cp INCAR POSCAR_$i POTCAR KPOINTS atoms.d latt-$i
 cd latt-$i
 mv POSCAR_$i POSCAR
 srun $VASP > log.out
 rm DOSCAR WAVECAR CHGCAR CHG vasprun.xml
 ENE=`grep "energy  without entropy=" OUTCAR | tail -1 | awk '{print $7}'`
 printf " %s %10.6f %10.6f \n"   "$i"  "$apa11"  "$ENE"   >> $OUT_ENERGY
 gzip *
 cd $CURDIR
 done
 echo "## 6. LATT end : ` date '+%m-%d-%y %H:%M'`"
 fi



## 7. NEB (static)
 if [ $irun -eq 7 ];then
 echo "## 7. NEB (static) start," "irun = $irun, NEB  : ` date '+%m-%d-%y %H:%M'`"
 for i in 0{0..8}
 do
 echo "## $i th structure with the displacements.. : ` date '+%m-%d-%y %H:%M'`"
 cp vdw_kernel.bindat INCAR POTCAR KPOINTS atoms.d $i/
 cd $i
 srun $VASP > log.out
 rm PROCAR DOSCAR WAVECAR CHGCAR CHG
 cd $CURDIR
 done
 echo "## 7. NEB (static) end  : ` date '+%m-%d-%y %H:%M'`"
 fi

## 8. SCF+BAND (HSE calculation, KPOINTS_HSE should be prepared in priori)
 if [ $irun -eq 8 ];then
 echo "## 1. SCF start " "irun = $irun, SCF+BAND : ` date '+%m-%d-%y %H:%M'`"
 srun $VASP > log.out
 echo "## 1. SCF end"
 rm vasprun* CHG

 mkdir BAND
 cp WAVECAR INCAR POSCAR POTCAR KPOINTS* $CURDIR/BAND/
 cd $CURDIR/BAND/
 cp KPOINTS_HSE KPOINTS
 sed 's/ISTART  =         0/ISTART  =         1/g' INCAR > temp;mv temp INCAR
 sed 's/ICHARG  =         2/ICHARG  =         1/g' INCAR > temp;mv temp INCAR
 sed 's/#  NELMIN  =        21/   NELMIN  =        21/g' INCAR > temp;mv temp INCAR
#sed 's/  EMIN /# EMIN /g'                 INCAR > temp;mv temp INCAR
#sed 's/  EMAX /# EMAX /g'                 INCAR > temp;mv temp INCAR
#sed 's/  NEDOS /# NEDOS /g'             INCAR > temp;mv temp INCAR
 sed 's/# LORBIT /  LORBIT /g'                   INCAR > temp;mv temp INCAR
 sed 's/LWAVE   =        .TRUE. /LWAVE   =        .FALSE./g' INCAR > temp;mv temp INCAR
 echo "## 2. BAND start : ` date '+%m-%d-%y %H:%M'`"
 srun $VASP > log.out
 rm CHG* vaspr* DOSCAR WAVECAR
 gzip PRO*
 echo "## 2. BAND end : ` date '+%m-%d-%y %H:%M'`"
 cd $CURDIR
 gzip PROCAR DOSCAR CHGCAR WAVECAR
 fi

## 9. OPT+SCF+BAND (HSE calculations, KPOINTS_HSE should be prepared in priori)
 if [ $irun -eq 9 ];then
 echo "## 1. OPTIMIZATION start" "irun = $irun, OPT+SCF+BAND : ` date '+%m-%d-%y %H:%M'`"
 srun $VASP > log.out
 mkdir OPTIMIZATION
 cp OUTCAR EI* POT* POS* INC* OS* CONT* IBZ* at* log* mpi* KPOINTS XDATCAR OPTIMIZATION
 echo "## 1. OPTIMIZATION end : ` date '+%m-%d-%y %H:%M'`"
 mv CONTCAR POSCAR
#sed 's/ICHARG  =         2/ICHARG  =         1/g' INCAR > temp;mv temp INCAR
#sed 's/LCHARG  =        .FALSE./LCHARG  =        .TRUE. /g' INCAR > temp;mv temp INCAR
 sed 's/   NSW /#  NSW /g' INCAR > temp;mv temp INCAR
 sed 's/   EDIFFG/#  EDIFFG/g'                     INCAR > temp;mv temp INCAR
 sed 's/   IBRION /#  IBRION /g' INCAR > temp;mv temp INCAR
#sed 's/   ISYM    =         2/   ISYM    =        -1/g' INCAR > temp;mv temp INCAR
 sed 's/   ISIF    =         3/#  ISIF    =         3/g' INCAR > temp;mv temp INCAR
#sed 's/# EMIN /  EMIN /g'                 INCAR > temp;mv temp INCAR
#sed 's/# EMAX /  EMAX /g'                 INCAR > temp;mv temp INCAR
#sed 's/# NEDOS /  NEDOS /g'             INCAR > temp;mv temp INCAR
 sed 's/LWAVE   =        .FALSE./LWAVE   =        .TRUE. /g' INCAR > temp;mv temp INCAR
 sed 's/# LORBIT /  LORBIT /g'                   INCAR > temp;mv temp INCAR
 echo "## 2. SCF start : ` date '+%m-%d-%y %H:%M'`"
 srun $VASP > log.out
 rm vasprun* CHG
 echo "## 2. SCF end : ` date '+%m-%d-%y %H:%M'`"
 mkdir BAND
 cp WAVECAR INCAR POSCAR POTCAR KPOINTS* $CURDIR/BAND/
 cd $CURDIR/BAND/
 cp KPOINTS_HSE KPOINTS
 sed 's/ISTART  =         0/ISTART  =         1/g' INCAR > temp;mv temp INCAR
 sed 's/ICHARG  =         2/ICHARG  =         1/g' INCAR > temp;mv temp INCAR
 sed 's/#  NELMIN  =        21/   NELMIN  =        21/g' INCAR > temp;mv temp INCAR
#sed 's/  EMIN /# EMIN /g'                 INCAR > temp;mv temp INCAR
#sed 's/  EMAX /# EMAX /g'                 INCAR > temp;mv temp INCAR
#sed 's/  NEDOS /# NEDOS /g'             INCAR > temp;mv temp INCAR
 sed 's/# LORBIT /  LORBIT /g'                   INCAR > temp;mv temp INCAR
 sed 's/LWAVE   =        .TRUE. /LWAVE   =        .FALSE./g' INCAR > temp;mv temp INCAR
 echo "## 3. BAND start : ` date '+%m-%d-%y %H:%M'`"
 srun $VASP > log.out
 rm CHG* vaspr* DOSCAR WAVECAR
 gzip PRO*
 echo "## 3. BAND end : ` date '+%m-%d-%y %H:%M'`"
 cd $CURDIR
 gzip PROCAR DOSCAR CHGCAR WAVECAR
 fi


  CURDATE=`date '+%m-%d-%y %H:%M' `
  SED_ARG=`echo "-i 's/jobid_${SLURM_JOBID}.running/jobid_${SLURM_JOBID}.end at $CURDATE/g'"` ; eval sed $SED_ARG $HOME/job_status.dat
