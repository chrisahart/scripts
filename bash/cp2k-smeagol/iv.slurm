#!/usr/bin/env bash

#SBATCH -n 8
#SBATCH -C  "scarf21"
#SBATCH -t 06:00:00

module purge
module load GCC/10.3.0 OpenMPI/4.1.1-GCC-10.3.0 OpenBLAS/0.3.15-GCC-10.3.0

# Sergey branch
#cp2k=/home/vol07/scarf1157/code/smeagol/cp2k-smeagol/scarf21/smgl/cp2k-private/exe/local/cp2k.psmp

# Chris branch 1
#cp2k=/home/vol07/scarf1157/code/smeagol/cp2k-smeagol/scarf21/smgl/cp2k-private-chris-1/exe/local/cp2k.psmp
cp2k=/home/vol07/scarf1157/code/smeagol/cp2k-smeagol/scarf21/smgl/cp2k-private-chris-1/exe/local/cp2k.pdbg

find . -maxdepth 1 ! -name 'run.slurm' -type f -exec rm  {} +
export OMP_NUM_THREADS=1
cp input/* .

mpirun $cp2k  -i 1_bulkLR.inp -o log_1-bulkLR.out 
mpirun $cp2k  -i 2_dft_wfn.inp -o log_2_dft_wfn.out 
#mpirun $cp2k  -i 3_0V.inp -o log_3_0V.out 

rm IV.log; rm VI.log
echo "IV curve" > IV.log
for ii in $(seq -2 0.1 2) ; do
        echo ${ii} >> VI.out
        sed -e "s/PROJECT_REPLACE/V-${ii}/g" 4_V.inp > V_${ii}.inp
        sed -i -e "s/V_REPLACE/${ii}/g"  V_${ii}.inp
        mpirun $cp2k -i V_${ii}.inp -o log_V_${ii}.out
        cat V-${ii}.CUR >> IV.log
done
echo "done"  >> VI.out
