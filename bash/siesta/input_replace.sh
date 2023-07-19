#!/bin/bash

files=("1_bulkLR.fdf"  "2-1_dft_wfn.fdf"  "2-2_dft_wfn.fdf"  "3-1_0V.fdf"  "3-2_0V.fdf")
input="MaxSCFIterations*"
output="MaxSCFIterations                10000"

for folder in */; do
                work_dir=${folder}
                echo $work_dir

                for k in "${!files[@]}"; do

			echo ${files[$k]}
			sed -i -e "s/$input/$output/g"  $folder/input/${files[$k]}
        done
done

