#!/bin/bash
  
input=("log_2-1_dft_wfn.out" "log_2-2_dft_wfn.out" "log_3-1_0V.out" "log_3-2_0V.out")
output=("data_2-1.out" "data_2-2.out" "data_3-1.out" "data_3-2.out")

for folder in */ ; do
                work_dir=${folder}
                echo $work_dir

                for k in "${!input[@]}"; do

                        atoms="$(grep 'NumberOfAtoms' $work_dir/${input[$k]} | awk '{print $2}')"
                        lines=$((atoms+2))
                        echo $lines

                        rm $work_dir/${output[$k]}

			grep -A $lines 'siesta: Atomic forces' $work_dir/${input[$k]} >> $work_dir/temp.out
			grep -B $lines 'siesta:  Tot' $work_dir/${input[$k]} >> $work_dir/temp.out
                        tail -n $lines $work_dir/temp.out > $work_dir/${output[$k]}
        	 	sed -i -e :a -e '$d;N;1,2ba' -e 'P;D'  $work_dir/${output[$k]}
	done
done


