#!/bin/bash
  
param_1="0.0 0.1 0.3 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.45 1.5 1.55 1.6 1.7 1.8 1.9"

job_directory=bias
input_file=opt.fdf
job_file=run.slurm
template_folder=template

log=log_opt.out
data=data.out
rm $data

for ii in $param_1 ; do
                work_dir=$job_directory/${ii}
                
		echo $work_dir
		cp traj.sh $work_dir
		cd $work_dir
                ./traj.sh
                cd ../..

		 grep 'siesta: E_KS(eV)' $work_dir/$log | tail -n1 | awk '{print $4}' >> $data
done

wait
