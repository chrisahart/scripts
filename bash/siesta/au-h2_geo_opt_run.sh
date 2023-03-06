#!/bin/bash
  
#param_1="0.0 0.1 0.3 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.45 1.5 1.55 1.6 1.7 1.8 1.9"
#job_directory=bias
#input_file=opt.fdf
#job_file=run.slurm
#template_folder=template

param_1="1.7 1.8"
job_directory=bias-rs
input_file=opt.fdf
job_file=run.slurm
template_folder=template
geom_dir=/home/mmm1169/Scratch/work/transport/AuH2/siesta-smeagol/geo_opt/bias/1.9
lead='%block AtomicCoordinatesAndAtomicSpecies'
tail='%endblock AtomicCoordinatesAndAtomicSpecies'

for ii in $param_1 ; do
                work_dir=$job_directory/${ii}
                if [ -d $work_dir ] ; then
                        rm -r $work_dir
                fi

		cp -r $template_folder $work_dir

                sed -i -e "s/V_REPLACE/${ii}/g" \
                        $work_dir/input/$input_file

		geom_file=$geom_dir/final.siesta
                echo $geom_file
		sed -i -e "/$lead/,/$tail/{ /$lead/{p; r $geom_file
                                }; /$tail/p; d }"  $work_dir/input/$input_file
done

for ii in $param_1 ; do
                work_dir=$job_directory/${ii}
                
		cd $work_dir
                qsub $job_file
                cd ../..
done

wait
