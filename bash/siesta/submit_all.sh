#! /bin/bash

for folder in */ ; do
        echo "$folder"
        cd $folder
        qsub run.slurm
        cd ..
done
