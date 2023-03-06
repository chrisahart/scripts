#! /bin/bash

for d in */ ; do
    echo "$d"
    cd $d
    pot.sh
    cd ..
done
