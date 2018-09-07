#!/bin/bash -f

for ff in genericCube100x100x100.5051.pts line100000.5.pts \
    offsetCube79x79x79.1603.pts randomCube2000000.pts \
    big/vena_HARMINGGEE.pts.gz \
    big/PowerPlantM.pts.gz \
    big/Barcelona_fullmodel_16m_20M_ascii.pts.gz  \
    big/Barcelona_fullmodel_8m_100M_ascii.pts.gz  \
    big/Barcelona_4_370M_ascii.pts.gz  \
    big/Barcelona_4_370MFluid_ascii.pts.gz  \
    big/ConfidencialBoeing777_all.pts.gz
do 
    echo "Doing $ff ..."
    pp=`../build/SearchLab.exe -s -n 1 $ff |grep " gnuplot "`
    $pp
    echo "... $pp done"
done
