PREFIX_VALID=validation
EXE=../../build/SearchLab.exe
files=" \
    ../genericCube100x100x100.5051.pts                  \
    ../genericCube10x10x10.55556.pts		        \
    ../genericCube10x10x10_removedPoints.55556.pts	\
    ../genericCube2x2x2.500000.pts			\
    ../line100000.5.pts				        \
    ../offsetCube79x79x79.1603.pts			\
    ../randomCube2000000.pts                            \
"
files_big=" \
    ../big/Barcelona_4_370M_ascii.pts.gz		\
    ../big/Barcelona_4_370MFluid_ascii.pts.gz	        \
    ../big/Barcelona_fullmodel_16m_20M_ascii.pts.gz	\
    ../big/Barcelona_fullmodel_8m_100M_ascii.pts.gz	\
    ../big/ConfidencialBoeing777_all.pts.gz		\
    ../big/ConfidencialBoeing777_all_sorted.pts.gz	\
    ../big/PowerPlantM.pts.gz			        \
    ../big/PowerPlantM_sorted.pts.gz		        \
    ../big/PowerPlantM_sorted_unique.pts.gz		\
    ../big/vena_HARMINGGEE.pts.gz                       \
"
for ff in $files $files_big
do
    xpath=${ff%/*}
    xbase=${ff##*/}
    xfext=${xbase##*.}
    xpref=${xbase%.*}
    # xpref=$(basename "${ff%.*}")
    # eventually adde date: $(date +%F_%T
    outfile=${PREFIX_VALID}_${xpref}.val
    # filter out running times which may vary...
    echo creating validation file ${outfile}
    # $EXE -type bin -n 1 $ff | grep -v Times | grep -v Creating | tee "${outfile}"
    $EXE -type bin -n 1 $ff | grep -v Times | grep -v Creating > "${outfile}"
done
