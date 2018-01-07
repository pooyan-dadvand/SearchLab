export EXE=/c/newGiD/SearchLab/build/Release/SearchLab.exe.exe
for ff in genericCube100x100x100.5051.pts line100000.5.pts offsetCube79x79x79.1603.pts randomCube2000000.pts \
    big/Barcelona_4_370M_ascii.pts big/Barcelona_4_370MFluid_ascii.pts \
    big/Barcelona_fullmodel_16m_20M_ascii.pts big/Barcelona_fullmodel_8m_100M_ascii.pts \
    big/PowerPlantM.pts big/vena_HARMINGGEE.pts
do
    echo ===== $ff ===== | tee -a cases.txt
    $EXE -n 1 -s $ff | grep -C 1 -i occup|grep -v "=== Bins"|grep -v time | tee -a cases.txt
done;
