cd makefiles
pwd
for filename in *; do
    cd $filename
    make clean -C ../../ && make -C ../../
    cd ..
done