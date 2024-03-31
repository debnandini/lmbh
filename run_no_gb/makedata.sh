#!/bin/sh

#make full data
gcc -o segmentSangria_t segmentSangria_t.c -lhdf5 -lgsl
./segmentSangria_t 37.0 
#rm AET_f.dat AET_t.dat AET_seg1*.dat AET_seg2*.dat AET_seg3*.dat AET_seg4*.dat AET_seg5*.dat AET_seg6*.dat AET_seg7*.dat AET_seg8*.dat AET_seg9*.dat

#back up full data
mkdir bk/
cp ./AET*.dat ./bk/
rm AET_seg*_f*.dat

#seg_dgb
mkdir seg_dgb
mv ./segmentSangria_dgb.c ./seg_dgb/
mv ./LDC2_sangria_training_v2.h5 ./seg_dgb/
cp Constants.h Header.h ./seg_dgb
cd ./seg_dgb
gcc -o segmentSangria_dgb segmentSangria_dgb.c -lhdf5 -lgsl
./segmentSangria_dgb 37.0
#rm AET_f.dat AET_t.dat AET_seg1*.dat AET_seg2*.dat AET_seg3*.dat AET_seg4*.dat AET_seg5*.dat AET_seg6*.dat AET_seg7*.dat AET_seg8*.dat AET_seg9*.dat
cd ../

#seg_vgb
mkdir seg_vgb
mv ./segmentSangria_vgb.c ./seg_vgb/
mv ./seg_dgb/LDC2_sangria_training_v2.h5 ./seg_vgb/
cp Constants.h Header.h ./seg_vgb/
cd ./seg_vgb
gcc -o segmentSangria_vgb segmentSangria_vgb.c -lhdf5 -lgsl
./segmentSangria_vgb 37.0
#rm AET_f.dat AET_t.dat AET_seg1*.dat AET_seg2*.dat AET_seg3*.dat AET_seg4*.dat AET_seg5*.dat AET_seg6*.dat AET_seg7*.dat AET_seg8*.dat AET_seg9*.dat
cd ../

#seg_igb
mkdir seg_igb
mv ./segmentSangria_igb.c ./seg_igb/
mv ./seg_vgb/LDC2_sangria_training_v2.h5 ./seg_igb/
cp Constants.h Header.h ./seg_igb/
cd ./seg_igb
gcc -o segmentSangria_igb segmentSangria_igb.c -lhdf5 -lgsl
./segmentSangria_igb 37.0
#rm AET_f.dat AET_t.dat AET_seg1*.dat AET_seg2*.dat AET_seg3*.dat AET_seg4*.dat AET_seg5*.dat AET_seg6*.dat AET_seg7*.dat AET_seg8*.dat AET_seg9*.dat
cd ../

python seg.py AET_seg*_t.dat
mv ./seg_igb/LDC2_sangria_training_v2.h5 ./
gcc -o segmentSangria_f segmentSangria_f.c -lhdf5 -lgsl
./segmentSangria_f 37.0
