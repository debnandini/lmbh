if [ "$#" -ne 1 ]; then
    echo "Request segment length in  days"
    echo "run.sh X"
fi
gcc -o segmentSangria segmentSangria.c -lhdf5 -lgsl
clang -Xpreprocessor -I/opt/local/include/libomp -L/opt/local/lib/libomp -lomp -w -o SpecFit SpecFit.c Utilities.c -lgsl -lm -O3
clang -Xpreprocessor -I/opt/local/include/libomp -L/opt/local/lib/libomp -lomp -w -o PTMCMC PTMCMC.c Utils.c Response.c IMRPhenomD_internals.c IMRPhenomD.c -lgsl -lgslcblas -lm -O3

#get data segment
./segmentSangria $1

#get PSD for E channel
./SpecFit AET_seg0_t.dat 1; 
mv specfit.dat specfit_1_0.dat; 

#get PSD for A channel
./SpecFit AET_seg0_t.dat 0; 
mv specfit.dat specfit_0_0.dat

#run Sampler
./PTMCMC 0 0
