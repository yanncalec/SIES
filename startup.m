cd +tools;
mex -largeArrayDims fulltosparse.c;
mex -largeArrayDims tensorplus_mex.c;
mex -largeArrayDims tensorprod_mex.c;

cd ../+wavelet/@OrthoWvl;
% mex Green2D_apprxcoeff_mex.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp";
% mex Green2D_detailcoeff_mex.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp";
% mex monomcoeff_mex.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp";
mex -largeArrayDims lookup_table_wvl_mex.c;

cd ../../;
