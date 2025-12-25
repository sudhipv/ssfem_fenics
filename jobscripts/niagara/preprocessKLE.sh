## Call KLE/PCE data generator and importer
## This take one command line input after number of dimensions (RVs)
## Enter the desired input and hit enter to run the code (or Ctr+C to exit)
cd ../data/klePceData/
gfortran KLE_PCE_Data_commandLineArg.F90
./a.out
cd -


## This code doesnot take any command line input everything must
## be specified in the code (hardcoding)
## Also, for special case, such as nDim=nOrd=KLEord = 1 :: Use this code
#cd ../data/klePceData/
#gfortran KLE_PCE_Data.F90
#./a.out
#cd -
