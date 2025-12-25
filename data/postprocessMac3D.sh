## KLE/PCE data generator and importer
## This take one command line input after number of dimensions (RVs)
cd vtkOutputs/

echo "1 - 3D Poisson"
echo "2 - 3D Elasticity"
read input1

if [ $input1 == 1 ]
then
    /Applications/MATLAB_R2017a.app/bin/./matlab -nojvm -nosplash -nodesktop -r writeVtk
elif [ $input1 ==  2 ]
then
    /Applications/MATLAB_R2017a.app/bin/./matlab -nojvm -nosplash -nodesktop -r writeVecVtk
else
echo "STOP: Wrong Input"
fi

exit
