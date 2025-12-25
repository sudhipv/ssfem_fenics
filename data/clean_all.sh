cd meshData
sh clean_meshData.sh
cd -

cd mallocData
sh clean_mallocData.sh
cd -

cd klePceData
sh clean_klePceData.sh
cd -

cd vtkOutputs
sh clean_vtkData.sh
cd -

cd ../external/dolfin/data
sh clean_fenicsData.sh
cd -

cd meshData
#rm PDE_Dim.txt
cd -

cd solution
#rm *.dat
cd -

