mkdir -p DATA
mkdir -p OUTPUT_FILES
cp src/Par_file DATA
cp src/SOURCE DATA
cp src/model.geo DATA
(cd DATA && gmsh model.geo -2 -o model.msh)
(cd DATA && python LibGmsh2Specfem_convert_Gmsh_to_Specfem2D_official.py model.msh -t F -l A -b A -r A)
xmeshfem2D
xspecfem2D