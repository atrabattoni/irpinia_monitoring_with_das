sh clean.sh
mkdir DATA
mkdir OUTPUT_FILES
cp src/Par_file DATA
cp src/SOURCE DATA
cp src/interfaces.dat DATA
xmeshfem2D
xspecfem2D
