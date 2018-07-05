cd Input
for %%X in (*.txt) do "C:\Program Files\VCG\MeshLab\meshlabserver.exe" -i "%%X" -o "../Output/%%X.stl" -s ../meshclean.mlx

pause