# MGndeAdjointFoam
A neutron space-time kinetics solver based on OpenFOAM 10.

It is developed based on multi-region solver chtMutiRegionFoam, which has good geometric adaptability. Users can model and divide the meshes through business software. It should be noted that topological relation should be specified for the boundary surfaces of different regions to ensure the continuity of the meshes.

It includes a steady-state eigenvalue solvers and a transient solvers based on direct methods.

Due to the limited programming ability, the code is for reference only and will be further improved in the future.
