# FEM---Heat-flux-calculations
This repo follows the thesis "Numerical calculation of the heat transfer through sleeping bag materials" and contain the files needed to calculate the heat flux through Sample 1 and Sample 2 for two different geometries. The testing procedures for children's sleeping bags have yet to be standardized. This is due to significant physiological differences, which make it difficult to establish a thermophysiological model for defining the utility range in the same way the ISO standard does for adult sleeping bags.
This thesis aims to explore different modes of heat transfer and material properties that are essential for understanding the sleeping bag system as a whole. Furthermore, it introduces a numerical finite element model for simulating dry heat loss, intended as a supplementary tool to be used alongside to the development of a new testing standard for children's sleeping bags. Two material samples are examined. Sample 1 is a 7.5 mm thick typical synthetic sleeping bag fabric, while sample 2 is a 5 mm thick neoprene material. 

Most meshes we have worked with symbolize only the sleeping bag material and will contain an inner face representing the thermal manikin and an outer face representing the exterior surface of the sleeping bag. It is not necessary to include the layer of air outside the bag as the surface effect from radiation and convection will be handled as boundary conditions. We introduce two geometries. first is a cylinder(in this repo called InfCyl) that replicate the shape of a thermal manikin. The second is a rounded cylinder mesh (Pill) that is also explored.

To utilize this repo, clone it to your and download the requirements file by typing the command below in the terminal.

pip install -r requirements.txt


Most files contain the exploration of a specific sample and shape. If you want to modify this, you can either use the shape proposed in the file and change the parameters, or you can upload your own mesh explore. In the latter case, just make sure to upload a mesh as a .msh file. Most meshes can easily be converted to this by importing it to Gmsh and exporting as a .msh.
