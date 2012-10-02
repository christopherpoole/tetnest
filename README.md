Create Nested Tetrahedral Meshes for Use with CADMesh
=====================================================

[CADMesh](http://code.google.com/p/cadmesh/) is a C++ library for loading CAD based geometry into GEANT4 - a more through description can be found [here](http://eprints.qut.edu.au/53299/).
In addition to this basic function, the fast tessellated solid navigation technique described [here](http://eprints.qut.edu.au/52696/) has been implemented.
Sometimes it is useful to apply this GEANT4 navigation acceleration technique to nested tessellated mesh geometry however, the creation of the requisite nested tetrahedral meshes can be difficult.
This simple (and in places, rough) utility performs the tetrahedral nesting of a list of .ply surface meshes automatically.

Usage
-----
    # Get the code
    git clone https://github.com/christopherpoole/tetnest.git

    cd tetnest/mesh
    python ../tetnest.py inner.ply outer.ply

This will dump the required .ele and .node files to the current working directory (in this case ./tetnest/mesh).
These files are suitable for loading into GEANT4 with CADMesh as individual assemblies of tetrahedra with each assembly having unique material properties.
The outputted files combined.1\_-1\_.\* describe the outer shell (for this example):

![Outer Sphere](https://raw.github.com/christopherpoole/tetnest/master/example/outer.png)

and the files combined.1\_-1\_0\_.\* describe the inner sphere:

![Inner Sphere](https://raw.github.com/christopherpoole/tetnest/master/example/inner.png)

Requirements
------------

This scripts calls [tetgen](http://tetgen.berlios.de/) directly, so you need it in your PATH somewhere, if you want to view the generated meshes, have a look at [tetview](http://tetgen.berlios.de/tetview.html).

Caveats/TODO's
-------

* This script will only work with .ply surface mesh files - they must be surface meshes, not just point clouds
* Needs to be reformed as a proper Python module
