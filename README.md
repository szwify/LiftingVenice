# PoroFEML-2D -   Matlab Finite Element routines for 2D poromechanics

This is a set of matlab routines for 1D (axi-symmetric only) and 2D (plane-strain & axi-symmetric) Finite Element analysis for elastic & coupled problems.

Weak form operators currently coded:
- Elasticity: stiffness,  pre stress field , boundary loads
- Laplacian
- Mass matrix
- poroelastic Coupling operator 

Problem type currently coded up:
- 2D plane-strain 
- 2D axi-symmetry

Available Element type:
- Linear 1D element : Seg2
- Linear Quadrilateral Element (2D) : Qua4
- Linear triangular Element (2D) : Tri3

The code uses operator overloading to assemble system matrices, vectors etc. At the high level, it allows to script the problem solution relatively easily. 

Test Examples can be found under /Tests

For mesh generation, either ones directly script it in matlab (for simple geometries), or use available matlab libraries e.g. MESH2D - Delaunay-based unstructured mesh-generation 
(https://ch.mathworks.com/matlabcentral/fileexchange/25555-mesh2d-delaunay-based-unstructured-mesh-generation)

Short description of the different classes 


FEmesh : define a mesh (FEnode object containing nodal coordinates) + connectivity table + material ID table

 

Operator classes: series of classes allowing the integration of element level contribution depending on the PDE operator (elasticity, laplacian, etc. , load vector, etc.) 
- Operator_Elasticity
- Operator_Laplacian
- Operator_Mass
- Operator_Load_Boundary
- Operator_Load_InitialStresses
- Operator_Load_FluxSource

Element classes - contained the B-matrices, isoparametric mapping etc.
- Element_Manifold_1 (for 1D element) with derived class
  + Element_Seg2
- Element_Manifold_2 (for 2D element) with derived class
  + Element_Qua4
  + Element_Tri3
 