Include "flexibleTubeParameters.geo";

outerR = innerR + tubeThickness;

Geometry.CopyMeshingMethod = 1;

Point(1) = {0, 0, 0};
Point(2) = {innerR, 0, 0};
Point(3) = {0, innerR, 0};
Point(4) = {outerR, 0, 0};
Point(5) = {0, outerR, 0};

Circle(1) = {2, 1, 3};
Circle(2) = {5, 1, 4};
Line(3) = {3, 5};
Line(4) = {2, 4};
Line Loop(5) = {3, 2, -4, 1};
Plane Surface(6) = {5};

Transfinite Line {3, 4} = refinementFactor*nElemThickness+1 Using Progression 1;
Transfinite Line {1, 2} = 2*refinementFactor*nElemEdge+1 Using Progression 1;
Transfinite Surface {6};
Recombine Surface {6};


Symmetry {1, 0, 0, 0} {
  Duplicata { Point{5}; Point{3}; Line{3}; Line{2}; Line{1}; Surface{6}; Point{4}; Point{2}; Line{4}; }
}

Symmetry {0, 1, 0, 0} {
  Duplicata { Point{5}; Point{3}; Line{8}; Line{3}; Line{9}; Line{1}; Line{2}; Surface{10}; Surface{6}; Point{2}; Point{4}; Line{4}; Point{14}; Point{17}; Line{13}; }
}

Extrude {0, 0, tubeLength} {
  Point{5}; Point{3}; Line{8}; Line{3}; Line{9}; Line{1}; Line{2}; Surface{10}; Surface{6}; Point{14}; Point{17}; Line{13}; Line{16}; Line{14}; Surface{19}; Point{19}; Point{18}; Line{15}; Line{17}; Line{18}; Surface{24}; Point{2}; Point{4}; Line{4}; Layers{nElemAxial}; Recombine;
}

Physical Surface("INLET_SOLID") = {6, 10, 19, 24};

Physical Surface("OUTLET_SOLID") = {68, 90, 150, 120};

Physical Surface("INTERFACE_SOLID") = {124, 42, 38, 94};

Physical Surface("OUTER_SOLID") = {30, 46, 98, 128};

Physical Volume("volume") = {1, 2, 3, 4};
