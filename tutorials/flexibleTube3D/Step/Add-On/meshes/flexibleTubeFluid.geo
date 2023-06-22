Include "flexibleTubeParameters.geo";

blHeight = innerR*1.0/4.0;
unstrR = 0.5*(innerR-blHeight); // unstructured node radius


pp = 1/(2*(-Tan(15*Pi/180) + 1)^2);
a = 0.5 - pp;
b = 2*pp + 1;
c = -pp - 0.5;
ratio = (-b + Sqrt(b*b - 4*a*c))/(2*a);


DeltaY1 = (innerR-blHeight-unstrR)/nElemEdge;
ratioH1_0 = (blHeight - DeltaY0)/(blHeight - DeltaY1);
nElemH1 = Ceil(Log(ratioH1_0*DeltaY1/DeltaY0)/Log(ratioH1_0))*refinementFactor;
ratioH1 = Exp((1/(nElemH1-1))*Log(DeltaY1/DeltaY0));

//cRadius = unstrR/ratio; //120 deg between 3 unstructured edges
cRadius = (innerR-blHeight)*(innerR-blHeight)/unstrR; // linear curvature variation

Point(1) = {0, 0, 0};
Point(2) = {innerR, 0, 0};
Point(3) = {0, innerR, 0};
Point(4) = {0, unstrR, 0};
Point(5) = {unstrR, 0, 0};
Point(6) = {0, unstrR-cRadius, 0};
Point(7) = {unstrR-cRadius, 0, 0};

r = cRadius - unstrR;
x = -0.5*r + 0.5*Sqrt(2*cRadius*cRadius - r*r);
Point(9) = {x, x, 0};
Point(10) = {innerR*0.5*Sqrt(2), innerR*0.5*Sqrt(2), 0};

Point(11) = {innerR-blHeight, 0, 0};
Point(12) = {0, innerR-blHeight, 0};
Point(13) = {(innerR-blHeight)*0.5*Sqrt(2), (innerR-blHeight)*0.5*Sqrt(2), 0};

Circle(1) = {4, 6, 9};
Circle(2) = {5, 7, 9};
Circle(3) = {3, 1, 10};
Circle(4) = {10, 1, 2};
Line(5) = {1, 5};
Line(6) = {1, 4};

Line(7) = {4, 12};
Line(8) = {12, 3};
Line(9) = {9, 13};
Line(10) = {13, 10};
Line(11) = {5, 11};
Line(12) = {11, 2};
Circle(13) = {12, 1, 13};
Circle(14) = {13, 1, 11};
Line Loop(15) = {6, 1, -2, -5};
Plane Surface(16) = {15};
Line Loop(17) = {1, 9, -13, -7};
Plane Surface(18) = {17};
Line Loop(19) = {9, 14, -11, 2};
Plane Surface(20) = {19};
Line Loop(21) = {14, 12, -4, -10};
Plane Surface(22) = {21};
Line Loop(23) = {10, -3, -8, 13};
Plane Surface(24) = {23};
Recombine Surface {16, 24, 18, 20, 22};
Transfinite Line {5, 1, 13, 3, 6, 2, 14, 4, 11, 9, 7} = nElemEdge*refinementFactor+1 Using Progression 1;
Transfinite Line {-12,-10,-8} = nElemH1+1 Using Progression ratioH1;
Transfinite Surface {16};
Transfinite Surface {18};
Transfinite Surface {20};
Transfinite Surface {22};
Transfinite Surface {24};

Geometry.CopyMeshingMethod = 1;

Symmetry {1, 0, 0, 0} {
  Duplicata { Point{1}; Point{13}; Point{12}; Point{11}; Point{10}; Point{9}; Point{5}; Point{4}; Point{3}; Point{2}; Line{5}; Line{14}; Line{13}; Line{12}; Line{11}; Line{10}; Line{9}; Line{8}; Line{7}; Line{6}; Line{4}; Line{3}; Line{2}; Line{1}; Surface{24}; Surface{22}; Surface{20}; Surface{18}; Surface{16}; }
}

Symmetry {0, 1, 0, 0} {
  Duplicata { Point{19}; Point{2}; Point{3}; Point{4}; Point{5}; Point{9}; Point{10}; Point{11}; Point{12}; Point{13}; Point{15}; Point{17}; Point{18}; Point{1}; Point{20}; Point{23}; Line{36}; Line{11}; Line{12}; Line{13}; Line{14}; Line{25}; Line{26}; Line{27}; Line{28}; Line{29}; Line{30}; Line{31}; Line{35}; Line{10}; Line{37}; Line{38}; Line{1}; Line{2}; Line{3}; Line{4}; Line{5}; Line{6}; Line{7}; Line{8}; Line{9}; Surface{16}; Surface{18}; Surface{20}; Surface{22}; Surface{24}; Surface{39}; Surface{44}; Surface{49}; Surface{54}; Surface{59}; }
}

Extrude {0, 0, tubeLength} {
  Point{17}; Point{86}; Point{84}; Point{83}; Point{81}; Point{80}; Point{78}; Point{23}; Point{20}; Point{19}; Point{18}; Point{87}; Point{15}; Point{13}; Point{12}; Point{11}; Point{10}; Point{9}; Point{5}; Point{4}; Point{3}; Point{2}; Point{1}; Point{90}; Point{88}; Line{70}; Line{82}; Line{81}; Line{79}; Line{78}; Line{77}; Line{76}; Line{75}; Line{74}; Line{73}; Line{72}; Line{71}; Line{83}; Line{67}; Line{66}; Line{64}; Line{63}; Line{60}; Line{38}; Line{37}; Line{36}; Line{35}; Line{31}; Line{30}; Line{84}; Line{9}; Line{29}; Line{28}; Line{27}; Line{26}; Line{25}; Line{14}; Line{13}; Line{12}; Line{11}; Line{10}; Line{8}; Line{7}; Line{6}; Line{5}; Line{4}; Line{3}; Line{2}; Line{1}; Surface{59}; Surface{130}; Surface{125}; Surface{120}; Surface{115}; Surface{110}; Surface{105}; Surface{100}; Surface{95}; Surface{90}; Surface{85}; Surface{54}; Surface{49}; Surface{44}; Surface{39}; Surface{24}; Surface{22}; Surface{20}; Surface{18}; Surface{16}; Layers{nElemAxial}; Recombine;
}

Physical Surface("INLET") = {16, 85, 95, 20, 22, 100, 105, 90, 125, 110, 115, 120, 130, 59, 18, 24, 39, 54, 49, 44};

Physical Surface("OUTLET") = {661, 683, 705, 639, 617, 595, 353, 749, 771, 727, 507, 529, 573, 375, 419, 441, 463, 397, 551, 485};

Physical Surface("INTERFACE_FLUID") = {239, 323, 319, 243, 171, 175, 227, 199};

Physical Volume("volume") = {16, 15, 19, 12, 17, 18, 14, 20, 13, 1, 11, 9, 2, 8, 4, 5, 10, 3, 7, 6};
