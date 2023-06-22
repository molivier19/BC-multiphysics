refinementFactor = 2;
nElemEdge = 2;
nElemThickness = 4; // nElem in the thickness direction of the outer shell

innerR = 0.005;
tubeLength = 0.05;
tubeThickness = 0.001;

DeltaY0 = 0.04*innerR; // inner tube first layer (at wall) thickness
nElemAxial = nElemEdge*refinementFactor*Round(tubeLength/(0.25*innerR*Pi));

