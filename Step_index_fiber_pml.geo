// Gmsh project created on Tue Jul 26 15:40:16 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.1};
//+
Circle(1) = {0, 0, 0, 4, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 7, 0, 2*Pi};
//+
Circle(3) = {0, 0, 0, 12, 0, 2*Pi};
//+
Circle(4) = {0, 0, 0, 15, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2};
//+
Curve Loop(3) = {3};
//+
Curve Loop(4) = {4};
//+
Plane Surface(1) = {1};
//+
Plane Surface(2) = {1,2};
//+
Plane Surface(3) = {2,3};
//+
Plane Surface(4) = {3,4};
//+
Physical Surface("Coeur", 5) = {1};
//+
Physical Surface("Gaine1", 6) = {2};
//+
Physical Surface("Gaine2", 7) = {3};
//+
Physical Surface("PML", 8) = {4};
//+
//Physical Curve("Interface", 6) = {1};
//+
Physical Curve("Bord", 7) = {4};
//+
Physical Point("Bord", 7) = {5};
//+
Physical Point("centre", 8) = {1};
//+
Transfinite Curve {1} = 60 Using Progression 1;
//+
Transfinite Curve {2} = 105 Using Progression 1;
//+
Transfinite Curve {3} = 180 Using Progression 1;
//+
Transfinite Curve {4} = 225 Using Progression 1;
//+
Point{1} In Surface{1};
//+
