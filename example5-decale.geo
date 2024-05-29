// Gmsh project created on Tue Jul 26 15:40:16 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {10, 0, 0, 0.1};
//+
Circle(1) = {10, 0, 0, 3.5, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 32, 0, 2*Pi};
//+
Circle(3) = {0, 0, 0, 35, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2};
//+
Curve Loop(3) = {3};
//+
Plane Surface(1) = {1};
//+
Plane Surface(2) = {1,2};
//+
Plane Surface(3) = {2,3};
//+
Physical Surface("Coeur", 5) = {1};
//+
Physical Surface("Gaine", 6) = {2};
//+
Physical Surface("PML", 7) = {3};
//+
Physical Point("centre", 8) = {1};
//+
Transfinite Curve {1} = 100 Using Progression 1;
//+
Transfinite Curve {2} = 250 Using Progression 1;
//+
Transfinite Curve {3} = 350 Using Progression 1;
//+
Point{1} In Surface{1};
//+
Physical Curve("Bord", 9) = {3};
//+
Physical Point("Bord", 9) = {4};