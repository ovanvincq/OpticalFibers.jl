// Gmsh project created on Tue Jul 26 15:40:16 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.05};
//+
Circle(1) = {0, 0, 0, 2};
//+
Circle(2) = {0, 0, 0, 10};
//+
Delete { Point{2}; }
//+
Curve Loop(1) = {2};
//+
Curve Loop(2) = {1};
//+
Plane Surface(1) = {1, 2};
//+
Curve Loop(3) = {1};
//+
Plane Surface(2) = {3};
//+
Physical Surface("Gaine", 4) = {1};
//+
Physical Surface("Coeur", 5) = {2};
//+
Physical Curve("Interface", 6) = {1};
//+
Physical Curve("Bord", 7) = {2};
//+
Physical Point("centre", 8) = {1};
//+
Transfinite Curve {1} = 100 Using Progression 1;
//+
Transfinite Curve {2} = 120 Using Progression 1;
//+
Point{1} In Surface{2};
//+
Physical Point("Bord", 7) = {3};
//+
Physical Point("Interface", 6) = {2};