// Gmsh project created on Tue Jul 26 15:40:16 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.1};
//+
Ellipse(1) = {0, 0, 0, 4, 2, 0, 2*Pi};
//+
Ellipse(2) = {0, 0, 0, 8, 4, 0, 2*Pi};
//+
Point(4) = {-12, -8, 0};
Point(5) = {-12, +8,0};
Point(6) = {12, 8, 0};
Point(7) = {12, -8, 0};
Line(11) = {4,5};              //Lines
Line(12) = {5,6};
Line(13) = {6,7};
Line(14) = {7,4};
Point(8) = {-15, -11, 0};
Point(9) = {-15, +11, 0};
Point(10) = {15, 11, 0};
Point(11) = {15, -11, 0};
Line(15) = {8,9};              //Lines
Line(16) = {9,10};
Line(17) = {10,11};
Line(18) = {11,8};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2};
//+
Curve Loop(3) = {11,12,13,14};
//+
Curve Loop(4) = {15,16,17,18};
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
Physical Curve("Bord", 7) = {15,16,17,18};
//+
Physical Point("Bord", 7) = {8,9,10,11};
//+
Physical Point("centre", 8) = {1};
//+
Transfinite Curve {1} = 64 Using Progression 1;
//+
Transfinite Curve {2} = 100 Using Progression 1;
//+
Transfinite Curve {11,12,13,14} = 25 Using Progression 1;
//+
Transfinite Curve {15,16,17,18} = 60 Using Progression 1;
//+
Point{1} In Surface{1};
//+

