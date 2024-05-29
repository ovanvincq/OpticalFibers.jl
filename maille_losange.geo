// Gmsh project created on Tue Jul 26 15:40:16 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.01};
//+
Circle(1) = {0, 0, 0, 0.25, 0, 2*Pi};
//+
Point(3) = {0.5, -1/2/Sqrt(3),0};
Point(4) = {0.5, 1/2/Sqrt(3), 0};
Point(5) = {0, 1/Sqrt(3), 0};
Point(6) = {-0.5, 1/2/Sqrt(3), 0};
Point(7) = {-0.5, -1/2/Sqrt(3), 0};
Point(8) = {0 , -1/Sqrt(3), 0};
Line(2) = {3,4};              //Lines
Line(3) = {4,5};
Line(4) = {5,6};
Line(5) = {6,7};
Line(6) = {7,8};
Line(7) = {8,3};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2,3,4,5,6,7};
//+
Plane Surface(1) = {1};
//+
Plane Surface(2) = {1,2};
//+
Physical Surface("Inclusion", 5) = {1};
//+
Physical Surface("Gaine", 6) = {2};
//+
Physical Curve("Interface", 7) = {1};
//+
Physical Point("Interface", 7) = {2};
//+
Physical Curve("Bord", 8) = {2,3,4};
//+
Physical Point("Bord", 8) = {3,4,5,6,7,8};
//+
Transfinite Curve {1} = 64 Using Progression 1;
//+
Transfinite Curve {2,3,4,5,6,7} = 25 Using Progression 1;
//+
Point{1} In Surface{1};
//+
Periodic Curve {2} = {5} Translate {1,0,0};
//+
Periodic Curve {3} = {6} Translate {0.5,Sqrt(3)/2,0};
//+
Periodic Curve {4} = {7} Translate {-0.5,Sqrt(3)/2,0};