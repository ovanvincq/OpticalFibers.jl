// Gmsh project created on Tue Jul 26 15:40:16 2022
SetFactory("OpenCASCADE");
//+
Pitch=6;
//+
Dnorm=0.5;
//+
Point(1) = {0, 0, 0, 0.01};
//+
Circle(1) = {0, 0, 0, Dnorm/2*Pitch, 0, 2*Pi};
//+
Point(4) = {-0.5*Pitch, -0.5*Pitch, 0};
Point(5) = {-0.5*Pitch, +0.5*Pitch, 0};
Point(6) = {0.5*Pitch, 0.5*Pitch, 0};
Point(7) = {0.5*Pitch, -0.5*Pitch, 0};
Line(2) = {4,5};              //Lines
Line(3) = {5,6};
Line(4) = {6,7};
Line(5) = {7,4};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2,3,4,5};
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
Physical Curve("Bord", 8) = {3,4};
//+
Physical Point("Bord", 8) = {4,5,6,7};
//+
Transfinite Curve {1} = 64 Using Progression 1;
//+
Transfinite Curve {2,3,4,5} = 25 Using Progression 1;
//+
Point{1} In Surface{1};
//+
Periodic Curve {2} = {4} Translate {-Pitch,0,0};
//+
Periodic Curve {5} = {3} Translate {0,-Pitch,0};
