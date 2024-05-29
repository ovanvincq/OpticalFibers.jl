// Gmsh project created on Tue Jul 26 15:40:16 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.1};
//+
Circle(1) = {0, 0, 0, 10, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 8, 0, 2*Pi};
//+
Pitch=2;
R=0.75;
//Ring1
x1={1,0.5,-0.5,-1,-0.5,0.5};
y1={0,0.866025403784439,0.866025403784439,0,-0.866025403784439,-0.866025403784439};
For i In {1:6}
	Circle(i+2) = {Pitch*x1[i-1], Pitch*y1[i-1], 0, 0.75, 0, 2*Pi};
EndFor
//Ring2
x2={2,1.5,1,0,-1,-1.5,-2,-1.5,-1,0,1,1.5};
y2={0,0.866025403784439,1.732050807568877,1.732050807568877,1.732050807568877,0.866025403784439,0,-0.866025403784439,-1.732050807568877,-1.732050807568877,-1.732050807568877,-0.866025403784439};
For i In {1:12}
	Circle(i+8) = {Pitch*x2[i-1], Pitch*y2[i-1], 0, 0.75, 0, 2*Pi};
EndFor	
//Ring3
x3={3,2.5,2,1.5,0.5,-0.5,-1.5,-2,-2.5,-3,-2.5,-2,-1.5,-0.5,0.5,1.5,2,2.5};
y3={0,0.866025403784439,1.732050807568877,2.598076211353316,2.598076211353316,2.598076211353316,2.598076211353316,1.732050807568877,0.866025403784439,0,-0.866025403784439,-1.732050807568877,-2.598076211353316,-2.598076211353316,-2.598076211353316,-2.598076211353316,-1.732050807568877,-0.866025403784439};
For i In {1:18}
	Circle(i+20) = {Pitch*x3[i-1], Pitch*y3[i-1], 0, 0.75, 0, 2*Pi};
EndFor	
//+
For i In {1:38}
	Curve Loop(i) = {i};
EndFor
//+
Plane Surface(1) = {2:38};
Physical Surface("Cladding", 1) = {1};
Transfinite Curve {1} = 160 Using Progression 1;
Plane Surface(2) = {1,2};
Physical Surface("PML", 2) = {2};
Transfinite Curve {2} = 200 Using Progression 1;
For i In {3:38}
	Plane Surface(i) = {i};
	Transfinite Curve {i} = 24 Using Progression 1;
EndFor
Physical Surface("Air", 3) = {3:38};
//+
Physical Curve("Edge", 1) = {1};
//+
Physical Point("Edge", 1) = {2};
//+
Physical Point("Center", 2) = {1};
//+
Point{1} In Surface{1};
//+
