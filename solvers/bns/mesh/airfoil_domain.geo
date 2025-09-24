
Include "airfoil.geo";
// Farfield box limits
xmin = -5.0;
xmax =  10.0;
ymin = -5.0;
ymax =  5.0;

// PML width
w = 2.0;

// --- Farfield rectangle ---
Point(132) = {xmin, ymin, 0};
Point(133) = {xmax, ymin, 0};
Point(134) = {xmax, ymax, 0};
Point(135) = {xmin, ymax, 0};

// PML Points
// upper left
Point(136) = {xmin-w, ymax+w, 0};
Point(137) = {xmin-w, ymax, 0};
Point(138) = {xmin, ymax+w, 0};
//lower left
Point(139) = {xmin-w, ymin-w, 0};
Point(140) = {xmin-w, ymin, 0};
Point(141) = {xmin, ymin-w, 0};
//upper right
Point(142) = {xmax+w, ymax+w, 0};
Point(143) = {xmax+w, ymax, 0};
Point(144) = {xmax, ymax+w, 0};
//lower right
Point(145) = {xmax+w, ymin-w, 0};
Point(146) = {xmax+w, ymin, 0};
Point(147) = {xmax, ymin-w, 0};

// Physical box
Line(3) = {132, 133};
Line(4) = {133, 134};
Line(5) = {134, 135};
Line(6) = {135, 132};

// PML region
Line(7) = {136, 137};
Line(8) = {136, 138};
Line(9) = {137, 135};
Line(10) = {138, 135};
Line(11) = {144, 142};
Line(12) = {142, 143};
Line(13) = {143, 134};
Line(14) = {134, 144};
Line(15) = {133, 146};
Line(16) = {146, 145};
Line(17) = {145, 147};
Line(18) = {147, 133};
Line(19) = {132, 141};
Line(20) = {141, 139};
Line(21) = {139, 140};
Line(22) = {140, 132};
Line(23) = {140, 137};
Line(24) = {138, 144};
Line(25) = {143, 146};
Line(26) = {147, 141};


Curve Loop(1) = {3, 4, 5, 6};	// Bounding Box
Curve Loop(2) = {1, 2};				// Airfoil
Plane Surface(3) = {1, 2};

// Upper left pml
Curve Loop(3) = {7, 9, -10, -8};
Plane Surface(4) = {3}; 
// Lower left pml
Curve Loop(4) = {19, 20 , 21, 22};
Plane Surface(5) = {4}; 
// Upper right pml
Curve Loop(5) = {11, 12, 13, 14};
Plane Surface(6) = {5}; 
// Lower right pml
Curve Loop(6) = {15, 16, 17, 18};
Plane Surface(7) = {6};

// Left X pml
Curve Loop(7) = {23, 9, 6, -22};
Plane Surface(8) = {7};
// Right X pml
Curve Loop(8) = {4, -13, 25, -15};
Plane Surface(9) = {8};

// Upper Y pml
Curve Loop(9) = {10, -5, 14, -24};
Plane Surface(10) = {9};
// Lower Y pml
Curve Loop(10) = {19, -26, 18, -3};
Plane Surface(11) = {10};


MeshSize{1:131} = 0.005;   // airfoil points (small size)
MeshSize{132, 133, 134, 135} = 0.25; // farfield corners (coarser)
MeshSize{136:147} = 0.25; // PML corners (coarser)

// Physical groups
Physical Curve("Wall", 1)       = {1, 2};
Physical Curve("Outflow", 3)  	= {7,8,24,11,23,21,20,26,17,12,16,25};
//Physical Curve("Outflow", 3) 	= {12, 16, 25};
Physical Surface("Interior", 9) = {3};
Physical Surface("XPML", 100)	  = {8, 9};
Physical Surface("YPML", 200)	  = {10, 11};
Physical Surface("XYPML", 300)  = {4, 5, 6, 7};