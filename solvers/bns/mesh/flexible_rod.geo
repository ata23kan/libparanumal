lc_box = 0.25;    // Mesh size far field
lc_wall = 0.025;  // Mesh size near structure

// Geometry Dimensions
D = 1;   // Square side length
L = 4;   // Rod length
H = 0.2;  // Rod thickness

// Farfield box limits
xmin = -5*D;
xmax =  10*D;
ymin = -5*D;
ymax =  5*D;

// PML width
w = 2*D;

// Square
Point(1) = {-D/2, -D/2, 0, lc_wall};
Point(2) = { D/2, -D/2, 0, lc_wall};
Point(3) = { D/2,  D/2, 0, lc_wall};
Point(4) = {-D/2,  D/2, 0, lc_wall};

// Rod
Point(5) = {D/2, -H/2, 0, lc_wall};
Point(6) = {D/2 + L, -H/2, 0, lc_wall};
Point(7) = {D/2 + L,  H/2, 0, lc_wall};
Point(8) = {D/2,  H/2, 0, lc_wall};

Line(1) = {1, 2};  // Square bottom
Line(2) = {2, 5};  // Square-rod intersect bottom
Line(3) = {5, 6};  // Rod bottom
Line(4) = {6, 7};  // Rod right
Line(5) = {7, 8};  // Rod top
Line(6) = {8, 3};  // Square-rod intersect top
Line(7) = {3, 4};  // Square top
Line(8) = {4, 1};  // Square left

// Bounding Box
Point(9)  = {xmin, ymin, 0, lc_box};
Point(10) = {xmax, ymin, 0, lc_box};
Point(11) = {xmax, ymax, 0, lc_box};
Point(12) = {xmin, ymax, 0, lc_box};

Line(9)  = {9, 10};  // Domain bottom
Line(10) = {10,11};  // Domain right
Line(11) = {11,12};  // Domain right
Line(12) = {12, 9};  // Domain right

Curve Loop (1) = {9, 10, 11, 12};          // Fluid domain
Curve Loop (2) = {1, 2, 3, 4, 5, 6, 7, 8}; // Object
Plane Surface(3) = {1, 2};

// PML Points
// upper left
Point(13) = {xmin-w, ymax+w, 0};
Point(14) = {xmin-w, ymax, 0};
Point(15) = {xmin, ymax+w, 0};
//lower left
Point(16) = {xmin-w, ymin-w, 0};
Point(17) = {xmin-w, ymin, 0};
Point(18) = {xmin, ymin-w, 0};
//upper right
Point(19) = {xmax+w, ymax+w, 0};
Point(20) = {xmax+w, ymax, 0};
Point(21) = {xmax, ymax+w, 0};
//lower right
Point(22) = {xmax+w, ymin-w, 0};
Point(23) = {xmax+w, ymin, 0};
Point(24) = {xmax, ymin-w, 0};

MeshSize{13:24} = lc_box;

// Upper Left PML
Line(13) = {13, 14};
Line(14) = {14, 12};
Line(15) = {12, 15};
Line(16) = {15, 13};
Curve Loop (3) = {13, 14, 15, 16};
Plane Surface (4) = {3};

// Lower left PML
Line(17) = {16, 17};
Line(18) = {17, 9};
Line(19) = {9, 18};
Line(20) = {18, 16};
Curve Loop (4) = {17, 18, 19, 20};
Plane Surface (5) = {4};

// Upper right PML
Line(21) = {19, 20};
Line(22) = {20, 11};
Line(23) = {11, 21};
Line(24) = {21, 19};
Curve Loop (5) = {21, 22, 23, 24};
Plane Surface (6) = {5};

// Lower right PML
Line(25) = {22, 23};
Line(26) = {23, 10};
Line(27) = {10, 24};
Line(28) = {24, 22};
Curve Loop (6) = {25, 26, 27, 28};
Plane Surface (7) = {6};

// Left X pml
Line(29) = {14, 17};
Curve Loop(7) = {29, 18, -12, -14};
Plane Surface(8) = {7};

// Upper Y pml
Line(30) = {15, 21};
Curve Loop(8) = {11, 15, 30, -23};
Plane Surface(9) = {8};

// Right X pml
Line(31) = {20, 23};
Curve Loop(9) = {10, -22, 31, 26};
Plane Surface(10) = {9};

// Lower Y pml
Line(32) = {24, 18};
Curve Loop(10) = {9, 27, 32, -19};
Plane Surface(11) = {10};


// Physical Groups
Physical Line("Wall", 1)    = {3, 4, 5};
//Physical Line("Wall", 1)    = {1, 2, 3, 4, 5, 6, 7, 8};
Physical Line("Square", 10) = {2, 6, 7, 8, 1};
Physical Line("Outflow", 3) = {13, 16, 30, 24, 21, 31, 25, 28, 32, 20, 17, 29};

Physical Surface("Interior", 9)  = {3};
Physical Surface("XPML", 100)    = {8, 10};
Physical Surface("YPML", 200)    = {9, 11};
Physical Surface("XYPML", 300)   = {4, 5, 6, 7};
