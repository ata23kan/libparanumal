//cl__1 = 0.05;
//Point(1) = {0, 0, 0, cl__1};
//Point(2) = {1, 0, 0, cl__1};
//Point(3) = {1, 1, 0, cl__1};
//Point(4) = {0, 1, 0, cl__1};
//Line(1) = {1, 2};
//Line(2) = {2, 3};
//Line(3) = {3, 4};
//Line(4) = {4, 1};
//Line Loop(6) = {4, 1, 2, 3};
//Plane Surface(6) = {6};
//Physical Surface("Domain", 9) = {6};
//Physical Line("Dirichlet", 1) = {1, 2, 3, 4};
//Physical Line("Stationary", 1) = {1, 2, 3, 4};
//Physical Line("Moving", 2) = {3};


cl__1 = 0.1;

// Bottom face (z=0)
Point(1) = {0, 0, 0, cl__1};
Point(2) = {1, 0, 0, cl__1};
Point(3) = {1, 1, 0, cl__1};
Point(4) = {0, 1, 0, cl__1};

// Top face (z=1)
Point(5) = {0, 0, 1, cl__1};
Point(6) = {1, 0, 1, cl__1};
Point(7) = {1, 1, 1, cl__1};
Point(8) = {0, 1, 1, cl__1};

// Bottom edges
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Top edges
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Vertical edges
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// 6 faces
Line Loop(101) = {1, 2, 3, 4}; Plane Surface(101) = {101}; // z=0
Line Loop(102) = {5, 6, 7, 8}; Plane Surface(102) = {102}; // z=1
Line Loop(103) = {1, 10, -5, -9}; Plane Surface(103) = {103}; // y=0
Line Loop(104) = {2, 11, -6, -10}; Plane Surface(104) = {104}; // x=1
Line Loop(105) = {3, 12, -7, -11}; Plane Surface(105) = {105}; // y=1
Line Loop(106) = {4, 9, -8, -12}; Plane Surface(106) = {106}; // x=0

Surface Loop(200) = {101, 102, 103, 104, 105, 106};
Volume(300) = {200};

// Physical groups
Physical Volume("Domain", 9) = {300};
Physical Surface("Dirichlet", 1) = {101, 102, 103, 104, 105, 106};