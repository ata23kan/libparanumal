cl__1 = 0.1;
Point(1) = {0, 0, 0, cl__1};
Point(2) = {5, 0, 0, cl__1};
Point(3) = {5, 2, 0, cl__1};
Point(4) = {0, 2, 0, cl__1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {4, 1, 2, 3};
Plane Surface(6) = {6};
Transfinite Surface {6};
/Physical Line("Wall",1) = {1, 3, 4};
/Physical Line("Outflow",2) = {2};

Physical Line("Walls",1) = {1, 4};
Physical Line("Outflow",2) = {2};
Physical Line("Wallm",3) = {3};
Physical Surface("Domain") = {6};
