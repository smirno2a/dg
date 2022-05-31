lcx = 0.3;
Point(1) = {-3., -3.,0,lcx};
Point(2) = {3.,-3.,0,lcx};
Point(3) = {3.,3.,0,lcx};
Point(4) = {-3.,3.,0,lcx};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(8) = {2,3,4,1};
Plane Surface(9) = {8};

Physical Surface (0) = {9};
Physical Line (30000) = {1,2,3,4};

