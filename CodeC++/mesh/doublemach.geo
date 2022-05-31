lc = 0.08;
Point(1) = {0,1.,0,lc};
Point(2) = {0,0,0,lc};
Point(3) = {0.166666666666666666,0,0,lc};
Point(4) = {3.5,0,0,lc};
Point(5) = {3.5,1,0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,1};

Line Loop(8) = {2,3,4,5,1};
Plane Surface(9) = {8};

Physical Surface (0) = {9};
Physical Line (40000) = {3};
Physical Line (2) = {4}; 
Physical Line (30000) = {5,1,2};
