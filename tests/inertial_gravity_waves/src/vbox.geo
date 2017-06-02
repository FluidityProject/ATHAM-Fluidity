Point(1) = {0, 0, 0, 50};
Point(2) = {300000, 0, 0, 50};
Point(3) = {300000, 10000, 0, 50};
Point(4) = {0, 10000, 0, 50};

Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};

Physical Line(9) = {5};
Physical Line(10) = {6};
Physical Line(11) = {7};
Physical Line(12) = {8};

Transfinite Line{5, 7} = 601 Using Progression 1;
Transfinite Line{6, 8} = 21 Using Progression 1;

Line Loop(13) = {5, 6, 7, 8};
Plane Surface(14) = {13};
Transfinite Surface{14} Alternate;
Physical Surface(15) = {14};
