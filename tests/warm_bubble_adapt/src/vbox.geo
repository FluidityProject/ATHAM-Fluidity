Point(1) = {0, 0, 0, 50};
Point(2) = {1000, 0, 0, 50};
Point(3) = {1000, 1500, 0, 50};
Point(4) = {0, 1500, 0, 50};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


Physical Line(8) = {1};
Physical Line(9) = {3};
Physical Line(10) = {4};
Physical Line(11) = {2};


Line Loop(5) = {1, 2, 3, 4};


Plane Surface(6) = {5};
Physical Surface(7) = {6};


Transfinite Line{1, 3} = 200 Using Progression 1;
Transfinite Line{2, 4} = 301 Using Progression 1;

Transfinite Surface{6} Alternate;
