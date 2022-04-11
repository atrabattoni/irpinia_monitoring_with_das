lc_1 = 50;
lc_2 = 2.4;
lc_3 = 1.2;

Point(1) = {-150, 0, 0, lc_3};
Point(2) = {150, 0, 0, lc_3};
Point(3) = {-135, -7, 0, lc_3};
Point(4) = {135, -7, 0, lc_3};
Point(5) = {95, -21, 0, lc_2};
Point(6) = {-95, -21, 0, lc_2};
Point(7) = {-300, -0, 0, lc_1};
Point(8) = {300, -0, 0, lc_1};
Point(9) = {300, -200, 0, lc_1};
Point(10) = {-300, -200, 0, lc_1};
Point(13) = {-30, -25, -1.2, lc_2};
Point(14) = {30, -25, -0.2, lc_2};
Point(15) = {90, -10, -0.6, lc_3};
Point(16) = {-90, -10, -0.6, lc_3};

Line(1) = {8, 2};
Line(2) = {2, 4};
Line(3) = {4, 5};
Line(4) = {6, 3};
Line(5) = {3, 1};
Line(6) = {1, 7};
Line(7) = {7, 10};
Line(8) = {10, 9};
Line(9) = {9, 8};
Line(10) = {2, 1};
Line(13) = {5, 14};
Line(16) = {13, 6};
Line(17) = {4, 15};
Line(18) = {15, 16};
Line(19) = {16, 3};
Line(20) = {14, 13};

Curve Loop(1) = {7, 8, 9, 1, 2, 3, 13, 20, 16, 4, 5, 6};
Plane Surface(1) = {1};
Curve Loop(2) = {17, 18, 19, -4, -16, -20, -13, -3};
Plane Surface(2) = {2};
Curve Loop(3) = {10, -5, -19, -18, -17, -2};
Plane Surface(3) = {3};

Mesh.Algorithm = 6;
Mesh.ElementOrder = 2;
Mesh.RecombinationAlgorithm = 3;
Mesh.RecombineAll = 1;
Mesh.Smoothing = 100;
Mesh.MshFileVersion = 2.2;

Physical Surface("M1") = {1};
Physical Surface("M2") = {2};
Physical Surface("M3") = {3};
Physical Curve("Right") = {9};
Physical Curve("Top") = {1, 10, 6};
Physical Curve("Left") = {7};
Physical Curve("Bottom") = {8};
