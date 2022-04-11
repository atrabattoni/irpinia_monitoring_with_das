lc_1 = 50;

Point(7) = {-300, -0, 0, lc_1};
Point(8) = {300, -0, 0, lc_1};
Point(9) = {300, -200, 0, lc_1};
Point(10) = {-300, -200, 0, lc_1};

Line(6) = {8, 7};
Line(7) = {7, 10};
Line(8) = {10, 9};
Line(9) = {9, 8};

Curve Loop(1) = {6, 7, 8, 9};
Plane Surface(1) = {1};
Transfinite Surface {1};

Physical Surface("M1") = {1};
Physical Curve("Right") = {9};
Physical Curve("Top") = {6};
Physical Curve("Left") = {7};
Physical Curve("Bottom") = {8};

Mesh.Algorithm = 6;
Mesh.ElementOrder = 2;
Mesh.RecombinationAlgorithm = 3;
Mesh.RecombineAll = 1;
Mesh.Smoothing = 100;
Mesh.MshFileVersion = 2.2;
