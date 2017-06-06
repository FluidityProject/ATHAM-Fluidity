//Box 1:
Point(1) = {0, 0, 0, 100};
Point(2) = {2000, 0, 0, 100};
Point(3) = {2000, 0, 2000, 100};
Point(4) = {0, 0, 2000, 100};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

out1[]=
Extrude {0,300,0} {
  Surface{6};
  Layers{3};
};

//For k In {0:5}
//    Printf("out1[%g] = %g",k,out1[k]);
//EndFor

//Box 2:
Point(101) = {2000, 0, 0, 100};
Point(102) = {4000, 0, 0, 100};
Point(103) = {4000, 0, 2000, 100};
Point(104) = {2000, 0, 2000, 100};

Line(101) = {101, 102};
Line(102) = {102, 103};
Line(103) = {103, 104};
Line(104) = {104, 101};

Line Loop(105) = {101, 102, 103, 104};
Plane Surface(106) = {105};

out2[]=
Extrude {0,300,0} {
  Surface{106};
  Layers{3};
};

//For k In {0:5}
//    Printf("out2[%g] = %g",k,out2[k]);
//EndFor

//Box 3:
Point(201) = {4000, 0, 0, 100};
Point(202) = {8000, 0, 0, 100};
Point(203) = {8000, 0, 2000, 100};
Point(204) = {4000, 0, 2000, 100};

Line(201) = {201, 202};
Line(202) = {202, 203};
Line(203) = {203, 204};
Line(204) = {204, 201};

Line Loop(205) = {201, 202, 203, 204};
Plane Surface(206) = {205};

out3[]=
Extrude {0,300,0} {
  Surface{206};
  Layers{3};
};

//For k In {0:5}
//    Printf("out3[%g] = %g",k,out3[k]);
//EndFor

//Domain:
Coherence;
Physical Surface(400) = {6}; //side front left
Physical Surface(401) = {out1[0]}; //side back left
Physical Surface(402) = {out1[2]}; //floor left
Physical Surface(403) = {out1[4]}; //roof left
Physical Surface(404) = {out1[5]}; //inlet
Physical Surface(405) = {106}; //side front middle
Physical Surface(406) = {out2[0]}; //side back middle
Physical Surface(407) = {out2[2]}; //floor middle
Physical Surface(408) = {out2[4]}; //roof middle
Physical Surface(409) = {206}; //side front right
Physical Surface(410) = {out3[0]}; //side back right
Physical Surface(411) = {out3[2]}; //floor right
Physical Surface(412) = {out3[3]}; //outlet
Physical Surface(413) = {out3[4]}; //roof right

Physical Volume(500) = {out1[1],out2[1],out3[1]};
