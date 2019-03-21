SetFactory("OpenCASCADE");

Mesh.CharacteristicLengthMin = 3e-3;
Mesh.CharacteristicLengthMax = 3e-3;
Geometry.OCCTargetUnit = "M";

DefineConstant[
  angle = {90, Name "Parameters/wedge angle"}
  extrude = {0.01, Name "Parameters/extrusion length (with mesh)"}
];

a() = ShapeFromFile("component8.step");

Cylinder(2) = {0,0.15,0, 0,0.2,0, 0.04, angle*2*Pi/360};

BooleanIntersection{ Volume{a(0)}; Delete; }{ Volume{2}; Delete; }

Extrude {0, extrude, 0} {
  Surface{1}; Layers{5}; Recombine;
}
