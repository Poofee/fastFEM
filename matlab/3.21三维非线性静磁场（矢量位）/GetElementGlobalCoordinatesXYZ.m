function ElementNodeXYZ = GetElementGlobalCoordinatesXYZ(...
    eln,...
    ElementConnectivityTable,...
    GlobalNodeNumbering...
    )
%{
    Returns the global coordinates of the nodes associated with the
    argument element (eln)
    
    eln - Element Number
%}
A = ElementConnectivityTable(eln,:);
ElementNodeXYZ = GlobalNodeNumbering(A,:);


