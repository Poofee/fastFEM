function DisplayElements(...
    BrickElementConnectivityTable,...
    GlobalNodeNumbering,...
    TNEL...
    )
% Displays Nodes
DisplayNodes(GlobalNodeNumbering);

% Displays Elements
for eln = 1:TNEL
    ElementNodeXYZ = GetElementGlobalCoordinatesXYZ(...
        eln,...
        BrickElementConnectivityTable,...
        GlobalNodeNumbering);
    
    X = ElementNodeXYZ(:,1);    % Local node global x coordinates
    Y = ElementNodeXYZ(:,2);    % Local node global y coordinates
    Z = ElementNodeXYZ(:,3);    % Local node global z coordinates
    %{
    LineX1 = [1 2];  % Local nodes that correspond to an x line
    LineX2 = [3 4];
    LineX3 = [5 6];
    LineX4 = [7 8];
    LineY1 = [1 4];  % Local nodes that correspond to a y line
    LineY2 = [2 3];
    LineY3 = [5 8];
    LineY4 = [6 7];
    LineZ1 = [1 5];  % Local nodes that correspond to a z line
    LineZ2 = [2 6];
    LineZ3 = [3 7];
    LineZ4 = [4 8];
    line(X(LineX1),Y(LineX1),Z(LineX1),'Color',[0 0 0]);
    line(X(LineX2),Y(LineX2),Z(LineX2),'Color',[0 0 0]);
    line(X(LineX3),Y(LineX3),Z(LineX3),'Color',[0 0 0]);
    line(X(LineX4),Y(LineX4),Z(LineX4),'Color',[0 0 0]);
    line(X(LineY1),Y(LineY1),Z(LineY1),'Color',[0 0 0]);
    line(X(LineY2),Y(LineY2),Z(LineY2),'Color',[0 0 0]);
    line(X(LineY3),Y(LineY3),Z(LineY3),'Color',[0 0 0]);
    line(X(LineY4),Y(LineY4),Z(LineY4),'Color',[0 0 0]);
    line(X(LineZ1),Y(LineZ1),Z(LineZ1),'Color',[0 0 0]);
    line(X(LineZ2),Y(LineZ2),Z(LineZ2),'Color',[0 0 0]);
    line(X(LineZ3),Y(LineZ3),Z(LineZ3),'Color',[0 0 0]);
    line(X(LineZ4),Y(LineZ4),Z(LineZ4),'Color',[0 0 0]);
    %}
%     Lines = [1 3 5 7 1 2 5 6 1 2 3 4;...
%         2 4 6 8 4 3 8 7 5 6 7 8];
    line(X([1:end,1]),Y([1:end,1]),Z([1:end,1]),'Color',[0 0 0]);
%     text(...
%         mean(X),...             % Places text at the average x location
%         mean(Y),...             % Places text at the average y location
%         mean(Z),...             % Places text at the average z location
%         num2str(eln),...
%         'EdgeColor',[0 0 0],...
%         'FontWeight','bold');
    
    axis equal
end
