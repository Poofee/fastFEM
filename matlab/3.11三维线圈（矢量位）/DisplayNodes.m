function DisplayNodes(GlobalNodeNumbering)
x = GlobalNodeNumbering(:,1);
y = GlobalNodeNumbering(:,2);
z = GlobalNodeNumbering(:,3);
% plot3(x,y,z,'.r','MarkerSize',14)
hold on;

% for node = 1:length(GlobalNodeNumbering(:,1))
%     text(...
%         GlobalNodeNumbering(node,1),...
%         GlobalNodeNumbering(node,2),...
%         GlobalNodeNumbering(node,3),...
%         num2str(node),...
%         'BackgroundColor',[0 1 1],...
%         'FontWeight','bold');
% end
plot3(x,y,z,'.r','MarkerSize',14)
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
title('Nodes');