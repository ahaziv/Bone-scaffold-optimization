function [Vertices] = AdjustBorders(Vertices, BorderVal)
% this function detects vertices moving from one boundary to the other,
% changing the Vertex boundary accordingly.
% BorderVal is an array containing the minimum and maximum values of the
% geometry along the three axes:
% BorderVal(1)= minimun on x.
% BorderVal(2)= maximum on x.
% BorderVal(3)= minimun on y axis. etc...
% eps - a certain distance in which the vertex is manually moved to the boundary
BorderVertices = Vertices(sum(Vertices(:,4:6),2)~=0,:);
for ii = 1:3
    BorderVertices(BorderVertices(:,ii) < BorderVal(2*ii-1),ii) = 1;
    BorderVertices(BorderVertices(:,ii) < BorderVal(2*ii-1),ii) = BorderVal(2*ii-1);
    BorderVertices(BorderVertices(:,ii) > BorderVal(2*ii),ii) = 1;
    BorderVertices(BorderVertices(:,ii) > BorderVal(2*ii),ii) = BorderVal(2*ii);
end
Vertices(sum(Vertices(:,4:6),2)~=0,:) = BorderVertices;