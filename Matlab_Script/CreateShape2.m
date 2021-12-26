function [Faces, Vertices, NormRatio, Flag] = CreateShape2(ShrinkFactor)
Faces = xlsread('Bartolo_Elements_Data','C:E');
Vertices = xlsread('Bartolo_Nodes_Data','B:D');
Vertices = Vertices.*ShrinkFactor;
[Vertices, NormRatio] = NormPoints(Vertices);
Flag = 1;
end

function [Vertices, NormRatio] = NormPoints(Vertices)
% this function normalizes a vetex array to a unit sphere
% NormRatio - the normalization ratio (
Vertices = Vertices - mean(Vertices); % translating the array so surround the origin
NormRatio = 1/(max(max(abs(Vertices))));
Vertices = Vertices.*NormRatio;
end
