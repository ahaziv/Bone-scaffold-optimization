function [RefFaces, RefVertices] = ReflectBoundaryMesh(Faces, Vertices, PlaneParams, eps)
% PlaneParams(1) - the boundary type (XY,XZ,YZ)
% PlaneParams(2) - the boundary value (Xval, Yval, Zval)
% eps - the maximal distance from the boundary value to detect boundary vertices 

if length(PlaneParams) < 4
    % finding the appropriate boundary vertices 
    [BoundVerIndex,~] = find(Vertices(:,PlaneParams(1)+3));
%     BoundVerIndex(Vertices(BoundVerIndex,PlaneParams(1))-PlaneParams(2)~=0) = [];
    BoundVerIndex(abs(Vertices(BoundVerIndex,PlaneParams(1))-PlaneParams(2))>eps) = [];
    % finding the appropriate Faces to reflect and the the indeces to reflect
    [BoundFacesIdx,~] = find(sum(ismember(Faces,BoundVerIndex),2));
    RefFaces = Faces(BoundFacesIdx,:);
    RefVerIndex = RefFaces(~ismember(RefFaces,BoundVerIndex));
    RefVerIndex = unique(RefVerIndex);
    RefVertices = Vertices(RefVerIndex,:);
    % reflecting the vertices
    RefVertices(:,PlaneParams(1)) = RefVertices(:,PlaneParams(1))...
        - 2*(RefVertices(:,PlaneParams(1)) - PlaneParams(2));
    RefVertices(:,4:7) = 0;
    % building the Faces array
    for ii = 1:length(RefVerIndex)
        RefFaces(RefFaces == RefVerIndex(ii)) = length(Vertices(:,1)) + ii;
    end
    % flipping the 2nd and third columns of RefFaces so inexing will be aligned to the original array 
    temp = RefFaces(:,3);
    RefFaces(:,3) = RefFaces(:,2);
    RefFaces(:,2) = temp;

else
    % Case 2
    % PlaneParams - a 4 parameter set defining the plane as followes:
    % Ax+By+Cz-D=0. PlaneParams =[A B C D]

    % distance is calculated through: (A*x1+B*x2+C*x3+D)/(A^2+B^2+c^2)^0.5
    % where the sign of the value represent the direction to the plane
    % (positive is...
    Distance = (sum(PlaneParams(1:3).*Vertices(1:3))+PlaneParams(4))...
        /(sum(PlaneParams(1:3).^2)^(1/2));

    % The normal to the plane is calculated by: Normal = (A,B,C)/(A^2+B^2+c^2)^0.5
    Normal = PlaneParams(1:3)/(sum(PlaneParams(1:3).^2)^(1/2));
end
