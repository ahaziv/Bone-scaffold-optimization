function [Avertex] = CalcVertexAreas(Faces,Vertices)
%% Summary
% Author: Itzik Ben Shabat , Editor: Ziv Aharoni
% Last Update: march 2019

% summary: CalcVertexAreas calculates the voronoi areas at each vertex
% INPUT:
% Faces - Faces array containing vertex pointers
% Vertices - an array containing vertex coordinates

% OUTPUT:
% Avertex - [NvX1] voronoi area at each vertex

%% Code
% Get all edge vectors
e0 = Vertices(Faces(:,3),:)-Vertices(Faces(:,2),:);
e1 = Vertices(Faces(:,1),:)-Vertices(Faces(:,3),:);
e2 = Vertices(Faces(:,2),:)-Vertices(Faces(:,1),:);

%normalization procedure
%calculate face Area
%edge lengths
de0 = sqrt(e0(:,1).^2+e0(:,2).^2+e0(:,3).^2);
de1 = sqrt(e1(:,1).^2+e1(:,2).^2+e1(:,3).^2);
de2 = sqrt(e2(:,1).^2+e2(:,2).^2+e2(:,3).^2);
l2 = [de0.^2 de1.^2 de2.^2];

% using ew to calulate the cot of the angles for the voronoi area
% calculation. ew is the triangle barycenter, I later check if its inside or
% outide the triangle
ew = [l2(:,1).*(l2(:,2)+l2(:,3)-l2(:,1)) l2(:,2).*(l2(:,3)+l2(:,1)-l2(:,2)) l2(:,3).*(l2(:,1)+l2(:,2)-l2(:,3))];

s = (de0+de1+de2)/2;
% Af - face area vector
Af = sqrt(s.*(s-de0).*(s-de1).*(s-de2));%herons formula for triangle area, could have also used 0.5*norm(cross(e0,e1))

% calculate weights
Acorner = zeros(size(Faces,1),3);
Avertex = zeros(size(Vertices,1),1);

for i=1:size(Faces,1)
    % Calculate areas for weights according to Meyer et al. [2002]
    % check if the tringle is obtuse, right or acute
    
    if ew(i,1) <= 0
        Acorner(i,2) = -0.25*l2(i,3)*Af(i)/(e0(i,:)*e2(i,:)');
        Acorner(i,3) = -0.25*l2(i,2)*Af(i)/(e0(i,:)*e1(i,:)');
        Acorner(i,1) = Af(i)-Acorner(i,2)-Acorner(i,3);
    elseif ew(i,2) <= 0
        Acorner(i,3) = -0.25*l2(i,1)*Af(i)/(e1(i,:)*e0(i,:)');
        Acorner(i,1) = -0.25*l2(i,3)*Af(i)/(e1(i,:)*e2(i,:)');
        Acorner(i,2) = Af(i)-Acorner(i,1)-Acorner(i,3);
    elseif ew(i,3) <= 0
        Acorner(i,1) = -0.25*l2(i,2)*Af(i)/(e2(i,:)*e1(i,:)');
        Acorner(i,2) = -0.25*l2(i,1)*Af(i)/(e2(i,:)*e0(i,:)');
        Acorner(i,3) = Af(i)-Acorner(i,1)-Acorner(i,2);
    else
        ewscale = 0.5*Af(i)/(ew(i,1)+ew(i,2)+ew(i,3));
        Acorner(i,1) = ewscale*(ew(i,2)+ew(i,3));
        Acorner(i,2) = ewscale*(ew(i,1)+ew(i,3));
        Acorner(i,3) = ewscale*(ew(i,2)+ew(i,1));
    end
    Avertex(Faces(i,1)) = Avertex(Faces(i,1)) + Acorner(i,1);
    Avertex(Faces(i,2)) = Avertex(Faces(i,2)) + Acorner(i,2);
    Avertex(Faces(i,3)) = Avertex(Faces(i,3)) + Acorner(i,3);
end

end