function [VertexNormals] = BoundaryNormal(Vertices,VertexNormals)
% this function adjusts the normal movement direction along the boundaries.
% it does so by first detecting the boundaries and then projecting the
% surface normals at each boundary vertex upon the appropriate boundary.
for ii=1:3
     VertexNormals(Vertices(:,ii+3)~=0,ii) = 0;
end
VertexNormals(sum(Vertices(:,:),2)~=0,:) = VertexNormals(sum(Vertices(:,:),2)~=0,:)./sqrt(sum(VertexNormals(sum(Vertices(:,:),2)~=0,:).^2,2));
for ii=1:3
     VertexNormals(isnan(VertexNormals(:,ii)),ii) = 0;
end