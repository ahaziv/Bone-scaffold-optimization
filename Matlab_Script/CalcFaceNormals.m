function [FaceNormals]=CalcFaceNormals(Faces, Vertices)
%% Summary
%Author: Itzik Ben Shabat
%Last Update: July 2014

%CalcFaceNormals recives a list of vrtexes and faces in FV structure
% and calculates the normal at each face and returns it as FaceNormals
%INPUT:
% Faces - an array containing each faces' vertices numbers
% Vertices - an array containing the verices coordinates
%OUTPUT:
%FaceNormals - an nX3 matrix (n = number of faces) containng the norml at each face
%% Code
% Get all edge vectors
e0 = Vertices(Faces(:,3),:)-Vertices(Faces(:,2),:);
e1 = Vertices(Faces(:,1),:)-Vertices(Faces(:,3),:);
% Calculate normal of face
FaceNormals = cross(e0,e1);
FaceNormals = normr(FaceNormals);
end

