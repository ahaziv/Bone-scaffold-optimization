function [Faces, Vertices, NormRatio, Flag] = CreateShape(ShrinkFactor)
% this function imports a .dat file of a certain mesh and returns:
% Faces - the element array of vertices triangles
% Vertices - the array of node coordinates
% NormRatio - a normalizing ratio used to speed calculations
% Flag - a signal telling whether the input mesh is of acceptable form
% (Flag = 0 usually means a corrupt geometry)
addpath(fullfile(cd, '..'))

Flag = 1;
InitData = readtable('SYS.dat','Delimiter','\n','ReadVariableNames',false); 
idx = ismember(InitData.Var1, '-1');
[idxnum,~] = find(idx);

opts = delimitedTextImportOptions(...
        'Delimiter', ' ',...
        'ConsecutiveDelimitersRule', 'join','LeadingDelimitersRule','ignore',...
        'MissingRule', 'omitvar','LineEnding','\n','ImportErrorRule','omitvar');

InitData = readtable('SYS.dat',opts);
% in the cases where there are multipile element tables (idxnum>2)
% or volumetric elements (width(InitData)>16) send failure signal.
if length(idxnum) > 2 || width(InitData) > 16
    Flag = 0;
end
% extracting the Elements and Node tables
NoTable = InitData(26:idxnum(1)-1,2:4);
ElTable = InitData(idxnum(1)+7:idxnum(2)-1,12:15);
%% converting the tables to arrays
Nodes = table2array(NoTable);
Vertices = str2double(Nodes).*ShrinkFactor;
Elements = table2array(ElTable);
Faces = str2double(Elements);
% in case of quadrilateral elements send failure signal
if ~isempty(find(Faces(:,3)~=Faces(:,4),1))
    Flag = 0;
end
Faces = Faces(:,1:3);
%% normalizing
Vertices = Vertices - mean(Vertices); % translating the array so surround the origin
NormRatio = 1/(max(max(abs(Vertices))));
Vertices = Vertices.*NormRatio;