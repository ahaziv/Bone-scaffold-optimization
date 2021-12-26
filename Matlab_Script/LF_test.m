%% This is the main code for simulating CCTG 
clear all; close all; clc;
FilePath = mfilename('fullpath');
idcs = strfind(FilePath,'\');
ParDir = FilePath(1:idcs(end-1)-1);
addpath(FilePath);
addpath(ParDir);

% addpath(fullfile(cd, '..'))
% addpath(fullfile(cd))

%% input parameters
Radius = 325*10^-6;
VUnitCell = 8*Radius^3;     % the total volume of the unit cell

ShrinkFactor = 10^-3;
aPar = 8.534*10^-8;                       % the true parameter is 7*10^-8 but was normalized
bPar = 50*10^-6;                % the true parameter is 50*10^-6 but was normalized

%% general parameters
MaxNumIter = 500;              % maximal number of iteration steps
TimeLimit = 1;
LFactor = 0.24;             % Lambda multiplication factor
DeltaTheta = 0.1;          % the square of the minimal angle for remeshing 
FlagBoundary = 0;           % this flag tells the program wether to impose boundary conditions
FlagAnaSol = 0;
AnaSolution = 1;
if FlagAnaSol
    AnalyticParam = struct('GenEps', 0.5, 'SphereEps', 0.2, 'CylinderEps', 0.5);
end

Faces = cell(2,1);
Vertices = cell(2,1);
CurveCalcMethod = 'Average';

GrwothRateBreak = 0.1;          % a growth rate factor which tells the algorithm when to stop
Vgrowth = zeros(MaxNumIter,1); % this array stores the filled up volume at each iteration
VGrowthRate = zeros(MaxNumIter,1);
Lambda = zeros(MaxNumIter,1);
TotTime = 0;                % the total time of the simulation

% VVacantTot = real(CalcVolume(FW, FD, ST));        % the total vacant volume
%% Generate the example
[Faces{1}, Vertices{1}, NormRatio, GoodGeom] = CreateShape(ShrinkFactor);  
temp = Faces{1}(:,2);
Faces{1}(:,2) = Faces{1}(:,3);
Faces{1}(:,3) = temp;
Vertices{1}(1,4:7) = 0;
[emag] = ecalc(Faces{1},Vertices{1}(:,1:3)); 
DeltaRemesh = min(min(emag));
DeltaRemeshMax = 1.3*max(max(emag));
maxis = 1.3*max(max(abs(Vertices{1}(:,1:3))));

%% detect mesh edges for later use
if FlagBoundary
    [Vertices{1}, BoundaryVal, BoundStrVer] = DetectBoundaries(Faces{1},Vertices{1}, 0.5*DeltaRemesh);
    if length(BoundStrVer) ~= 4
        GoodGeom = 0;
    end
end

% in case the geometry is corrupted assign 10^6 to the Total loss and abort
if ~GoodGeom
    CurveLoss = 0; PoreSizeLoss = 0; TotLoss = 10^6;
    fileID = fopen(strcat(fullfile(cd, '..'),'\Loss.txt'), 'w+');
%     fileID = fopen('D:\Documents\2nd degree\Mod Frontier analysis\Curvature Model\Loss.txt', 'w+');
    fprintf(fileID,'%4.4f\n',CurveLoss);
    fprintf(fileID,'%4.4f\n',PoreSizeLoss);
    fprintf(fileID,'%4.4f\n',TotLoss);
    fclose(fileID);
    return
end

%% open a figure
fig_h = figure('name','Scaffold','numbertitle','off','color',[0.75 0.75 0.75]);
ax = axes('DataAspectRatio', [1,1,1]);
set(gcf, 'OuterPosition', get(0, 'Screensize'));
axis tight
%% run the algorithm
TotVgrowth = 0;
ii = 1;
try
    while (TotTime < TimeLimit) && (ii < MaxNumIter)% && (TotVgrowth < MaxVolFillRate * Vtot)
        if FlagBoundary % reflecting the boundaries to calculate the curvatures and normals
            OrigFacesLength = length(Faces{1}(:,1));
            OrigVertiLength = length(Vertices{1}(:,1));
            for ll = 1:6
                [RefFaces, RefVertices] = ReflectBoundaryMesh(Faces{1}, Vertices{1}, [ceil(0.5*ll) BoundaryVal(ll)], 0.5*DeltaRemesh);
                Faces{1} = [Faces{1}; RefFaces];
                Vertices{1} = [Vertices{1}; RefVertices];
            end
        end
        % calculating the curvatures and normals
        [FaceNormals] = CalcFaceNormals(Faces{1}(:,1:3),Vertices{1}(:,1:3));
        [VertexNormals, Avertex, Acorner,up,vp] = CalcVertexNormals(Faces{1}(:,1:3),Vertices{1}(:,1:3),FaceNormals);
        [VertexSFM, wfp] = CalcCurvature(Faces{1}(:,1:3),Vertices{1}(:,1:3),VertexNormals,FaceNormals,Avertex,Acorner,up,vp);

        if FlagBoundary % after calculating the curvatures and normals, cut off the reflected entities
            Faces{1} = Faces{1}(1:OrigFacesLength,:);
            Vertices{1} = Vertices{1}(1:OrigVertiLength,:);
            VertexNormals = VertexNormals(1:OrigVertiLength,:);
            VertexSFM = VertexSFM(:,:,1:OrigVertiLength);
            Avertex = Avertex(1:OrigVertiLength);
        end

        %% Calculating the main curvature in a point and making a step in the direction
        PartVertices = [];
        if FlagAnaSol
            [CurvatureVec, ColorVec, AnaSolution] = CurveMethod(CurveCalcMethod, VertexSFM, AnalyticParam, PartVertices);
            if AnaSolution
                break
            end
        else
            [CurvatureVec, ColorVec, AnaSolution] = CurveMethod(CurveCalcMethod, VertexSFM);
        end
        if ii == 1
            minmaxColor = [min(ColorVec) max(abs(ColorVec))];
        end

        Lambda(ii) = -LFactor*DeltaRemesh/min(min(CurvatureVec),-max(CurvatureVec)); % setting a value of lambda according to the maximal avg curvature
        Step = Lambda(ii)*VertexNormals.*CurvatureVec;
        Vertices{2} = Vertices{1} - [Step, zeros(size(Vertices{1},1),4)];
        %% Calculating the total volume filled up in this iteration
        AvertexNew = CalcVertexAreas(Faces{1}(:,1:3),Vertices{2}(:,1:3));
        Vgrowth(ii) = sum(((Avertex+AvertexNew)/2).*sqrt(sum(Step.^2,2)));
        TotVgrowth = TotVgrowth + Vgrowth(ii)/(NormRatio^3);
        VGrowthRate(ii) = Vgrowth(ii)/(NormRatio*Lambda(ii));
        if VGrowthRate(ii) < GrwothRateBreak*VGrowthRate(1)
            break
        end
        TotTime = TotTime + Lambda(ii);
        %% Remeshing the shape to avoid clustering of vertices
        if FlagBoundary
            [Faces{2}, Vertices{2}, BoundStrVer] = Remesh(Faces{1}, Vertices{2}, DeltaRemesh, DeltaTheta, DeltaRemeshMax, BoundStrVer);
            % imposing special boundary conditions
            [Vertices{2}, BoundStrVer] = SpecialBoundaryCond(Vertices{2} ,BoundStrVer);
        else
            [Faces{2}, Vertices{2}] = Remesh(Faces{1}, Vertices{2}, DeltaRemesh, DeltaTheta, DeltaRemeshMax);
        end    
        %% Plotting the mesh  
        PlotMesh(Faces{1}, Vertices{1}, ColorVec, maxis, minmaxColor, FlagBoundary, [], [])
        drawnow

        ii = ii+1;
        Faces{1} = Faces{2};
        Vertices{1} = Vertices{2};
    end
catch
    TotVgrowth = sqrt(-1);
end
%% Reorganizing the output
if AnaSolution
    if AnaSolution(1) == 1
        Vvacant = 4*pi*AnaSolution(2)^3/3;
        Lambda(ii) = (AnaSolution(2)^2)/4;
    elseif AnaSolution(1) == 2
        % in the case of a cylinder the length of the cylinder is
        % estimated using shapes' surface
        SurfArea = sum(AvertexNew);
        Vvacant = SurfArea*AnaSolution(2)/2;
        Lambda(ii) = (AnaSolution(2)^2)/2;
    end
    Vgrowth(ii) = Vvacant;
    Vgrowth = Vgrowth(1:ii);
    Lambda = Lambda(1:ii);
else
    Vgrowth = Vgrowth(1:ii-1);
    Lambda = Lambda(1:ii-1);
end
%% calculating the curvature loss function
Vgrowth = Vgrowth/(NormRatio^3);
Lambda = Lambda/(NormRatio^2);
TrueTotTime = sum(Lambda);
CurveLoss = TrueTotTime*(VUnitCell/TotVgrowth);

%% calculating the pore size loss function and total loss function
PoreSizeLoss = -aPar*log(abs(Radius - bPar));
TotLoss = CurveLoss + PoreSizeLoss;
if ~isreal(TotLoss)
    CurveLoss = 0; PoreSizeLoss = 0; TotLoss = 10^5;
end
% fileID = fopen(strcat(ParDir,'\Loss.txt'), 'w+');
% fprintf(fileID,'%4.7f\n',CurveLoss);
% fprintf(fileID,'%4.7f\n',PoreSizeLoss);
% fprintf(fileID,'%4.7f\n',TotLoss);
% fclose(fileID);