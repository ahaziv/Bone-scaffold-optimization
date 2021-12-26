%% This is the main code for simulating CCTG 
clear all; close all; clc;

addpath(fullfile(cd, '..'))
addpath(fullfile(cd))

%% input parameters
FW = 5*10^-4;              % fillament width
FD = 5.5*10^-4;              % fillament distance
ST = 3*10^-4;               % slice thickness

aPar = 7;                       % the true parameter is 7*10^-8 but was normalized
bPar = 50*10^-4;                % the true parameter is 50*10^-6 but was normalized

%% general parameters
NumIter = 500;              % number of iteration steps
TimeLimit = 0.3;
LFactor = 0.24;              % Lambda multiplication factor
DeltaTheta = 0.08;          % the square of the minimal angle for remeshing 
FlagBoundary = 1;           % this flag tells the program wether to impose boundary conditions
FlagAnaSol = 1;
AnaSolution = 1;
if FlagAnaSol
    AnalyticParam = struct('GenEps', 0.5, 'SphereEps', 0.2, 'CylinderEps', 0.5);
end

Faces = cell(NumIter,1);
Vertices = cell(NumIter,1);
CurveCalcMethod = 'Average';

Vgrowth = zeros(NumIter,1); % this array stores the filled up volume at each iteration
Lambda = zeros(NumIter,1);
TotTime = 0;                % the total time of the simulation
Vtot = 2*ST*FD^2;
%% Generate the example
[Faces{1}, Vertices{1}, NormRatio] = CreateShape();

Vertices{1}(1,4:7) = 0;
[emag] = ecalc(Faces{1},Vertices{1}(:,1:3)); 
DeltaRemesh = min(min(emag));
DeltaRemeshMax = 1.3*max(max(emag));
maxis = 1.3*max(max(abs(Vertices{1}(:,1:3))));

%% detect mesh edges for later use
if FlagBoundary
    [Vertices{1}, BoundaryVal, BoundStrVer] = DetectBoundaries(Faces{1},Vertices{1}, 0.5*DeltaRemesh);
end
%% saving the initial geometry for display
if FlagDispInitGeom
    InitGeom = {Faces{1} Vertices{1}};
else
    InitGeom = [];
end
%% open a figure
fig_h = figure('name','Scaffold','numbertitle','off','color',[0.75 0.75 0.75]);
ax = axes('DataAspectRatio', [1,1,1]);
set(gcf, 'OuterPosition', get(0, 'Screensize'));
axis tight
%% run the algorithm
TotVgrowth = 0;
ii = 1;
while TotTime < TimeLimit 
    if FlagBoundary % reflecting the boundaries to calculate the curvatures and normals
        OrigFacesLength = length(Faces{ii}(:,1));
        OrigVertiLength = length(Vertices{ii}(:,1));
        for ll = 1:6
            [RefFaces, RefVertices] = ReflectBoundaryMesh(Faces{ii}, Vertices{ii}, [ceil(0.5*ll) BoundaryVal(ll)], 0.5*DeltaRemesh);
            Faces{ii} = [Faces{ii}; RefFaces];
            Vertices{ii} = [Vertices{ii}; RefVertices];
        end
    end
    % calculating the curvatures and normals
    [FaceNormals] = CalcFaceNormals(Faces{ii}(:,1:3),Vertices{ii}(:,1:3));
    [VertexNormals, Avertex, Acorner,up,vp] = CalcVertexNormals(Faces{ii}(:,1:3),Vertices{ii}(:,1:3),FaceNormals);
    [VertexSFM, wfp] = CalcCurvature(Faces{ii}(:,1:3),Vertices{ii}(:,1:3),VertexNormals,FaceNormals,Avertex,Acorner,up,vp);
            
    if FlagBoundary % after calculating the curvatures and normals, cut off the reflected entities
        Faces{ii} = Faces{ii}(1:OrigFacesLength,:);
        Vertices{ii} = Vertices{ii}(1:OrigVertiLength,:);
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
    Vertices{ii+1} = Vertices{ii} - [Step, zeros(size(Vertices{ii},1),4)];
    TotTime = TotTime + Lambda(ii);
    %% Calculating the total volume filled up in this iteration
    AvertexNew = CalcVertexAreas(Faces{ii}(:,1:3),Vertices{ii+1}(:,1:3));
    Vgrowth(ii) = sum(((Avertex+AvertexNew)/2).*sqrt(sum(Step.^2,2)));
    TotVgrowth = TotVgrowth + Vgrowth(ii)/(NormRatio^3);
    if TotVgrowth > 0.97*Vtot
        TotVgrowth = Vtot;
        break
    end
    %% Remeshing the shape to avoid clustering of vertices
    if FlagBoundary
        [Faces{ii+1}, Vertices{ii+1}, BoundStrVer] = Remesh(Faces{ii}, Vertices{ii+1}, DeltaRemesh, DeltaTheta, DeltaRemeshMax, BoundStrVer);
        % imposing special boundary conditions
        [Vertices{ii+1}, BoundStrVer,P] = SpecialBoundaryCond(Vertices{ii+1} ,BoundStrVer);
    else
        [Faces{ii+1}, Vertices{ii+1}] = Remesh(Faces{ii}, Vertices{ii+1}, DeltaRemesh, DeltaTheta, DeltaRemeshMax);
    end  
    %% Plotting the mesh
    if ~mod(ii,1)   
        PlotMesh(Faces{ii}, Vertices{ii}, ColorVec, maxis, minmaxColor, FlagBoundary, FlagPartView, InitGeom)
        drawnow
    end

    ii = ii+1;
end
%% Reorganizing the output
Faces = Faces(1:ii-1);
Vertices = Vertices(1:ii-1);
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
CurveLoss = TrueTotTime*(Vtot/TotVgrowth);

%% calculating the pore size loss function and total loss function
PoreRad = 0.25*(ST^2+(FD-FW)^2)^(1/2);
PoreSizeLoss = -aPar*log(PoreRad - bPar);
TotLoss = CurveLoss + PoreSizeLoss;
if ~isreal(TotLoss)
    TotLoss = -1;
end