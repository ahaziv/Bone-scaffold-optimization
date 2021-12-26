function [CurvatureVec, ColorVec, AnaSolution] = CurveMethod(CurveCalcMethod, VertexSFM, AnalyticParam, PartVertices)
% this function calculates an invariant of the curvature at each vertex given a curvature matrix at each point
% CurveCalcMethod - the curvature calculation method, (there are many different invariants to choose from)
% VertexSFM - An array containing the curvature matrix at each vertex
% PartVertices - The list of participating Vertices in the partial view case
% CurvatureVec - the output curvature vector
% AmaSolution - in case an analytic solution exists this parameter is an
% array containing: [case, Shape Radius]
switch CurveCalcMethod
    case 'Average'
        AvgCurvature = squeeze(VertexSFM(1,1,:) + VertexSFM(2,2,:)); % not average as we use to see
        ModVerNorm = AvgCurvature < 0;
        CurvatureVec = AvgCurvature.*ModVerNorm;
        ColorVec = AvgCurvature.';
        if nargin > 2 % this checks if an analytic solution for the case exists
           if ~isempty(PartVertices)
               AnaSolution = AnalyticSolution(VertexSFM(:,:,PartVertices), AnalyticParam);
           else
               AnaSolution = AnalyticSolution(VertexSFM, AnalyticParam);
           end
        else
            AnaSolution = 0;
        end
    case 'Normal'  % based upon the equation: k=sqrt(a1*|a1|+a2*|a2|)
        NormCurvature = zeros(length(VertexSFM(1,1,:)),1);
        for ii = 1:length(NormCurvature)
            eigenvals = eig(VertexSFM(:,:,ii));
            NormCurvature(ii) = eigenvals(1)*abs(eigenvals(1))+eigenvals(2)*abs(eigenvals(2));
        end
        ModVerNorm = NormCurvature < 0;
        CurvatureVec = -(-NormCurvature.*ModVerNorm).^0.5;
        ColorVec = (CurvatureVec+(NormCurvature.*~ModVerNorm).^0.5).';
        if nargin > 2 % this checks if an analytic solution for the case exists
           if ~isempty(PartVertices)
               AnaSolution = AnalyticSolution(VertexSFM(:,:,PartVertices), AnalyticParam);
           else
               AnaSolution = AnalyticSolution(VertexSFM, AnalyticParam);
           end
        else
            AnaSolution = 0;
        end
    case 'ThirdMainInv' % NOT WORKING!
        AvgCurvature = (VertexSFM(1,:) + VertexSFM(2,:));
        GausianCurvature = abs(VertexSFM(1,:)).*VertexSFM(2,:);    
        ThirdInvariantCurvature = (sum(VertexSFM(:,:).^3,1)+3*GausianCurvature);
        ThirdInvariantCurvature = min(ThirdInvariantCurvature,0);                   % killing concave normals 
        CurvatureVec = (-ThirdInvariantCurvature).^(1/3);
        ColorVec = AvgCurvature.';
end
end

function [AnaSolution] = AnalyticSolution(VertexSFM, AnalyticParam)
% AnaSolution - an (1,2) array, the first component is the solution type
% (1==Sphere, 2==Cylinder) and the second component the Shape radius
AvgCurvature = squeeze(VertexSFM(1,1,:) + VertexSFM(2,2,:)); % not average as we use to see   
STDAvgCurve = std(AvgCurvature);                            % checking if the AvgCurvature has sexual transmitted disease
if STDAvgCurve < AnalyticParam.GenEps
    GausCurvature = zeros(length(AvgCurvature),1);
    for ii = 1:length(AvgCurvature)
        PrincipalCurv = eig(VertexSFM(:,:,ii));
        GausCurvature(ii) = PrincipalCurv(1,:).*PrincipalCurv(2,:);    
    end
    TotAvgCurvature = sum(abs(AvgCurvature))/length(AvgCurvature);
    TotGausCurvature = sum(abs(GausCurvature))/length(AvgCurvature);
    % checking which of the two cases if met:
    if abs(TotAvgCurvature/2 - sqrt(TotGausCurvature)) < AnalyticParam.SphereEps      % sphere
        Radius = 2/TotAvgCurvature;
        AnaSolution = [1 Radius];
%         AnaSolution = [1 (Radius^2)/4];
%         AnaSolution = Radius/(4*AnalyticParam.LFactor);
        disp('A Spherical analytic solution was achieved m`lord.')
    elseif TotGausCurvature/TotAvgCurvature < AnalyticParam.CylinderEps               % cylinder
        Radius = 1/TotAvgCurvature;
        AnaSolution = [2 Radius];
%         AnaSolution = [2 (Radius^2)/2];
%         AnaSolution = Radius/(4*AnalyticParam.LFactor);
        disp('A Cylindrical analytic solution was achieved m`lord.')
    else
        error('An analytic solution was reached, yet unable to discern which case was met. Please adjust the analytic parameters.')
    end
else
    AnaSolution = 0;
end
end