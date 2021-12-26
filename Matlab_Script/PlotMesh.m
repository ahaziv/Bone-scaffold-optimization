function [] = PlotMesh(Faces, Vertices, ColorVec, maxis, minmaxColor, FlagBoundary, FlagPartView, InitGeom)
% this function plots a given mesh
% InitGeom - a 1x2 cell containing the initial Faces and Vertices arays to display ate the background

if FlagPartView
    if ~isempty(InitGeom)
        trisurf(InitGeom{1}(InitGeom{1}(:,4)==1,1:3),InitGeom{2}(:,1),InitGeom{2}(:,2),InitGeom{2}(:,3),'facecolor','c','FaceAlpha',0.2,'edgecolor','none');
        hold on
    end
    if isempty(ColorVec)
        mesh1 = trisurf(Faces(Faces(:,4)==1,1:3),Vertices(:,1),Vertices(:,2),Vertices(:,3),'edgecolor','interp');
    else
        mesh1 = trisurf(Faces(Faces(:,4)==1,1:3),Vertices(:,1),Vertices(:,2),Vertices(:,3),'edgecolor','interp','FaceVertexCdata',-ColorVec');
    end
else
    if ~isempty(InitGeom)
        trisurf(InitGeom{1}(:,1:3),InitGeom{2}(:,1),InitGeom{2}(:,2),InitGeom{2}(:,3),'facecolor','c','FaceAlpha',0.2,'edgecolor','none');
        hold on
    end
    if isempty(ColorVec)
        mesh1 = trisurf(Faces,Vertices(:,1),Vertices(:,2),Vertices(:,3),'edgecolor','interp');
    else
        mesh1 = trisurf(Faces,Vertices(:,1),Vertices(:,2),Vertices(:,3),'edgecolor','interp','FaceVertexCdata',-ColorVec');
    end
    if FlagBoundary
        EdgeVertices = find(sum(Vertices(:,4:6),2) ~= 0);
        hold on
        plot3(Vertices(EdgeVertices,1),Vertices(EdgeVertices,2),Vertices(EdgeVertices,3),'.','Color','r','markersize',10);
        plot3(Vertices(sum(Vertices(:,4:6),2)==2,1),Vertices(sum(Vertices(:,4:6),2)==2,2),Vertices(sum(Vertices(:,4:6),2)==2,3),'.','Color','g','markersize',10)
%         plot3(P(:,1),P(:,2),P(:,3),'.','Color','b','markersize',10) 
    end
end

hold off
axis off
set(mesh1,'ambientstrength',0.35);
if ~isempty(ColorVec)
    caxis([minmaxColor(1) minmaxColor(2)])
end
camlight();
lighting gouraud       
xlim([-maxis maxis]);
ylim([-maxis maxis]);
zlim([-maxis maxis]);
view([1 1 1])