function [Vertices, BoundaryVal, BoundStrVer] = DetectBoundaries(Faces, Vertices, eps)
% this function goes over the entire faces array and detects edge faces,
% the faces' edge vertices are then marked.
% Vertices(:,4) is a flag column where 0 marks regular vertex and 1 marks an
% edge vertex.
BoundaryVal = zeros(6,1);
BoundaryVal(1) = min(Vertices(:,1));
BoundaryVal(2) = max(Vertices(:,1));
BoundaryVal(3) = min(Vertices(:,2));
BoundaryVal(4) = max(Vertices(:,2));
BoundaryVal(5) = min(Vertices(:,3));
BoundaryVal(6) = max(Vertices(:,3));

Edges = zeros(3*length(Faces(:,1)),2);
Edges(1:3:3*length(Faces(:,3))-2,1:2) = [Faces(:,3) Faces(:,2)];
Edges(2:3:3*length(Faces(:,3))-1,1:2) = [Faces(:,1) Faces(:,3)];
Edges(3:3:3*length(Faces(:,3))  ,1:2) = [Faces(:,2) Faces(:,1)];
%% setting the appropriate boundary for each vertex according to it's XYZ values
% BoundaryEdges = zeros(length(Edges(:,1)),2);
%% detect boundaries - new version 26/10
SortedEdges = sort(Edges,2);
[UniqueEdges,~,RowIdx] = unique(SortedEdges,'rows','stable');
[C,~] = hist(RowIdx,unique(RowIdx));
[~, BoundEdgesIdx] = find(C==1);
BoundaryEdges = UniqueEdges(BoundEdgesIdx.',:);
% reorganizing the data in BoundaryEdges to match the original Edges array
for ii = 1:length(BoundaryEdges(:,1))
    [~,RowTempIdx] = ismember(BoundaryEdges(ii,:),SortedEdges,'rows');
    BoundaryEdges(ii,:) = Edges(RowTempIdx,:);
end
BoundVerIdx = unique(reshape(BoundaryEdges, [numel(BoundaryEdges) 1]));

for ii = 1:length(BoundVerIdx)
    for jj = 1:6
        if abs(Vertices(BoundVerIdx(ii),ceil(0.5*jj)) - BoundaryVal(jj)) < eps
            Vertices(BoundVerIdx(ii),ceil(0.5*jj)) = BoundaryVal(jj);
            Vertices(BoundVerIdx(ii),ceil(0.5*jj)+3) = 1;
        end
    end
end

%% detect boundaries - old version before 26/10
% while ii <= length(Edges(:,1))
%     flag = 1;
%     for jj = ii+1:length(Edges(:,1))
%         if Edges(jj,1) == Edges(ii,2) && Edges(jj,2) == Edges(ii,1)
% %             Edges = [Edges(1:jj-1,:); Edges(jj+1:end,:)];
%             Edges(jj,:) = [];
%             flag = 0;
%             break
%         end
%     end
%     if flag == 1
%         kk = kk+1;
%         BoundaryEdges(kk,:) = Edges(ii,:);
%         for jj=1:3
%             if (Vertices(Edges(ii,1),jj) == BoundaryVal(2*jj-1) || Vertices(Edges(ii,1),jj) == BoundaryVal(2*jj))
%                 Vertices(Edges(ii,1),jj+3) = 1;
%             end
%         end
%     end
%     ii = ii+1;
% end
% BoundaryEdges = BoundaryEdges(1:kk,:);
%% creating a cell array containing arrays of the different boundaries
BoundStrVer = cell(10,1);
ii = 0;
while ~isempty(BoundaryEdges)
    ii=ii+1;
    BoundStrVer{ii} = zeros(200,1);
    BoundStrVer{ii}(1:2) = BoundaryEdges(1,:);
    BoundaryEdges = BoundaryEdges(2:end,:);
    jj = 2;
    while BoundStrVer{ii}(jj) ~= BoundStrVer{ii}(1)
        jj=jj+1;
        temp = find(BoundaryEdges(:,1) == BoundStrVer{ii}(jj-1));
        if ~temp
            temp = find(BoundaryEdges(:,2) == BoundStrVer{ii}(jj-1));
        end
        BoundStrVer{ii}(jj) = BoundaryEdges(temp,2);
        BoundaryEdges = [BoundaryEdges(1:temp-1,:); BoundaryEdges(temp+1:end,:)];
    end
    BoundStrVer{ii} = BoundStrVer{ii}(1:jj);
end
BoundStrVer = BoundStrVer(1:ii);