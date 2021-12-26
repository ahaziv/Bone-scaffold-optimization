function [Faces, Vertices, BoundStrVer] = Remesh(Faces, Vertices, MinDeltaL, DeltaTheta, MaxDeltaL, BoundStrVer)
% this function adapts a given mesh by merging vertices that are closer
% than a certain distance.
% Faces - an array containing the vertices numbers of each face (triangles)
% Vertices - an array cotaining all the vertices coordinates
% MInDeltaL - the distance under which merging of vertices is executed
% MInDeltaL - the angle square value under which merging of vertices is executed
% MaxDeltaL - the distance above which splitting of vertices is executed

%% lets do it! 
[emag, theta] = ecalc(Faces(:,1:3),Vertices(:,1:3));
%% merging/adding vertices according to length/angle constraints
IndexArray = [2 3; 3 1; 1 2];
ii = 1;
if nargin < 6
    BoundStrVer = 0;
end
while ii < length(emag(:,1))
    for jj = 1:3
        % deleting faces and vertices
        if emag(ii,jj) < MinDeltaL || theta(ii,jj) < DeltaTheta
            [Vertices, Faces, BoundStrVer] = DeleteFaces(Vertices, Faces, IndexArray(jj,1) ,IndexArray(jj,2), ii, BoundStrVer);
            [emag, theta] = ecalc(Faces(:,1:3), Vertices(:,1:3));
            break
        end
        % adding faces and vertices
        if emag(ii,jj) > MaxDeltaL
            [Vertices, Faces] = AddFaces(Vertices, Faces, BoundStrVer, IndexArray(jj,1), IndexArray(jj,2), ii);
            [emag, theta] = ecalc(Faces(:,1:3),Vertices(:,1:3));
            break
        end
    end
    ii=ii+1;
end

%% deleting faces in case of congruent/protruding faces
while 1
    Flag = 0;
    Edges = zeros(3*length(Faces(:,1)),2);
    Edges(1:3:3*length(Faces(:,3))-2,1:2) = [Faces(:,3) Faces(:,2)];
    Edges(2:3:3*length(Faces(:,3))-1,1:2) = [Faces(:,1) Faces(:,3)];
    Edges(3:3:3*length(Faces(:,3))  ,1:2) = [Faces(:,2) Faces(:,1)];
    % deleting degenerate faces (zero surface)
    temp = find(Edges(:,1)==Edges(:,2));
    if temp
        DeletedFaces = Faces(ceil(temp/3),:);
        Faces(ceil(temp/3),:) = [];
        % if face deletion cast out a Vertex - delete it
        for ii=1:length(DeletedFaces(:,1))
            for jj = 1:3
                if ~find(Faces(:,:)==Faces(ii,jj))
                    [Vertices, Faces, BoundStrVer] = DeleteVertex(Vertices, Faces, Faces(ii,jj), BoundStrVer);
                end
            end
        end
        continue
    end
    % deleting degenerate faces (zero Volume)
    SortedEdges = sort(Edges,2);
    UniqueEdges = SortedEdges;
    [~, UniqueIndeces,~] = unique(UniqueEdges(:,:),'rows');
    UniqueEdges(UniqueIndeces,:) = [];
    [~, UniqueIndeces,~] = unique(UniqueEdges(:,:),'rows');
    UniqueEdges(UniqueIndeces,:) = [];
    if ~isempty(UniqueEdges)
        %% deleting the second out of each two identical faces
        SortedFaces = sort(Faces(:,1:3),2);
        [u,I,~] = unique(SortedFaces, 'rows');
        hasDuplicates = size(u,1) < size(SortedFaces,1);
        if hasDuplicates
            ixDupRows = setdiff(1:size(SortedFaces,1), I);
            DeletedFaces = Faces(ixDupRows,:);
            Faces(ixDupRows,:)= [];
            % if face deletion cast out a Vertex - delete it
            for ii = 1:length(DeletedFaces(:,1))
                for jj = 1:3
                    if ~find(Faces(:,:)==Faces(ii,jj))
                        [Vertices, Faces, BoundStrVer] = DeleteVertex(Vertices, Faces, Faces(ii,jj), BoundStrVer);
                    end
                end
            end
            continue
        end
        %% deleting all the single protruding faces
        for ii = 1:length(UniqueEdges(:,1))
            for jj = 1:length(SortedEdges(:,1))
                if SortedEdges(jj,1) == UniqueEdges(ii,1) && SortedEdges(jj,2) == UniqueEdges(ii,2)
                    SneakyFace = ceil(jj/3); 
                    ThirdVertex = setdiff(Faces(SneakyFace,1:3),SortedEdges(jj,:));
                    if sum(sum(Faces(:,1:3) == ThirdVertex)) == 1
                        Faces(SneakyFace,:) = [];
                        % in the case an Evil Vertex is secluded - delete it
                        if ~find(Faces(:,:)==ThirdVertex)
                            [Vertices, Faces, BoundStrVer] = DeleteVertex(Vertices, Faces, ThirdVertex, BoundStrVer);
                        end
                        Flag = 1;
                        break
                    end
                end
            end
            if Flag
                break
            end
        end
        if Flag
            continue
        end 
    end
    break
end
Vertices = Vertices(:,1:7);
end

%% DeleteFaces function
function [Vertices, Faces, BoundStrVer] = DeleteFaces(Vertices, Faces, KeptIndex ,DiscardedIndex, ii, BoundStrVer)
VerKeptIndex = Faces(ii,KeptIndex);
VerDiscardedIndex = Faces(ii,DiscardedIndex);
if VerKeptIndex == VerDiscardedIndex
    return
end
Vertices(VerKeptIndex,end+1) = 1; % this vertex is kept
Vertices(VerKeptIndex,1:3) = (Vertices(VerKeptIndex,1:3) + Vertices(VerDiscardedIndex,1:3))/2;

Vertices(VerKeptIndex,4:6) = max(Vertices(VerKeptIndex,4:6),Vertices(VerDiscardedIndex,4:6));

Faces(Faces(:,1:3) == VerDiscardedIndex) = VerKeptIndex;
% deleting the two faces
Faces(ii,:) = [];
for jj = 1:length(Faces(:,1))
    if sum(Vertices(Faces(jj,1:3),end)) == 2
        Faces(jj,:) = [];
        break
    end
end
% deleting the Vertex out of the system
[Vertices, Faces, BoundStrVer] = DeleteVertex(Vertices, Faces, VerDiscardedIndex, BoundStrVer);

Vertices(:,end) = []; % reseting the array
end

%% AddFaces function
% this functio adds faces in cases of long vertices
function [Vertices, Faces, BoundStrVer] = AddFaces(Vertices, Faces, BoundStrVer, Index1 ,Index2, Face1Index)
% Face1Index - the row index in Faces containing the edge.
% Index1/2 - the j indices pointing to the long edge.
KeptSideIndex = Faces(Face1Index,Index1);
NewSideIndex  = Faces(Face1Index,Index2);
Vertices(KeptSideIndex,7) = 1; Vertices(NewSideIndex,7) = 1;
% finding the other Face index
Face2Index = 0;
for jj = 1:length(Faces(:,1))
    if sum(Vertices(Faces(jj,1:3),7)) == 2
        if jj ~= Face1Index
            Face2Index = jj;
            break
        end
    end
end
% defining a new vertex and two new faces
Vertices(length(Vertices(:,1))+1,:) = (Vertices(KeptSideIndex,:) + Vertices(NewSideIndex,:))/2;
% if the cut edge is made of two border vertices, mark the new vertex as a
% border one and add it to BoundStrVer
Vertices(end,4:6) = floor(Vertices(end,4:6));
% if sum(Vertices(end,4:6))
%     for ii = 1:length(BoundStrVer)
%         index = find(BoundStrVer{ii}==KeptSideIndex);
%         if index
%             if index == 1
%                 BoundStrVer{ii} = [BoundStrVer{ii}(:); length(Vertices(:,1))];
%             elseif index == length(BoundStrVer{ii})
%                 BoundStrVer{ii} = [length(Vertices(:,1)); BoundStrVer{ii}(:)]; 
%             else
%                 BoundStrVer{ii} = [BoundStrVer{ii}(1:index-1); length(Vertices(:,1)); BoundStrVer{ii}(index:end)];
%             end
%         end
%     end
% end
if Face2Index                           % regular edge case
    NewFaces = [Faces(Face1Index,:); Faces(Face2Index,:)];

    Faces(Face1Index,Faces(Face1Index,1:3) == NewSideIndex) = length(Vertices(:,1));
    Faces(Face2Index,Faces(Face2Index,1:3) == NewSideIndex) = length(Vertices(:,1));
    NewFaces(1,NewFaces(1,:) == KeptSideIndex) = length(Vertices(:,1));
    NewFaces(2,NewFaces(2,:) == KeptSideIndex) = length(Vertices(:,1));
else                                    %  border edge case
    NewFaces = Faces(Face1Index,:);
    
    Faces(Face1Index,Faces(Face1Index,1:3) == NewSideIndex) = length(Vertices(:,1));
    NewFaces(1,NewFaces(1,:) == KeptSideIndex) = length(Vertices(:,1));
end

Faces = [Faces(1:Face1Index,:); NewFaces; Faces(Face1Index+1:end,:)];

Vertices(:,7) = 0; % reseting the array
end

%% DeleteVertex function
function [Vertices, Faces, BoundStrVer] = DeleteVertex(Vertices, Faces, VerIndex, BoundStrVer)
Vertices(VerIndex,:) = [];
Faces(Faces(:,1:3)>VerIndex) = Faces(Faces(:,1:3)>VerIndex)-1;
if iscell(BoundStrVer)
    for jj=1:length(BoundStrVer)
        BoundStrVer{jj}(BoundStrVer{jj}(:,:)==VerIndex) = [];
        BoundStrVer{jj}(BoundStrVer{jj}(:,:)>VerIndex) = BoundStrVer{jj}(BoundStrVer{jj}(:,:)>VerIndex)-1;
    end
end
end