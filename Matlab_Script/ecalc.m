% this function calculates the magnitude of e0, e1, e2 out of a given set
% of faces and vertices
function [emag,theta] = ecalc(Faces,Vertices)

e0 = Vertices(Faces(:,3),:)-Vertices(Faces(:,2),:);
e1 = Vertices(Faces(:,1),:)-Vertices(Faces(:,3),:);
e2 = Vertices(Faces(:,2),:)-Vertices(Faces(:,1),:);
emag = [sqrt(sum(e0.^2,2)) sqrt(sum(e1.^2,2)) sqrt(sum(e2.^2,2))];

% in case we want to estimate the angles between vectors
% according to the formula: theta^2=sin^2(theta)=1-(dot(A,B)/|A|*|B|)^2
if nargout == 2
    theta = zeros(size(emag));
    theta(:,1) = ones(1,length(e0))-(dot(e1.',e2.')./((emag(:,2).*emag(:,3)).')).^2;
    theta(:,2) = ones(1,length(e0))-(dot(e0.',e2.')./((emag(:,1).*emag(:,3)).')).^2; 
    theta(:,3) = ones(1,length(e0))-(dot(e0.',e1.')./((emag(:,1).*emag(:,2)).')).^2;
end


% old script
% if nargout == 2
%     theta = zeros(size(emag));
%     theta(:,1) = ones(1,length(e0))+dot(e1.',e2.')./((emag(:,2).*emag(:,3)).');
%     theta(:,2) = ones(1,length(e0))+dot(e0.',e2.')./((emag(:,1).*emag(:,3)).'); 
%     theta(:,3) = ones(1,length(e0))+dot(e0.',e1.')./((emag(:,1).*emag(:,2)).');
% end