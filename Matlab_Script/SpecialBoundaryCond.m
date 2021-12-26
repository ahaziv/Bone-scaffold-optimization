function [Vertices, BoundStrVer] = SpecialBoundaryCond(Vertices, BoundStrVer)
% this function checks for special boundary conditions on the function and
% adjusts boundaries/locations if necessary.
%% checking for vanishing boundaries
ii=1;
while ii <= length(BoundStrVer)
    if length(BoundStrVer{ii}) < 4
        Vertices(BoundStrVer{ii},4:6) = 0;
        BoundStrVer = [BoundStrVer(1:ii-1); BoundStrVer(ii+1:end)];
        continue
    end
    ii=ii+1;
end