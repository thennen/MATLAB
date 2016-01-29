function cellplot(xcell,ycell)
%CELLPLOT plot ycell{i} vs xcell{i} for each i
% Assume that ycell{i} and xcell{i} are arrays of the same length.

figure;
hold all;
for i=1:length(xcell)
    plot(xcell{i},ycell{i})
end
hold off;
end

