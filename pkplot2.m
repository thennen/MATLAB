function pkplot(PKdata,X,Y,varargin)
% Plots cell arrays contained in a structure PKData  PKData.(X){i} vs
% PKData.(Y){i} for all i
% varargin{1} can specify colormap
%
% T Hennen 2013

if(~ishold)
    washeld = false;
    axes1 = axes;
    box(axes1,'on');
    hold(axes1,'all');
    xlabel(X);
    ylabel(Y);
else
    washeld = true;
end

if ~isfield(PKdata,X)|| ~isfield(PKdata,Y)
    error('reference to nonexistent field')
elseif length(PKdata.(X)) ~= length(PKdata.(Y))
    error('cell arrays X and Y must have the same length')
end
nplots = size(PKdata.(X),2);

% Default colormap
colormap = [0 0 1];

if nargin == 4 %colormap specified
    colormap = varargin{1};
    if size(colormap,2) == 3
        colormap = varargin{1};
    else
        error('invalid colormap')
    end    
elseif nargin > 4
    error('too many input arguments');
end
ncolors = size(colormap,1);

%Interpolate color list for nplots
if ncolors > 1
    interpvals = linspace(1,ncolors,nplots);
    interpR = interp1(colormap(:,1),interpvals);
    interpG = interp1(colormap(:,2),interpvals);
    interpB = interp1(colormap(:,3),interpvals);
else
    PlotColor = colormap;
end

hold all;
for i=1:nplots
    if ncolors > 1
        PlotColor = [interpR(i) interpG(i) interpB(i)];
    end
    plot(PKdata.(X)(:,i),PKdata.(Y)(:,i),'color',PlotColor,'linewidth',2,'markersize',20);
end

if(~washeld)
    hold(axes1,'off');
end

end
