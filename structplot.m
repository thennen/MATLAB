function structplot(searchString,xvect,yvect,varargin)
% Plots struct.(xvect) vs struct.(yvect) for all structures containing searchString in their name
% if xvect and yvect have length one, all values will be plotted as one
% varargin can specify a color map

if(~ishold)
    washeld = false;
    fig = figure();
    %axes1 = axes;
    %box(axes1,'on');
    %hold(axes1,'all');
    xlabel(PlotLabel(xvect));
    ylabel(PlotLabel(yvect));
    title(strrep(searchString,'_',' '));
else
    washeld = true;
end

% Default colormap
colormap = [0 0 1; 0 1 0];  % Blue to green

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

varnames = findbasevars(searchString,'struct');
%% Look for valid plots
Xmatrix = cell(length(varnames), 1);
Ymatrix = cell(length(varnames), 1);
for i = length(varnames):-1:1
    if isbasefield(varnames{i},xvect) && isbasefield(varnames{i},yvect)
        Xmatrix{i} = evalin('base',[varnames{i} '.' xvect]);
        Ymatrix{i} = evalin('base',[varnames{i} '.' yvect]);

        if iscell(Xmatrix{i}) % Assume Ymatrix is also cell
        % If X and Y fields are cell arrays (assume of 1d arrays), convert to matrix
        % padded with NaN
            xlength = cellfun(@length,Xmatrix{i});
            ylength = cellfun(@length,Ymatrix{i});
            if length(Xmatrix{i}) == length(Ymatrix{i}) && all(xlength == ylength)
                
                maxlength = max(xlength);
                
                for j=1:length(Xmatrix{i})
                    Xmatrix{i}{j}(end+1:maxlength) = NaN;
                    Ymatrix{i}{j}(end+1:maxlength) = NaN;
                end
                Xmatrix{i} = cell2mat(Xmatrix{i}.').';
                Ymatrix{i} = cell2mat(Ymatrix{i}.').';
                
            else
                warning(['Cell array size mismatch in' varnames{i}])
                varnames(i) = [];
                Xmatrix(i) = [];
                Ymatrix(i) = [];
            end
        end
        
        
        % Ensure xvect and yvect have the same size
        if size(Xmatrix{i}) ~= size(Ymatrix{i})
            warning(['xvect and yvect have different lengths in' varnames{i}]);
            varnames(i) = [];
            Xmatrix(i) = [];
            Ymatrix(i) = [];
        end
        
    else
        warning(['Specified field does not exist in' varnames{i}]);
        varnames(i) = [];
        Xmatrix(i) = [];
        Ymatrix(i) = [];
    end
end

nplots = length(varnames);
if nplots == 0;
    error('No valid plots matching search criteria')
end

%%Interpolate color list for nplots
if ncolors > 1
    interpvals = linspace(1,ncolors,nplots);
    interpR = interp1(colormap(:,1),interpvals);
    interpG = interp1(colormap(:,2),interpvals);
    interpB = interp1(colormap(:,3),interpvals);
    PlotColor = [interpR(1) interpG(1) interpB(1)];
else
    PlotColor = colormap;
end

if 	all(cellfun(@length,Ymatrix) == 1)
    % if length of all yvects are 1, then combine the points from all the
    % matching structures into one plot
    plot(cell2mat(Xmatrix),cell2mat(Ymatrix),'.-','DisplayName',searchString,'Color',PlotColor(1,:));
else
    for i = 1:nplots
        if ncolors > 1
            PlotColor = [interpR(i) interpG(i) interpB(i)];
        end
        
        % Decide display name
        if isbasefield(varnames{i},'DisplayName')
            DisplayName = evalin('base',[varnames{i} '.' DisplayName]); %DisplayName needs to be a string
        else
            DisplayName = strrep(varnames{i},'_','');
        end
        
        plot(Xmatrix{i},Ymatrix{i},'.-','DisplayName',DisplayName,'Color',PlotColor,'LineWidth',2);
    end
end

%legend('toggle');
if(~washeld)
    hold(axes1,'off');
end

end