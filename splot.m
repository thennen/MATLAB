function splot(searchString,xvect,yvect,varargin)
% Plots struct.(xvect) vs struct.(yvect) for all structures containing 
% searchString in their name.  if xvect and yvect have length one, all
% values will be plotted as one.
%
% Varargin can specify a color map
% e.g.
% RGB fade:
% splot('PKDynamic', 'Field', 'Kerr', [1 0 0; 0 1 0; 0 0 1])
% All blue:
% splot('PKDynamic', 'Field', 'Kerr', [1 0 0; 0 1 0; 0 0 1])
% Builtin map:
% splot('PKDynamic', 'Field', 'Kerr', jet)
%
%
% Can input struct directly

% TO DO:
% xvect and yvect can refer to nested structures
%
% xvect keywords:
% 'ii' : plot y (if length(y)=1) vs variable number

% Default colormap
colormap = [0 0 1; 0 1 0];  % Blue to green

if nargin == 4
    %colormap specified
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

if isstruct(searchString)
    % set varnames to keyword 'Direct' so that the rest of the function 
    % knows what to do
    varnames = {'Direct'};
    Direct = true;
elseif ischar(searchString)
    varnames = findbasevars(searchString,'struct');
    Direct = false;
elseif iscell(searchString)
    % Matlab sucks
    Direct = false;
    varnames = {};
    for i=1:size(searchString,1)
        if evalin('base',['exist(''' searchString{i} ''')'])
            varnames{end+1} = searchString{i};
        end
    end
else
    error('First argument must be either a search string, cell of strings, or a struct.')
end


%% Look for valid plots
Xmatrix = cell(length(varnames), 1);
Ymatrix = cell(length(varnames), 1);
% Delete varnames for which the xvect and yvect fields don't both exist
for i = length(varnames):-1:1
    if Direct && isfield(searchString,xvect) && isfield(searchString,yvect)
        Xmatrix{i} = searchString.(xvect);
        Ymatrix{i} = searchString.(yvect);
    elseif ~Direct && strcmp(xvect,'ii') && isbasefield(varnames{i},yvect)
        Ymatrix{i} = evalin('base',[varnames{i} '.' yvect]);
        Xmatrix{i} = i*ones(size(Ymatrix{i},1),1);
    elseif ~Direct && isbasefield(varnames{i},xvect) && isbasefield(varnames{i},yvect)
        Xmatrix{i} = evalin('base',[varnames{i} '.' xvect]);
        Ymatrix{i} = evalin('base',[varnames{i} '.' yvect]);
    else
        disp(['A specified field does not exist in ' varnames{i}]);
        varnames(i) = [];
        Xmatrix(i) = [];
        Ymatrix(i) = [];
    end
end
% Delete varnames for which there are size mismatches
for i = length(varnames):-1:1
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
                warning(['Cell array size mismatch in ' varnames{i}])
                varnames(i) = [];
                Xmatrix(i) = [];
                Ymatrix(i) = [];
            end
        end
        
        % Ensure xvect and yvect have the same size
        if size(Xmatrix{i}) ~= size(Ymatrix{i})
            warning(['xvect and yvect have different lengths in ' varnames{i}]);
            varnames(i) = [];
            Xmatrix(i) = [];
            Ymatrix(i) = [];
        end
end

nplots = length(varnames);
if nplots == 0;
    error('No valid plots matching search criteria')
elseif nplots == 1;
    % If there's only one struct to plot from, then split it up and apply
    % color map to columns, if any
    nplots = size(Xmatrix{1},2);
    for w=nplots:-1:1
        Xmatrix{w} = Xmatrix{1}(:,w);
        Ymatrix{w} = Ymatrix{1}(:,w);
        % Not naming the curves in this case, can do later.
        varnames{w} = varnames{1};
    end
end

if(~ishold)
    washeld = false;
    %axes1 = axes;
    axes1 = gca;
    cla(gca);
    %box(axes1,'on');
    hold(axes1,'all');
    xlabel(PlotLabel(xvect));
    ylabel(PlotLabel(yvect));
    %title(strrep(searchString,'_',' '));
    title(strrep(inputname(1), '_', ' '));
else
    washeld = true;
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
    [xplot, order] = sort(cell2mat(Xmatrix));
    yplot = cell2mat(Ymatrix);
    yplot = yplot(order);
    plot(xplot, yplot, '.-','DisplayName',searchString,'Color',PlotColor(1,:));
else
    for i = 1:nplots
        if ncolors > 1
            PlotColor = [interpR(i) interpG(i) interpB(i)];
        end
        % Decide display name
        if ~Direct && isbasefield(varnames{i},'DisplayName')
            % DisplayName needs to be a string
            DisplayName = evalin('base',[varnames{i} '.' DisplayName]);
        elseif ~Direct
            DisplayName = strrep(varnames{i},'_','');
        elseif isfield(searchString,'Info') && isfield(searchString.Info,'DiskID')
            DisplayName = strrep(searchString.Info.DiskID,'_','');
        end
        plot(Xmatrix{i},Ymatrix{i},'.-','DisplayName',DisplayName,'Color',PlotColor,'LineWidth',2);
    end
end

%legend('toggle');
if(~washeld)
    hold('off');
end

end
