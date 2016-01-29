function structplot2(searchString,x,y,Param)
% 

if(~ishold)
    washeld = false;
    figure1 = figure;
    axes1 = axes('Parent',figure1);
    %axes1 = axes;
    box(axes1,'on');
    hold(axes1,'all');
    xlabel(x);
    ylabel(y);
    legend();
    title(strrep(searchString,'_',' '));
else
    washeld = true;
end

StartColor = [0 0 0];
EndColor = [0 .3 1];

varnames = FindBaseVars(searchString,'struct');
numvars = length(varnames);

XArray = zeros(numvars,1);
YArray = zeros(numvars,1);
ParamArray = zeros(numvars,1);
for i = numvars:-1:1
    if isbasefield(varnames{i},x) && isbasefield(varnames{i},y) && isbasefield(varnames{i},Param)
        Xi = evalin('base',[varnames{i} '.' x]);
        Yi = evalin('base',[varnames{i} '.' y]);
        Parami = evalin('base',[varnames{i} '.' Param]);
        if length(Xi) == 1 && length(Yi) == 1 && length(Parami) == 1
            XArray(i) = Xi;
            YArray(i) = Yi;
            ParamArray(i) = Parami;
        else
            disp(['Specified field in ' varnames{i} ' has length > 1']);
            XArray(i) = [];
            YArray(i) = [];
            ParamArray(i) = [];
        end
    else
        disp(['Specified field in ' varnames{i} ' does not exist']);
        XArray(i) = [];
        YArray(i) = [];
        ParamArray(i) = [];
    end
end

PlotGroups = unique(ParamArray);
NGroups = length(PlotGroups);

for j=1:NGroups
    XGroup{j} = XArray(ParamArray == PlotGroups(j));
    YGroup{j} = YArray(ParamArray == PlotGroups(j));
    
    % Sort by ascending X
    [XGroup{j},Orderj] = sort(XGroup{j});
    YGroup{j} = YGroup{j}(Orderj);
    
    PlotColor = (NGroups-j).*StartColor./NGroups + j.*EndColor./NGroups;
    DisplayName = num2str(PlotGroups(j));
    
    plot(XGroup{j},YGroup{j},'-s','DisplayName',DisplayName,'MarkerSize',5,'MarkerFaceColor',PlotColor,'Color',PlotColor,'LineWidth',2);
        
end

legend('show');

if(~washeld)
    hold(axes1,'off');
end

end