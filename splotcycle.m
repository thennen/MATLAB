function [varnames, inputs] = splotcycle(searchstring, x, y)
% Cycle through variables, calling splot on each.  Values entered into
% prompt can be output.  Prompt keywords are
% 
% b           move back one plot
% !command    execute command
% s #         skips # plots
% q           quit splotting

    varnames = findbasevars(searchstring, 'struct');
    inputs = cell(length(varnames), 1);
    i = 1;
    
    fig = figure;
    xlabel(PlotLabel(x));
    ylabel(PlotLabel(y));
    while i <= length(varnames)        
        %delete(gca)
        
        items = get(gca, 'Children');
        if ~isempty(items)
            delete(items(:));
        end
        hold all
        
        evalin('base',['splot(''' varnames{i} ''', ''' x ''', ''' y ''')'])
        title(strrep(varnames{i},'_',' '))
        %xlabel('Field (Oe)')
        %ylabel('Kerr')
        stringin = input('do something: ','s');
        inputs{i} = stringin;
        
        while strncmp(inputs{i},'!',1)
            eval(stringin(2:end))
            %get another input
            stringin = input('do something: ','s');
            inputs{i} = stringin;
        end
        
        if strcmp(inputs{i},'b')
            %move back one plot
            i = i - 1;
            continue;
        end
        
        if strncmp(inputs{i},'s',1)
            %skip a number of plots specified by number after s
            % ex: 's 8' skips 8
            s = str2num(stringin(2:end));
            i = i + s;
            continue;
        end
            
        if strcmp(inputs{i},'q')
            break;
        end
        i = i + 1;
    end
    %close
end