function [multiplier, inputs] = visualsubtraction(x1, y1, x2, y2, multiplier, varargin)

    if nargin == 6
        % extra argument gives prompt message
        prompt = varargin{1};
    else
        prompt = 'Do Something: ';
    end

    inputs = cell(length(multiplier), 1);
    i = 1;

    if ~ishold
        figure;
    end
    
    %if ~all(x1 == x2) % all X must match
    if length(x1) == length(x2) && ~all(abs(x1 - x2) < 5)  % All x must be within 5 units
        y2interp = interp1(x2, y2, x1, 'linear', 'extrap'); % error if not strictly increasing..
    else
        y2interp = y2;
    end        

    ploti = [];
    while i <= length(multiplier)
        % Should delete previous lines from the axis, leaving others
        delete(ploti)
        ploti = plot(x1, y1-multiplier(i).*y2interp, '-.');
        title(num2str(multiplier(i)))
        xlabel('X')
        %ylabel('y1 - m y2')
        stringin = input([prompt ' '],'s');
        inputs{i} = stringin;
        
        if strncmp(inputs{i},'!',1)
            try
                evalin('base',inputs{i}(2:end));
            catch
                disp('didn''t work')
            end
            continue;
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
            for k = i:length(multiplier)
                inputs{k} = ''; %match output size
            end
            break;
        end
        i = i + 1;
    end
    %close
    
    % output multiplier corresponding to 'g' if there's only one
    if sum(strcmpi(inputs, 'g')) == 1
        multiplier = multiplier(strcmpi(inputs, 'g'));
    end
    
end
