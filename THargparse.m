function [ArgStructOut, ParamsSet] = THargparse(argin, defaults)
% THARGPARSE parses input arguments in Tyler Hennen fashion
% Parses arguments in a simple way that does not require values
% for each parameter, assuming a value of true by default.
% Argument values cannot be strings, all parameter names will 
% be converted to lowercase.
% No value checking of any kind

if ~isstruct(defaults)
    defaults = struct();
end

ArgStructOut = defaults;
% Keep track of the parameters that are set by THargparse
ParamsSet = {};

nargs = length(argin);

for i=1:nargs
    argi = argin{i};
    if ischar(argi)
        argi = lower(argi);
        ParamsSet{end + 1} = argi;
        % If the following argument is not a string
        % then set it as the parameter value
        if i~=nargs && ~ischar(argin{i+1})
            ArgStructOut.(argi) = argin{i+1};
        else
            ArgStructOut.(argi) = true;
        end
    end
end

end
