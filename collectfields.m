function Result = collectfields(varin,varargin)
% COLLECTFIELDS Collect fields from structs determined by varin
% e.g. PKParams = CollectFields('PKRecoil','Hc','Hn','Hs')

if ischar(varin)
    varnames = findbasevars(varin,'struct');
    % build table of variable names matching search string
    Result = varnames;
    % shift down 1 cell
    Result(2:end+1,1) = Result(1:end);                          
    Result{1} = 'Name';
    
elseif iscell(varin)                                            
    %varnames = cell(size(varin,1),1);
    %varcounter = 0;
    if ~isempty(varin{1}) && ~isempty(findbasevars(varin{1}))
        %if varin(1,1) matches a workspace variable, there's probably no header.  make one.
       Result = varin;
       Result(2:end+1,1:end) = varin;
       Result(1,:) = {''};
    else
        Result = varin; 
    end
    
else
    error('varin must be a string or cell array')
end

StartSize = size(Result,2);

% Loop through columns, skipping any that existed in the input cell array
for i=StartSize+1:StartSize+nargin-1
    % Loop through rows
    for j=1:size(Result,1)
        MatchingVars = findbasevars(Result{j});
        % Fill first row with requested field names (labels)
        if j==1
            Result{j,i} = varargin{i-StartSize};
        elseif ~isempty(Result{j}) && ~isempty(MatchingVars)
            % support for nested structures
            dotlocation = strfind(varargin{i-StartSize},'.');                                
            % if requested field has a '.' in it, assume everything before
            % the '.' is a nested struct
            if isempty(dotlocation)
                structname = MatchingVars{1};
                fieldname = varargin{i-StartSize}(1:end);
            else
                structname = [MatchingVars{1} '.' varargin{i-StartSize}(1:dotlocation(end)-1)];
                fieldname = varargin{i-StartSize}(dotlocation(end)+1:end);
            end  
            
            bracketlocation = strfind(fieldname,'(');
            if isempty(bracketlocation)
                bracketlocation = length(fieldname)+1;
            end
                                                                    
            % If the requested field exists..
            if isbasefield(structname,fieldname(1:bracketlocation-1))
            Result{j,i} = evalin('base',[structname '.' fieldname]);
            end
        end
   end
end
