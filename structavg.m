function avgstruct = structavg(searchstring)
% For all structs in the workspace whose name contains searchstring,
% average all double array fields that the structs have in common

% could extend to average cell arrays of double arrays ..

matchingstructs = FindBaseVars(searchstring,'struct');
if isempty(matchingstructs)
    error('No fields match your search string')
end


%% Find matrix fields that structs have in common
structi{1} = evalin('base',matchingstructs{1});
CommonFields = fieldnames(structi{1});
for i = 2:length(matchingstructs)
    %save cell array of matching structs
    structi{i} = evalin('base',matchingstructs{i});
    %save list of field names for each struct
    fields{i} = fieldnames(structi{i});
    %find fields that are in every struct
    CommonFields = intersect(CommonFields,fields{i});
end

for i = 1:length(matchingstructs)
    for j = 1:length(CommonFields)
        %save each common field class
        Class{i,j} = class(structi{i}.(CommonFields{j}));
        %and size
        Size{i,j} = size(structi{i}.(CommonFields{j}));
    end
end

for j= length(CommonFields):-1:1
    for i = 1:length(matchingstructs)
        if ~strcmp(Class{i,j},'double');
            %If not double array, don't average
            CommonFields(j) = [];
            break;
        end    
    end
end

% for now, will error if fields are double arrays with different sizes
avgstruct = struct();
for j = 1:length(CommonFields)    
    avgstruct = setfield(avgstruct,CommonFields{j},0);
    for i = 1:length(matchingstructs)
        avgstruct.(CommonFields{j}) = avgstruct.(CommonFields{j}) + structi{i}.(CommonFields{j})/length(matchingstructs);
    end
end

end
