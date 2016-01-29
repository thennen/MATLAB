function FitRange = FindFitRange(DataIn)
% Searches DataIn fields and base workspace for a matrix FitRange =
% [a b] OR [a b c d], and if found returns its value, otherwise returns the
% default value [0.75*maxfield, maxfield].  DataIn fields have priority
% over base workspace variable
FitRangeField = [];
FitRangeBase =[];
FitRange = [];
% Look in fields of DataIn
fields = fieldnames(DataIn);
fitrangematch = fields(strcmpi(fields,'fitrange'));
if ~isempty(fitrangematch)
    FitRangeField = DataIn.(fitrangematch{1});
    if length(fitrangematch) > 1
        warning(['Multiple fitrange fields, using' fitrangematch{1} '=' mat2str(DataIn.(fitrangematch{1}))])
    end
end
% Look in base workspace
string = 'who (''-regexp'' ,''(?i)fitrange'')';
FitRangeVar = evalin('base', string);
if ~isempty(FitRangeVar)
    FitRangeBase = evalin('base', FitRangeVar{1});
end
% Check format of fitrange(s), use default if no valid alternatives,
% warn of ambiguities
if ~isempty(FitRangeField)
    ismat = ismatrix(FitRangeField);
    len = length(FitRangeField);
    if ismat && len == 2 && FitRangeField(1) < FitRangeField(2)
        FitRange = FitRangeField;
    elseif ismat && len == 4 && all(diff(FitRangeField) > 0)
        FitRange = FitRangeField;
    else
        warning('Fitrange field has wrong format and will be ignored.')
    end
end
if ~isempty(FitRangeBase) && ~isempty(FitRange) && ~isempty(FitRangeField)
    disp('Found FitRange variable in struct and in base workspace.  Using Struct value.')
elseif ~isempty(FitRangeBase)
    ismat = ismatrix(FitRangeBase);
    len = length(FitRangeBase);
    if ismat && len == 2 && FitRangeBase(1) < FitRangeBase(2)
        FitRange = FitRangeBase;
    elseif ismat && len == 4 && all(diff(FitRangeBase) > 0)
        FitRange = FitRangeBase;
    else
        warning('Fitrange base variable has wrong format and will be ignored.')
    end
end
% If FitRange still hasn't been assigned, give it default value
if isempty(FitRange)
    if iscell(DataIn.Field)
        MaxField = max(cell2mat(DataIn.Field));
    elseif ismatrix(DataIn.Field)
        MaxField = max(max(DataIn.Field));
    else
        error('Something is wrong with the ''Field'' field')
    end
    FitRange = [0.75*MaxField MaxField];
end
end
