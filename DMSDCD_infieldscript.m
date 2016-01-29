% notably awful script for correcting IN FIELD dcd


vars = findbasevars('4_DCD','struct');
for i=1:length(vars)
    DCDvar = evalin('base',vars{i});
    DCDvar = findperppar(DCDvar);
    
    DCDvar.Perpinfield = DCDvar.Perpinfield *1000000;
    
    uppermask = DCDvar.Field > 400;
    lowermask = DCDvar.Field < 25;
    
    upperfit = polyfit(DCDvar.Field(uppermask),DCDvar.Perpinfield(uppermask),1);
    lowerfit = polyfit(DCDvar.Field(lowermask),DCDvar.Perpinfield(lowermask),1);
    
    avgslope = (upperfit(1) + lowerfit(1))/2;
    
    offset = (upperfit(2) + lowerfit(2))/2;
    
    
    DCDvar.infieldMs = ( upperfit(2)-lowerfit(2) )/2;
    DCDvar.infieldlowersat = lowerfit(2);
    DCDvar.infieldoffset = offset;
    
    lowersati(i) = lowerfit(2);
    offseti(i) = offset;
    
    %avgoffset = -0.413;
    %avglowersat = -197.795 -avgoffset;
    
    %DCDvar.CorrPerpinfield = DCDvar.Perpinfield -avgslope.*DCDvar.Field - avgoffset;
    DCDvar.CorrPerpinfield = DCDvar.Perpinfield -upperfit(1).*DCDvar.Field - avgoffset;
    
    slope = (mean(DCDvar.CorrPerpinfield(end-1:end)) + avglowersat)/DCDvar.Field(end);
    
    DCDvar.CorrPerpinfield = DCDvar.CorrPerpinfield - slope.*DCDvar.Field;
    
    
    
    assignin('base',vars{i},DCDvar);
    
    avgoffset = mean(offseti);  %run twice
    avglowersat = mean(lowersati) - avgoffset;
    
    
    
end