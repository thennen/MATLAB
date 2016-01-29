function pksetup(filename,RFields)
% Generates a recoil macro recipe for the polar kerr tool
%   Recipe will perform minor loops, in order, for the reverse fields specified in
%   the input double array 'RFields'.  RFields may be positive or negative,
%   but -abs(RFields) will be used (no support for positive reverse
%   fields).  
%
%   Example 1 (standard reverse fields):  
%   pkrecoilrecipe('G:\NewRecoilRecipe.man',[0:240:12000])
%   Example 2 (finer reverse fields, with major loop at the end)
%   pkrecoilrecipe('G:\NewRecoilRecipe.man',[0:120:12000 15000])
%

MaxField = 15000;
SweepTime = 6;                       %Time to sweep from positive to negative MaxField
SweepRate = 2*MaxField/SweepTime;
% Calculate total number of lines
nlines = 4 + 2*length(RFields) + 1;  % 4 lines for header/major loop, Two lines per minor loop, plus one line for no apparent reason

fhandle = fopen(filename,'w');       % Create and open a new file, destroying existing file


%Header and major loop
Header = {
    '[File Summary]';
    'MacroDataFileName=C:\HeadWaferMappingSystem\Data\THenn\THRecoil-12T.mdd';
    'Version=13.5.0';
    'DoNotCreateMacroDataFile=False'
    ['NumOfSeqEntryLines=' num2str(nlines)];
    '[Sequence Entries]'
    'Line 0=Create New FileName For Data;False;No File;*';
    'Line 1=Set X Field Range;30K;*'
    ['Line 2=Set X Field;5;' num2str(MaxField) ';*']
    ['Line 3=X Field Sweep and Measure MOKE Signal;' num2str(SweepTime) ';' num2str(-MaxField) ';;True;0.03;10;*'];
    ['Line 4=X Field Sweep and Measure MOKE Signal;' num2str(SweepTime) ';' num2str(MaxField)  ';;True;0.03;10;*'];
    };
% Write header to file
for n = 1:length(Header);
    fprintf(fhandle,'%s\r\n',Header{n});
end


% Write Minor loops
LineNumber = 5;
for i = 1:length(RFields)
    SweepTime = (RFields(i)+MaxField)/SweepRate;
    fprintf(fhandle,'%s %u=Set X Field;%.3f;%d;*\r\n','Line',LineNumber,SweepTime,-abs(RFields(i)));  %Set Reverse Field
    LineNumber = LineNumber + 1;
    fprintf(fhandle,'%s %u=X Field Sweep and Measure MOKE Signal;%.3f;%d;;True;0.03;10;*\r\n','Line',LineNumber,SweepTime,MaxField); % Sweep and measure to MaxField
    LineNumber = LineNumber + 1;
end

fprintf(fhandle,'Line %d=Set X Field;4;%d;*',LineNumber,MaxField);  % Go to MaxField for no reason (for consistency with existing recipe)

fclose(fhandle);

end

