function pkrecoilrecipe(filename,MaxField,SweepRate,RFields)
% Generates a recoil macro recipe for the polar kerr tool
%   Recipe will perform minor loops, in order, for the reverse fields specified in
%   the input double array 'RFields'.  RFields may be positive or negative,
%   but -abs(RFields) will be used (no support yet for positive reverse
%   fields).  
%
%   Example 1 (Maximum field 15000 Oe, 5000 Oe/s sweep rate, standard reverse fields):  
%   pkrecoilrecipe('G:\NewRecoilRecipe.man',15000,5000,[0:240:12000])
%
%   Example 2 (finer reverse fields, with major loop at the end)
%   pkrecoilrecipe('G:\NewRecoilRecipe.man',15000,5000,[0:120:12000 15000])
%

% Calculate total number of lines
% nlines = 4 + 2*length(RFields) + 1;  % 4 lines for header/major loop, 2 lines per minor loop, plus one line for no apparent reason
% 5 lines for header/major loop, 2 lines per minor loop
nlines = 5 + 2*length(RFields);
% Create and open a new file, destroying existing file
fhandle = fopen(filename,'w');

% Field step between samples
Hstep = 150;
% Number of samples to average
Avgs = 10;
AvgTime = Hstep/SweepRate;

%Header and major loop
Header = {
    '[File Summary]';
    'MacroDataFileName=C:\HeadWaferMappingSystem\Data\Sylvia\Recoil-2T.mdd';
    'Version=13.5.0';
    'DoNotCreateMacroDataFile=False'
    ['NumOfSeqEntryLines=' num2str(nlines)];
    '[Sequence Entries]'
    'Line 0=Create New FileName For Data;False;No File;*';
    'Line 1=Set X Field Range;30K;*'
    'Line 2=Wait;60;*'
    ['Line 3=Set X Field;' num2str(MaxField/SweepRate) ';' num2str(MaxField) ';*']
    ['Line 4=X Field Sweep and Measure MOKE Signal;' num2str(2*MaxField/SweepRate) ';' num2str(-MaxField) ';;True;' num2str(AvgTime) ';' num2str(Avgs) ';*'];
    ['Line 5=X Field Sweep and Measure MOKE Signal;' num2str(2*MaxField/SweepRate) ';' num2str(MaxField)  ';;True;' num2str(AvgTime) ';' num2str(Avgs) ';*'];
    };
% Write header to file
for n = 1:length(Header);
    fprintf(fhandle,'%s\r\n',Header{n});
end

RFields = -abs(RFields);
% Write Minor loops
LineNumber = 6;
for i = 1:length(RFields)
    SweepTimei = (-RFields(i)+MaxField)/SweepRate;
    %Set Reverse Field
    fprintf(fhandle,'%s %u=Set X Field;%.3f;%.3f;*\r\n','Line',LineNumber,SweepTimei,RFields(i));
    LineNumber = LineNumber + 1;
    % Sweep and measure to MaxField
    fprintf(fhandle,'%s %u=X Field Sweep and Measure MOKE Signal;%.3f;%.3f;;True;%.4f;%u;*\r\n','Line',LineNumber,SweepTimei,MaxField,AvgTime,Avgs);
    LineNumber = LineNumber + 1;
end

% Go to MaxField for no reason (for consistency with existing recipe)
% fprintf(fhandle,'Line %d=Set X Field;4;%d;*',LineNumber,MaxField);  

fclose(fhandle);

end

