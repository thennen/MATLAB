function pkhenkelrecipe(filename,FieldsArray)
%   Generates a macro recipe for the polar kerr tool
%
%   Recipe will perform virgin loops (alternates between Happ and 0 field) followed by recoil loops
%   (after saturation, alternates between -Happ and 0). One major loop is added in the end.
%   FieldsArray may be positive or negative, but abs(FieldsArray) will be used
%   for virgin loops and -abs(FieldsArray) for recoil (no support for positive reverse fields).  
%   
%
%

MaxField = 15000;
SweepTime = 6;                       % Time to sweep from positive to negative MaxField
SweepRate = 2*MaxField/SweepTime;
% Calculate total number of lines
nlines = 2 + 4*length(FieldsArray) + 5;  % two lines for header (set field range and set field), two lines per virgin loop,
                                         % one line for field saturation, two lines per minor loop, three lines for two
                                         % major loops and finally one line for no apparent reason

fhandle = fopen(filename,'w');   % Create and open a new file, destroying existing file


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
    'Line 2=Set X Field;5;0;*'
    };
% Write header to file
for n = 1:length(Header);
    fprintf(fhandle,'%s\r\n',Header{n});
end

% Write Virgin loops
LineNumber = 3;
for i = 1:length(FieldsArray)
    SweepTime = abs(FieldsArray(i))/SweepRate;
    fprintf(fhandle,'%s %u=X Field Sweep and Measure MOKE Signal;%.3f;%d;;True;0.03;10;*\r\n','Line',LineNumber,SweepTime,abs(FieldsArray(i)));  % Sweep and measure up to Ha, the applied field
    LineNumber = LineNumber + 1;
    fprintf(fhandle,'%s %u=X Field Sweep and Measure MOKE Signal;%.3f;%d;;True;0.03;10;*\r\n','Line',LineNumber,SweepTime,0); % Sweep and measure to 0
    LineNumber = LineNumber + 1;
end

fprintf(fhandle,'%s %u=Set X Field;5;%d;*\r\n','Line',LineNumber,MaxField); % saturate field
LineNumber = LineNumber + 1;

% Write Minor loops
for i = 1:length(FieldsArray)
    SweepTime = abs(FieldsArray(i))/SweepRate;
    if i == 1
        SweepTime = (abs(FieldsArray(i))+MaxField)/SweepRate;
    end    
    fprintf(fhandle,'%s %u=X Field Sweep and Measure MOKE Signal;%.3f;%d;;True;0.03;10;*\r\n','Line',LineNumber,SweepTime,-abs(FieldsArray(i)));  %Set Reverse Field
    LineNumber = LineNumber + 1;
    if i == 1
        SweepTime = abs(FieldsArray(i))/SweepRate; 
    end    
    fprintf(fhandle,'%s %u=X Field Sweep and Measure MOKE Signal;%.3f;%d;;True;0.03;10;*\r\n','Line',LineNumber,SweepTime,0); % Sweep and measure to 0
    LineNumber = LineNumber + 1;
end

% Write Major loop
fprintf(fhandle,'%s %u=Set X Field;5;%d;*\r\n','Line',LineNumber,MaxField); % saturate field
fprintf(fhandle,'%s %u=X Field Sweep and Measure MOKE Signal;%.3f;%d;;True;0.03;10;*\r\n','Line',LineNumber+1,SweepTime,-MaxField); % Sweep and measure to -MaxField
fprintf(fhandle,'%s %u=X Field Sweep and Measure MOKE Signal;%.3f;%d;;True;0.03;10;*\r\n','Line',LineNumber+2,SweepTime,MaxField); % Sweep and measure to MaxField

fprintf(fhandle,'Line %d=Set X Field;4;%d;*',LineNumber+3,MaxField);  % Go to MaxField for no reason (for consistency with existing recipe)

fclose(fhandle);

end
