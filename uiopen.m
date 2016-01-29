function uiopen(filepath,direct)
% UIOPEN overloaded for custom Files.
% Executed when file is dragged into Command Window

DMStypes = {'.VHD','.VRD','.VDD'};                     %Define DMS types
PKtypes = {'.mdd','.svd'};                             %Define PK types

[~, filename, extension] = fileparts(filepath);
isDMS = any(strcmpi(DMStypes,extension));               %Determine if file extension matches a DMS type
isPK = any(strcmpi(PKtypes,extension));                 %or PK type

if (isDMS && (direct))
    varname = DMSImport(filepath);
    evalin('base',['dms(''' varname ''')']);
elseif (isPK && (direct))
    if strcmpi(extension,'.mdd')                                   % Determine the type of PK data by looking at extension and the end of the file name
        if strcmp(filename(end-2:end),'108')                % Dynamic files end in _108
            PKImportDyn(filepath);                          %
        elseif strcmp(filename(end-1:end),'_0')             % Recoil files end in _0
            PKImportRecoil(filepath);                       %
        end
    elseif strcmpi(extension,'.svd')
        PKImportSVD(filepath);
    end
else
    presentPWD = pwd;
    cd([matlabroot '\toolbox\matlab\uitools']);
    uiopen(filepath,direct);                        %If other file type, run uiopen as usual
    cd(presentPWD);
end