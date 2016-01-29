function uiimport(filepath)
% UIIMPORT overloaded for custom files.
%
% Executed when file is dragged to the workspace, will call custom import
% routines for recognized files, and import unrecognized files with
% standard uiimport.  A directory can also be dragged into the workspace,
% in which case all the recognized file types contained in the first level
% will be imported and unrecognized file types ignored.
%
% T Hennen 2013

%Define DMS types
DMStypes = {'.VHD','.VRD','.VDD'};
%Define PK types
PKtypes = {'.mdd','.svd','.svp'};

[dirpath, filename, extension] = fileparts(filepath);

%If the dragged file has an extension (not a folder)
if ~isempty(extension)
    Files(1).name = [filename extension];
    % Determine if file extension matches a DMS type
    isDMS = any(strcmpi(DMStypes,extension));
    % or PK type
    isPK = any(strcmpi(PKtypes,extension));
    
    % If other file type, run uiimport as usual
    if ~isDMS && ~isPK
        presentPWD = pwd;
        cd([matlabroot '\toolbox\matlab\codetools\']);
        uiimport(filepath);
        cd(presentPWD);
    end
else
    % User dragged in a folder, search for and import DMS and PK types in that folder, 
    % ignoring other types
    Files = dir(filepath);
    dirpath = [dirpath '\' filename];
end


for i=1:length(Files)
    [a, filename, extension] = fileparts(Files(i).name);
    % Determine if file extension matches a DMS type
    isDMS = any(strcmpi(DMStypes,extension));
    % or PK type
    isPK = any(strcmpi(PKtypes,extension));
    
    if (isDMS)
        % All DMS data uses same import function
        DMSImport([dirpath '\' Files(i).name]);
    elseif (isPK)
        % Determine the type of PK data by looking at extension and the end of the file name
        if strcmpi(extension,'.mdd')
            
            % If it's a recoil or recoil+dyn file, our only hope of identifying it is to 
            % find the word recoil in the filename.
            if regexpi(filename,'recoil')
                % Recoil files end in 0_0
                if strcmp(filename(end-2:end),'0_0')
                    PKImportRecoil([dirpath '\' Files(i).name]);
                    % Dynamic files end in whatever they please (eg. 0_108, 0_77)
                else
                    PKImportDyn([dirpath '\' Files(i).name]);
                end
            else
                % General PK import for .mdd with no "recoil" in the name
                PKImport([dirpath '\' Files(i).name]);
            end
            
        elseif strcmpi(extension,'.svd')
            PKImportSVD([dirpath '\' Files(i).name]);
        elseif strcmpi(extension,'.svp')
            PKImportSVP([dirpath '\' Files(i).name]);
        end
    end
end
