function funnames = DepSubFun(CS,ExcludeRoot);
% DEPSUBFUN  Locate dependent functions of an M-file.
% But, unlike DEPFUN, it will look *inside* m-files for
% *subfunctiosn* and the ability to excludes files that
% are in the MATLAB root (default).
%
% The intent of finding the subfunctions is to prevent
% repetition and to identify which functions need to be
% updated so if you fix one copy, then you will know there
% are other copies out there that need to be fixed.
%
% The intent of excluding files in the MATLAB root
% is to find which files must be included in a distribution
% so they may be uploaded to the MATLAB File Exchange.
%
% INPUTS:
% CS, a cellstring of file names (the files you want analysed)
% ExcludeRoot, 0 if you want all function names, 1 if you only
%    want the ones not in the MATLAB root.
%
% OUTPUTS
% funnames, a cellstring of functions and subfunctions
%
% USAGE
% >>d=dir('*.m');
% >> DepSubFun({d.name})
%
% See DEPFUN
%
% Keyworkds depfun, dependent, function, subfunction, name
%
% It's not fancy, but it works

% REVISION HISTORY
% Version 1.00
% 4/13/05 [/] Created
% 4/14/05 [/] Nested for loop variable name changed

% Michael Robbins
% robbins@bloomberg.net
% MichaelRobbins1@yahoo.com
% MichaelRobbinsUsenet@yahoo.com

% CONSTANTS
TRUE      = 1;
FALSE     = 0;
EMPTYCELL = {};

% INITIALIZE
if nargin<2 ExcludeRoot=TRUE; end;
funnames  = EMPTYCELL;
TESTREGEX = FALSE;

% FOR EACH FILE
for i=1:length(CS)
    
    % GET DEPENDENT FUNCTIONS' FILENAMES
    f = depfun(CS{i},'-quiet');
    
    % ADD FUNCTIONS THAT ARE NOT IN MATALBROOT
    if ExcludeRoot
        foofiles = f(setdiff([1:length(f)],strmatch(matlabroot,f)));
    else
        foofiles = f;
    end;
    
    % LOOP THROUGH EACH FOOFILE
    if ~isempty(foofiles)
        for j=1:length(foofiles)
            FID = fopen(foofiles{j},'r');
            if FID ~=-1
                foo = GetSubFunNames(FID);
                funnames={funnames{:} foo{:}};
            else
                warning('DepFunInFile :: Bad FID');
            end; % IF FID
        end; % FOR
    else
        warning('DepFunInFile :: foofiles is empty');
    end; % IF ~ISEMPTY
end;