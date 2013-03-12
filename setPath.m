%setPath
% Adds folders to MATLAB path and MATLAB Java path.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Created: 3/1/2013
function setPath()
if isdeployed
    return;
end

%% MATLAB path

%libraries
addpath(fullfile(pwd, 'lib/absolutepath'));              %full file name                  http://www.mathworks.com/matlabcentral/fileexchange/3857
addpath(fullfile(pwd, 'lib/glpkmex-2.9'));               %linear programming solver       http://glpkmex.sourceforge.net/ (Note win64 version is actually v2.7; maci version is actually v2.8)
addpath(fullfile(pwd, 'lib/matlab_xunit_3.0.1/xunit'));  %unit testing                    http://www.mathworks.com/matlabcentral/fileexchange/22846
addpath(fullfile(pwd, 'lib/util/strutil'));              %string utility functions        http://home.online.no/~pjacklam/matlab/software/util/
addpath(fullfile(pwd, 'lib/util/matutil'));              %matrix utility functions        http://home.online.no/~pjacklam/matlab/software/util/

%source code
addpath(fullfile(pwd, 'src_test'));
addpath(fullfile(pwd, 'src'));
addpath(pwd);
