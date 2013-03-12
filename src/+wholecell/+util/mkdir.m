%mkdir
% Creates directory as well as parent directories if necessary. Similar to
% linux mkdir -p.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
function mkdir(path, recursive)

%make parent directories
if nargin >= 2 && recursive
    path = absolutepath(path, '', false);
    idxs = find(path == filesep);
    for i = 1:numel(idxs)
        if ~exist(path(1:idxs(i) - 1), 'dir')
            mkdir(path(1:idxs(i) - 1))
        end
    end
end

%make directory
if ~exist(path, 'dir')
    mkdir(path);
end