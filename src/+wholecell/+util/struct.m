%struct
% Creates a scruct. Similar to built-in method, but creates a 1-dimensional
% struct even if the values are multi-dimensional cell arrays.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
function val = struct(varargin)
val = struct();
for i = 1:2:nargin
    val.(varargin{i}) = varargin{i+1};
end