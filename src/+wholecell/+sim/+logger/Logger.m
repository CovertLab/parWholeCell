%Logger
% Abstract class which defines the interface loggers expose to the
% simulation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
classdef Logger < handle
    methods (Abstract = true)
        %inialize -- called at beginning of simulation
        this = initialize(this, sim)
        
        %append -- called at each iteration of simulation
        this = append(this, sim)
        
        %finalize -- called at end of simulation
        this = finalize(this, sim)
    end
end