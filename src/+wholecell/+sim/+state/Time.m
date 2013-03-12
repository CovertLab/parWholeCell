%Time
% Time state variable. Represents the current time lapsed since the start
% of the simulation in seconds.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
classdef Time < wholecell.sim.state.State
    %metadata
    properties
        meta = wholecell.util.struct(...
            'id', 'Time', ...
            'name', 'Time', ...
            'dynamics', {'value'}, ...
            'units', struct('value', 's') ...
            )
    end
    
    %mass
    properties
        value %s
    end
    
    methods
        %constructor
        function this = Time(varargin)
            this = this@wholecell.sim.state.State(varargin{:});
        end
        
        %allocate memory
        function this = allocate(this)
            this.allocate@wholecell.sim.state.State();
            
            this.value = zeros(1, 1);
        end
        
        %calculate initial conditions
        function this = calcInitialConditions(this)
            this.value = 0;
        end
    end
end