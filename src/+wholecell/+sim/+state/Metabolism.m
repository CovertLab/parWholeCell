%Metabolism
% Metabolism state variable. Represents the instantaneous growth rate
% (fg/h) and metabolic reaction fluxes (reactions/s).
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
classdef Metabolism < wholecell.sim.state.State
    %metadata
    properties
        meta = wholecell.util.struct(...
            'id', 'Metabolism', ...
            'name', 'Metabolism', ...
            'dynamics', {'growth'; 'fluxes'}, ...
            'units', struct(...
                'growth', 'fg/h', ...
                'fluxes', 'reactions/s' ...
                ) ...
            )
        
        reactionIds
        reactionNames
        reactionIdx
    end
    
    %state, process references
    properties
        moleculeCounts %moleculeCounts states
        metabolism     %metabolism process
    end
    
    %mass
    properties
        growth %fg/h
        fluxes %reactions/s
    end
    
    methods
        %constructor
        function this = Metabolism(varargin)
            this = this@wholecell.sim.state.State(varargin{:});
        end
        
        %construct state-process graph
        function this = initialize(this, sim, kb)
            this.initialize@wholecell.sim.state.State(sim, kb);
            
            this.moleculeCounts = sim.getState('MoleculeCounts');
            this.metabolism = sim.getProcess('Metabolism');
            
            this.reactionIds = {kb.reactions.id}';
            this.reactionNames = {kb.reactions.name}';
        end
        
        %allocate memory
        function this = allocate(this)
            this.allocate@wholecell.sim.state.State();
            
            this.growth = zeros(1, 1);
            this.fluxes = zeros(numel(this.reactionIds), 1);
        end
        
        %calculate initial conditions
        function this = calcInitialConditions(this)
            mc = this.moleculeCounts;
            met = this.metabolism;
            
            bounds = this.metabolism.calcFluxBounds(...
                mc.counts(met.metabolite.mapping), mc.counts(met.enzyme.mapping));
            [this.growth, this.fluxes] = met.calcGrowthRate(bounds);
        end
    end
end