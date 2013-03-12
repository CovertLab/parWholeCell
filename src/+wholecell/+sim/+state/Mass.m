%Mass
% Mass state variable. Represents the total cellular mass.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
classdef Mass < wholecell.sim.state.State
    %metadata
    properties
        meta = wholecell.util.struct(...
            'id', 'Mass', ...
            'name', 'Mass', ...
            'dynamics', {'total'; 'cell'; 'cellDry'; 'metabolite'; 'rna'; 'protein'}, ...
            'units', struct(...
                'total', 'fg', ...
                'cell', 'fg', ...
                'cellDry', 'fg', ...
                'metabolite', 'fg', ...
                'rna', 'fg', ...
                'protein', 'fg' ...
                ) ...
            )
        
        compartments = [ %compartment ids, names
            struct('id', 'c', 'name', 'Cytosol')
            struct('id', 'e', 'name', 'Extracellular space')
            struct('id', 'm', 'name', 'Membrane')
            ];
        cIdx = struct('c', 1, 'e', 2, 'm', 3);
    end
    
    %references to other states
    properties
        moleculeCounts
    end
    
    %mass (fg)
    properties
        total
        cell
        cellDry
        metabolite
        rna
        protein
    end
    
    methods
        %constructor
        function this = Mass(varargin)
            this = this@wholecell.sim.state.State(varargin{:});
        end
        
        %construct object graph
        function this = initialize(this, sim, kb)
            this.initialize@wholecell.sim.state.State(sim, kb);
            
            this.moleculeCounts = sim.getState('MoleculeCounts');
        end
        
        %allocate memory
        function this = allocate(this)
            this.allocate@wholecell.sim.state.State();
            
            this.total = zeros(1, numel(this.compartments));
            this.cell = zeros(1, numel(this.compartments));
            this.cellDry = zeros(1, numel(this.compartments));
            this.metabolite = zeros(1, numel(this.compartments));
            this.rna = zeros(1, numel(this.compartments));
            this.protein = zeros(1, numel(this.compartments));
        end
        
        %calculate (and cache) any dependent properties
        function this = calculate(this)
            import wholecell.util.Constants;
            
            mc = this.moleculeCounts;
            
            %total
            this.total = (...
                + mc.mws' * mc.counts  ...
                ) / Constants.nAvogadro * 1e15;
            
            %cell
            this.metabolite = mc.mws(mc.types == mc.typeVals.metabolite, 1)' * mc.counts(mc.types == mc.typeVals.metabolite, :) / Constants.nAvogadro * 1e15;
            this.rna        = mc.mws(mc.types == mc.typeVals.rna,        1)' * mc.counts(mc.types == mc.typeVals.rna,        :) / Constants.nAvogadro * 1e15;
            this.protein    = mc.mws(mc.types == mc.typeVals.protein,    1)' * mc.counts(mc.types == mc.typeVals.protein,    :) / Constants.nAvogadro * 1e15;
            
            cIdxs = [this.cIdx.c; this.cIdx.m];
            
            this.cell(:) = 0;
            this.cell(cIdxs) = ...
                + this.metabolite(cIdxs) ...
                + this.rna(cIdxs) ...
                + this.protein(cIdxs);
            
            this.cellDry(:) = 0;
            this.cellDry(cIdxs) = this.cell(cIdxs) - mc.mws(mc.idx.h2o, 1)' * mc.counts(mc.idx.h2o, cIdxs) / Constants.nAvogadro * 1e15;
        end
    end
end