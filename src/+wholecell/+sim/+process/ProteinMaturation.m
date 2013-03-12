%ProteinMaturation
% Protein maturation sub-model. Encodes molecular simulation of protein
% maturation: processing, localization.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/4/2013
classdef ProteinMaturation < wholecell.sim.process.Process
    properties
        meta = wholecell.util.struct(...
            'id', 'ProteinMaturation', ...
            'name', 'Protein maturation' ...
            )
    end
    
    %references to states
    properties
        nascentProteinMonomer
        matureProteinMonomer
    end
    
    methods
        %constructor
        function this = ProteinMaturation()
            this = this@wholecell.sim.process.Process();
        end
        
        %construct object graph
        function this = initialize(this, sim, kb)
            this.initialize@wholecell.sim.process.Process(sim, kb);
            
            monomers = kb.proteins([kb.proteins.monomer]);
            
            this.nascentProteinMonomer = sim.getState('MoleculeCounts').addPartition(this, ...
                arrayfun(@(x) [x.id ':nascent[c]'], monomers, 'UniformOutput', false), ...
                @this.calcReqNascentProteinMonomer);
            this.matureProteinMonomer = sim.getState('MoleculeCounts').addPartition(this, ...
                arrayfun(@(x) [x.id ':mature[' x.compartment ']'], monomers, 'UniformOutput', false), ...
                @this.calcReqMatureProteinMonomer);
        end
        
        %calculate needed proteins
        function val = calcReqNascentProteinMonomer(this)
            val = ones(size(this.nascentProteinMonomer.fullCounts));
        end
        
        %calculate needed proteins
        function val = calcReqMatureProteinMonomer(this)
            val = zeros(size(this.matureProteinMonomer.fullCounts));
        end
        
        %calculate temporal evolution
        function this = evolveState(this)
            this.matureProteinMonomer.counts = ...
                + this.matureProteinMonomer.counts ...
                + this.nascentProteinMonomer.counts;
            this.nascentProteinMonomer.counts(:) = 0;
        end
    end
end