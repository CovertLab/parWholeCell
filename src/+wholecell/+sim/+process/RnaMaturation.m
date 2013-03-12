%RnaMaturation
% RNA maturation sub-model. Encodes molecular simulation of RNA
% maturation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/4/2013
classdef RnaMaturation < wholecell.sim.process.Process
    properties
        meta = wholecell.util.struct(...
            'id', 'RnaMaturation', ...
            'name', 'RNA maturation' ...
            )
    end
    
    %references to states
    properties
        nascentRna
        matureRna
    end
    
    methods
        %constructor
        function this = RnaMaturation()
            this = this@wholecell.sim.process.Process();
        end
        
        %construct object graph
        function this = initialize(this, sim, kb)
            this.initialize@wholecell.sim.process.Process(sim, kb);
            
            this.nascentRna = sim.getState('MoleculeCounts').addPartition(this, ...
                arrayfun(@(x) [x.id ':nascent[c]'], kb.rnas, 'UniformOutput', false), ...
                @this.calcReqNascentRna);
            this.matureRna = sim.getState('MoleculeCounts').addPartition(this, ...
                arrayfun(@(x) [x.id ':mature[c]'], kb.rnas, 'UniformOutput', false), ...
                @this.calcReqMatureRna);
        end
        
        %calculate needed RNA
        function val = calcReqNascentRna(this)
            val = ones(size(this.nascentRna.fullCounts));
        end
        
        %calculate needed RNA
        function val = calcReqMatureRna(this)
            val = zeros(size(this.matureRna.fullCounts));
        end
        
        %calculate temporal evolution
        function this = evolveState(this)
            this.matureRna.counts = ...
                + this.matureRna.counts ...
                + this.nascentRna.counts;
            this.nascentRna.counts(:) = 0;
        end
    end
end