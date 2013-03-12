%RnaDegradation
% RNA degradation sub-model. Encodes molecular simulation of RNA
% degradation as a Poisson process.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/4/2013
classdef RnaDegradation < wholecell.sim.process.Process
    properties
        meta = wholecell.util.struct(...
            'id', 'RnaDegradation', ...
            'name', 'RNA degradation' ...
            )
    end
    
    %references to states
    properties
        metabolite
        rna
        enzyme
    end
    
    %constants
    properties
        rnaLens     %RNA lengths
        rnaDegRates %RNA degradation rates (1/s)
        rnaDegSMat  %RNA degradation reaction stoichiometry matrix [metabolite x rna]
    end
    
    methods
        %constructor
        function this = RnaDegradation()
            this = this@wholecell.sim.process.Process();
        end
        
        %construct object graph
        function this = initialize(this, sim, kb)
            this.initialize@wholecell.sim.process.Process(sim, kb);
            
            %metabolites
            this.metabolite = sim.getState('MoleculeCounts').addPartition(this, {
                'AMP[c]'; 'CMP[c]'; 'GMP[c]'; 'UMP[c]'; 'H2O[c]'; 'H[c]'
                }, @this.calcReqMetabolites);
            this.metabolite.idx.nmps = this.metabolite.getIndex({'AMP[c]'; 'CMP[c]'; 'GMP[c]'; 'UMP[c]'});
            this.metabolite.idx.h2o = this.metabolite.getIndex('H2O[c]');
            this.metabolite.idx.h = this.metabolite.getIndex('H[c]');
            
            %rna
            this.rna = sim.getState('MoleculeCounts').addPartition(this, [
                arrayfun(@(x) [x.id ':nascent[c]'], kb.rnas, 'UniformOutput', false)
                arrayfun(@(x) [x.id ':mature[c]'], kb.rnas, 'UniformOutput', false)
                ], @this.calcReqRna, true);
            
            this.rnaDegRates = log(2) ./ [
                [kb.rnas.halfLife]'
                [kb.rnas.halfLife]'
                ];
            
            this.rnaLens = sum([[kb.rnas.ntCount] [kb.rnas.ntCount]], 1)';
            
            this.rnaDegSMat = zeros(numel(this.metabolite.ids), numel(this.rna.ids));
            this.rnaDegSMat(this.metabolite.idx.nmps, :) = [[kb.rnas.ntCount] [kb.rnas.ntCount]];
            this.rnaDegSMat(this.metabolite.idx.h2o, :)  = -(sum(this.rnaDegSMat(this.metabolite.idx.nmps, :), 1) - 1);
            this.rnaDegSMat(this.metabolite.idx.h, :)    =  (sum(this.rnaDegSMat(this.metabolite.idx.nmps, :), 1) - 1);
            
            %proteins
            this.enzyme = sim.getState('MoleculeCounts').addPartition(this, {
                'MG_104_MONOMER:mature[c]'; %ribonuclease R
                }, @this.calcReqEnzyme);
            this.enzyme.idx.rnaseR = this.enzyme.getIndex({
                'MG_104_MONOMER:mature[c]'; %ribonuclease R
                });
        end
        
        %calculate needed metabolites
        function val = calcReqMetabolites(this)
            val = zeros(size(this.metabolite.fullCounts));
            val(this.metabolite.idx.h2o) = this.rnaLens' * (this.rnaDegRates .* this.rna.fullCounts) * this.timeStepSec;
        end
        
        %calculate needed metabolites
        function val = calcReqRna(this)
            val = this.randStream.poissrnd(this.rnaDegRates .* this.rna.fullCounts * this.timeStepSec);
        end
        
        %calculate needed proteins
        function val = calcReqEnzyme(this)
            val = ones(size(this.enzyme.fullCounts));
        end
        
        %calculate temporal evolution
        function this = evolveState(this)
            %check if RNAse R expressed
            if this.enzyme.counts(this.enzyme.idx.rnaseR) == 0
                return;
            end
            
            %degrade RNA
            this.metabolite.counts = ...
                + this.metabolite.counts  ...
                + this.rnaDegSMat * this.rna.counts;
            this.rna.counts(:) = 0;
        end
    end
end