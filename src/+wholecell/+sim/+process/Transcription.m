%Transcription
% Transcription sub-model. Encodes molecular simulation of macromolecular
% bacterial transcription.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/4/2013
classdef Transcription < wholecell.sim.process.Process
    properties
        meta = wholecell.util.struct(...
            'id', 'Transcription', ...
            'name', 'Transcription' ...
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
        cellCycleLength = 9 * 3600 %s
        elngRate = 50              %nt/s
        rnaLens                    %RNA lengths
        rnaNtCounts                %RNA nucleotide counts [nt x RNA]
        rnaSynthProb               %relative RNA synthesis rates
    end
    
    
    methods
        %constructor
        function this = Transcription()
            this = this@wholecell.sim.process.Process();
        end
        
        %construct object graph
        function this = initialize(this, sim, kb)
            this.initialize@wholecell.sim.process.Process(sim, kb);
            
            %metabolites
            this.metabolite = sim.getState('MoleculeCounts').addPartition(this, {
                'ATP[c]'; 'CTP[c]'; 'GTP[c]'; 'UTP[c]'
                'PPI[c]'; 'H2O[c]'; 'H[c]'
                }, @this.calcReqMetabolites);
            this.metabolite.idx.ntps = this.metabolite.getIndex({'ATP[c]'; 'CTP[c]'; 'GTP[c]'; 'UTP[c]'});
            this.metabolite.idx.ppi = this.metabolite.getIndex({'PPI[c]'});
            this.metabolite.idx.h2o = this.metabolite.getIndex({'H2O[c]'});
            this.metabolite.idx.h = this.metabolite.getIndex({'H[c]'});
            
            %rna
            this.rna = sim.getState('MoleculeCounts').addPartition(this, ...
                arrayfun(@(x) [x.id ':nascent[c]'], kb.rnas, 'UniformOutput', false), ...
                @this.calcReqRna);
            this.rnaNtCounts = [kb.rnas.ntCount];
            this.rnaLens = sum(this.rnaNtCounts, 1)';            
            this.rnaSynthProb = [kb.rnas.exp]' .* (log(2) / this.cellCycleLength + 1 ./ [kb.rnas.halfLife]');
            this.rnaSynthProb = this.rnaSynthProb / sum(this.rnaSynthProb);
            
            %enyzmes
            this.enzyme = sim.getState('MoleculeCounts').addPartition(this, {
                'RNA_POLYMERASE:mature[c]'; %DNA-directed RNA polymerase
                }, @this.calcReqEnzyme);
            this.enzyme.idx.rnaPol = this.enzyme.getIndex({
                'RNA_POLYMERASE:mature[c]'; %DNA-directed RNA polymerase
                });
        end
        
        %calculate needed metabolites
        function val = calcReqMetabolites(this)
            val = zeros(size(this.metabolite.fullCounts));
            
            val(this.metabolite.idx.ntps) = min([
                this.enzyme.fullCounts(this.enzyme.idx.rnaPol) * this.elngRate * this.timeStepSec
                4 * min(this.metabolite.fullCounts(this.metabolite.idx.ntps))
                ]) / 4;
            
            val(this.metabolite.idx.h2o) = 1;
        end
        
        %calculate needed RNA
        function val = calcReqRna(this)
            val = zeros(size(this.rna.fullCounts));
        end
        
        %calculate needed enzymes
        function val = calcReqEnzyme(this)
            val = ones(size(this.enzyme.fullCounts));
        end
        
        %calculate temporal evolution
        function this = evolveState(this)
            %total synthesis rate
            totRate = min(...
                sum(this.metabolite.counts(this.metabolite.idx.ntps)), ...
                this.enzyme.counts(this.enzyme.idx.rnaPol) * this.elngRate * this.timeStepSec ...
                ) / (this.rnaLens' * this.rnaSynthProb);
            
            %Gillespie algorithm
            t = 0;
            while true
                %chose time step
                t = t + -log(this.randStream.rand()) / totRate;
                if t > this.timeStepSec
                    break;
                end
                
                %check if sufficient metabolic resources to make RNA
                newIdx = find(this.randStream.mnrnd(1, this.rnaSynthProb), 1, 'first');
                if ...
                        any(this.rnaNtCounts(:, newIdx) > this.metabolite.counts(this.metabolite.idx.ntps)) || ...
                        this.metabolite.counts(this.metabolite.idx.h2o) < 1
                    break;
                end
                
                %update metabolites
                this.metabolite.counts(this.metabolite.idx.ntps) = this.metabolite.counts(this.metabolite.idx.ntps) - this.rnaNtCounts(:, newIdx);
                this.metabolite.counts(this.metabolite.idx.h2o)  = this.metabolite.counts(this.metabolite.idx.h2o)  - 1;
                this.metabolite.counts(this.metabolite.idx.ppi)  = this.metabolite.counts(this.metabolite.idx.ppi)  + this.rnaLens(newIdx);
                this.metabolite.counts(this.metabolite.idx.h)    = this.metabolite.counts(this.metabolite.idx.h)    + 1;
                
                %increment RNA
                this.rna.counts(newIdx) = this.rna.counts(newIdx) + 1;
            end
        end
    end
end