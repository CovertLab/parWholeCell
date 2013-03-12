%Translation
% Translation sub-model. Encodes molecular simulation of bacterial
% translation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
classdef Translation < wholecell.sim.process.Process
    properties
        meta = wholecell.util.struct(...
            'id', 'Translation', ...
            'name', 'Translation' ...
            )
    end
    
    %references to states
    properties
        metabolite
        mrna
        protein
        enzyme
    end
    
    %constants
    properties
        elngRate = 16    %aa/s
        proteinAaCounts  %protein amino acid counts [aa X protein]
        proteinLens      %protein lengths
    end
    
    methods
        %constructor
        function this = Translation()
            this = this@wholecell.sim.process.Process();
        end
        
        %construct object graph, calculate constants
        function this = initialize(this, sim, kb)
            this.initialize@wholecell.sim.process.Process(sim, kb);
            
            %metabolites
            this.metabolite = sim.getState('MoleculeCounts').addPartition(this, {
                'ALA[c]'; 'ARG[c]'; 'ASN[c]'; 'ASP[c]'; 'CYS[c]'; 'GLU[c]'; 'GLN[c]'; 'GLY[c]'; 'HIS[c]'; 'ILE[c]';  'LEU[c]';
                'LYS[c]'; 'MET[c]'; 'PHE[c]'; 'PRO[c]'; 'SER[c]'; 'THR[c]'; 'TRP[c]'; 'TYR[c]'; 'VAL[c]';
                'FMET[c]';
                'GTP[c]'; 'GDP[c]'; 'PI[c]';  'H2O[c]'; 'H[c]'
                }, @this.calcReqMetabolites);
            this.metabolite.idx.aas = this.metabolite.getIndex({
                'ALA[c]'; 'ARG[c]'; 'ASN[c]'; 'ASP[c]'; 'CYS[c]'; 'GLU[c]'; 'GLN[c]'; 'GLY[c]'; 'HIS[c]'; 'ILE[c]';  'LEU[c]';
                'LYS[c]'; 'MET[c]'; 'PHE[c]'; 'PRO[c]'; 'SER[c]'; 'THR[c]'; 'TRP[c]'; 'TYR[c]'; 'VAL[c]';
                });
            this.metabolite.idx.gtp = this.metabolite.getIndex('GTP[c]');
            this.metabolite.idx.gdp = this.metabolite.getIndex('GDP[c]');
            this.metabolite.idx.pi = this.metabolite.getIndex('PI[c]');
            this.metabolite.idx.h2o = this.metabolite.getIndex('H2O[c]');
            this.metabolite.idx.h = this.metabolite.getIndex('H[c]');
            
            %mRNA, protein monomer
            mrnas = kb.rnas(strcmp({kb.rnas.type}, 'mRNA'));
            monomers = kb.proteins([kb.proteins.monomer]);
            this.mrna = sim.getState('MoleculeCounts').addPartition(this, ...
                arrayfun(@(x) [x.id ':mature[c]'], mrnas, 'UniformOutput', false), ...
                @this.calcReqMrna);
            this.protein = sim.getState('MoleculeCounts').addPartition(this, ...
                arrayfun(@(x) [x.monomerId ':nascent[c]'], mrnas, 'UniformOutput', false), ...
                @this.calcReqProtein);
            this.proteinAaCounts = [monomers.aaCount];
            this.proteinLens = sum(this.proteinAaCounts, 1)';
            
            %enyzmes
            this.enzyme = sim.getState('MoleculeCounts').addPartition(this, {
                'RIBOSOME_70S:mature[c]' %70S ribosome
                }, @this.calcReqEnzyme);
            this.enzyme.idx.ribosome70S = this.enzyme.getIndex('RIBOSOME_70S:mature[c]');
        end
        
        %calculate needed metabolites
        function val = calcReqMetabolites(this)
            val = zeros(size(this.metabolite.fullCounts));
            
            elng = this.enzyme.fullCounts(this.enzyme.idx.ribosome70S) * this.elngRate * this.timeStepSec;
            val(this.metabolite.idx.aas) = elng / numel(this.metabolite.idx.aas);
            val([this.metabolite.idx.gtp; this.metabolite.idx.h2o]) = 2 * elng;
        end
        
        %calculate needed mRNA
        function val = calcReqMrna(this)
            val = ones(size(this.enzyme.fullCounts));
        end
        
        %calculate needed protein monomers
        function val = calcReqProtein(this)
            val = zeros(size(this.enzyme.fullCounts));
        end
        
        %calculate needed enzymes
        function val = calcReqEnzyme(this)
            val = ones(size(this.enzyme.fullCounts));
        end
        
        %calculate temporal evolution
        function this = evolveState(this)
            if ~any(this.mrna.counts)
                return;
            end
            
            %total synthesis rate
            proteinSynthProb = this.mrna.counts / sum(this.mrna.counts);
            totRate = 1 / (this.proteinLens' * proteinSynthProb) * min([
                sum(this.metabolite.counts(this.metabolite.idx.aas))
                sum(this.metabolite.counts(this.metabolite.idx.gtp)) / 2
                this.enzyme.counts(this.enzyme.idx.ribosome70S) * this.elngRate * this.timeStepSec
                ]);
            
            %Gillespie algorithm
            t = 0;
            while true
                %chose time step
                t = t + -log(this.randStream.rand()) / totRate;
                if t > this.timeStepSec
                    break;
                end
                
                %check if sufficient metabolic resources to make protein
                newIdx = find(this.randStream.mnrnd(1, proteinSynthProb), 1, 'first');
                if ...
                        any(this.proteinAaCounts(:, newIdx) > this.metabolite.counts(this.metabolite.idx.aas)) || ...
                        2 * this.proteinLens(newIdx) > this.metabolite.counts(this.metabolite.idx.gtp) || ...
                        2 * this.proteinLens(newIdx) > this.metabolite.counts(this.metabolite.idx.h2o)
                    break;
                end
                
                %update metabolites
                this.metabolite.counts(this.metabolite.idx.aas) = this.metabolite.counts(this.metabolite.idx.aas) - this.proteinAaCounts(:, newIdx);
                this.metabolite.counts(this.metabolite.idx.h2o) = this.metabolite.counts(this.metabolite.idx.h2o) + (this.proteinLens(newIdx) - 1);
                
                this.metabolite.counts(this.metabolite.idx.gtp) = this.metabolite.counts(this.metabolite.idx.gtp) - 2 * this.proteinLens(newIdx);
                this.metabolite.counts(this.metabolite.idx.h2o) = this.metabolite.counts(this.metabolite.idx.h2o) - 2 * this.proteinLens(newIdx);
                this.metabolite.counts(this.metabolite.idx.gdp) = this.metabolite.counts(this.metabolite.idx.gdp) + 2 * this.proteinLens(newIdx);
                this.metabolite.counts(this.metabolite.idx.pi)  = this.metabolite.counts(this.metabolite.idx.pi)  + 2 * this.proteinLens(newIdx);
                this.metabolite.counts(this.metabolite.idx.h)   = this.metabolite.counts(this.metabolite.idx.h)   + 2 * this.proteinLens(newIdx);
                
                %increment protein monomer
                this.protein.counts(newIdx) = this.protein.counts(newIdx) + 1;
            end
        end
    end
end