%Complexation
% Macromolecular complexation sub-model. Encodes molecular simulation of
% Macromolecular complexation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/4/2013
classdef Complexation < wholecell.sim.process.Process
    properties
        meta = wholecell.util.struct(...
            'id', 'Complexation', ...
            'name', 'Macromolecular complexation' ...
            )
    end
    
    %references to states
    properties
        subunit
        complex
    end
    
    %constants
    properties
        sMat %complex subunit composition matrix [subunits x complex]
    end
    
    methods
        %constructor
        function this = Complexation()
            this = this@wholecell.sim.process.Process();
        end
        
        %construct object graph
        function this = initialize(this, sim, kb)
            this.initialize@wholecell.sim.process.Process(sim, kb);
            
            %complex
            complexes = kb.proteins(~[kb.proteins.monomer] & strcmp({kb.proteins.formationProcess}, this.meta.id));
            this.complex = sim.getState('MoleculeCounts').addPartition(this, ...
                arrayfun(@(x) [x.id ':mature[' x.compartment ']'], complexes, 'UniformOutput', false), ...
                @this.calcReqComplex);
            
            %subunits
            subunits = [];
            for iComplex = 1:numel(complexes)
                c = complexes(iComplex);
                subunits = [subunits; c.composition([c.composition.coeff] < 0)]; %#ok<AGROW>
            end
            subIdComps = unique(arrayfun(@(x) [x.molecule ':mature[' x.compartment ']'], subunits, 'UniformOutput', false));
            
            this.subunit = sim.getState('MoleculeCounts').addPartition(this, ...
                subIdComps, @this.calcReqSubunit);
            
            tmpSMat = cell(0, 3);
            for iComplex = 1:numel(complexes)
                c = complexes(iComplex);
                for iSubunit = 1:numel(c.composition)
                    s = c.composition(iSubunit);
                    if s.coeff > 0
                        continue;
                    end
                    
                    tmpSMat = [tmpSMat; {[s.molecule ':mature[' s.compartment ']'] iComplex -s.coeff}]; %#ok<AGROW>
                end
            end
            
            subIdx = this.subunit.getIndex(tmpSMat(:, 1));
            this.sMat = zeros(numel(subIdComps), numel(complexes));
            this.sMat(sub2ind(size(this.sMat), subIdx, cell2mat(tmpSMat(:, 2)))) = cell2mat(tmpSMat(:, 3));
        end
        
        %calculate needed proteins
        function val = calcReqSubunit(this)
            val = ones(size(this.subunit.fullCounts));
        end
        
        %calculate needed proteins
        function val = calcReqComplex(this)
            val = zeros(size(this.complex.fullCounts));
        end
        
        %calculate temporal evolution
        function this = evolveState(this)
            [this.subunit.counts, this.complex.counts] = ...
                this.calcNewComplexes(this.subunit.counts, this.complex.counts, 1);
        end
        
        %Gillespie algorithm
        function [subunits, complexes] = calcNewComplexes(this, subunits, complexes, leap)
            %TODO: implement tau leaping correctly
            while true
                %calculate rates
                rates = floor(min(...
                    subunits(:, ones(size(this.complex.ids))) ./ this.sMat ...
                    , [], 1))';
                totRate = sum(rates);
                rates = rates / totRate;
                if totRate <= 0
                    break;
                end
                
                %check if sufficient metabolic resources to make protein
                newCnts = this.randStream.mnrnd(leap, rates);
                newCnts = newCnts(:);
                if any(this.sMat * newCnts > subunits)
                    break;
                end
                
                %update subunits
                subunits = subunits - this.sMat * newCnts;
                
                %increment complex
                complexes = complexes + newCnts;
            end
        end
    end
end