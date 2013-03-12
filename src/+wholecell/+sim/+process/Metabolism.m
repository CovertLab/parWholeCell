%Metabolism
% Metabolism submodel. Encodes molecular simulation of microbial metabolism
% using flux-balance analysis [1,2].
%
% References
% 1. Orth JD, Thiele I, Palsson BO (2010). What is flux balance analysis?
%    Nat Biotechnol. 28(3):245-8.
% 2. Thiele I, Palsson BO (2010). A protocol for generating a
%    high-quality genome-scale metabolic reconstruction. Nat Protoc.
%    5(1): 93-121.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
classdef Metabolism < wholecell.sim.process.Process
    properties
        meta = wholecell.util.struct(...
            'id', 'Metabolism', ...
            'name', 'Metabolism', ...
            'options', {'lpSolver'; 'realMax'} ...
            )
    end
    
    %options
    properties
        lpSolver = 'glpk'
        realMax = 1e6
    end
    
    %references to states
    properties        
        metabolism
        metabolite
        enzyme
        mass
    end
    
    %constants
    properties
        avgCellInitMass = 13.1     %fg
        cellCycleLen = 9 * 3600    %s
        
        objective                  %FBA LP objective (max growth)
        sMat                       %stoichiometry matrix [met x rxn]
        eMat                       %enzyme catalysis matrix [rxn x enz]
        bounds                     %flux bounds data
        dConc_dt                   %dCont/dt (FBA LP RHS)
        metConstTypes              %metabolite constraint types
        rxnVarTypes                %rxn variable types
        
        metIds                     %metabolite ids
        metNames                   %metabolite names
        metIdx                     %metabolite indices
        
        rxnIds                     %reaction ids
        rxnNames                   %reaction names
        rxnIdx                     %reaction indices
    end
    
    methods
        %constructor
        function this = Metabolism()
            this = this@wholecell.sim.process.Process();
        end
        
        %construct object graph
        function this = initialize(this, sim, kb)
            this.initialize@wholecell.sim.process.Process(sim, kb);
            
            %list of metabolites, enzymes
            molIds = cell(0, 1);
            enzIds = cell(0, 1);
            for iReaction = 1:numel(kb.reactions)
                r = kb.reactions(iReaction);
                for iMolecule = 1:numel(kb.reactions(iReaction).stoichiometry)
                    s = r.stoichiometry(iMolecule);
                    molIds{end + 1, 1} = sprintf('%s:%s[%s]', s.molecule, s.form, s.compartment); %#ok<AGROW>
                end
                e = r.enzyme;
                if ~isempty(e)
                    enzIds{end + 1, 1} = sprintf('%s:%s[%s]', e.id, e.form, e.compartment); %#ok<AGROW>
                end
            end
            molIds = unique(molIds);
            enzIds = unique(enzIds);
            
            %partitions
            this.metabolism = sim.getState('Metabolism').addPartition(this);
            this.metabolite = sim.getState('MoleculeCounts').addPartition(this, molIds, @this.calcReqMetabolites);
            this.enzyme = sim.getState('MoleculeCounts').addPartition(this, enzIds, @this.calcReqEnzyme);
            this.mass = sim.getState('Mass').addPartition(this);
            
            this.metabolite.idx.ntps = this.metabolite.getIndex({'ATP[c]'; 'CTP[c]'; 'GTP[c]'; 'UTP[c]'});
            this.metabolite.idx.ndps = this.metabolite.getIndex({'ADP[c]'; 'CDP[c]'; 'GDP[c]'; 'UDP[c]'});
            this.metabolite.idx.nmps = this.metabolite.getIndex({'AMP[c]'; 'CMP[c]'; 'GMP[c]'; 'UMP[c]'});
            this.metabolite.idx.ppi  = this.metabolite.getIndex('PPI[c]');
            this.metabolite.idx.pi   = this.metabolite.getIndex('PI[c]');
            this.metabolite.idx.h2o  = this.metabolite.getIndex('H2O[c]');
            this.metabolite.idx.h    = this.metabolite.getIndex('H[c]');
            
            %indices
            nMet = numel(molIds) + 1;
            this.metIds = [this.metabolite.ids; 'biomass'];
            this.metNames = [this.metabolite.names; 'biomass'];
            this.metIdx = struct();
            this.metIdx.real = (1:numel(molIds))';
            this.metIdx.biomass = this.metIdx.real(end) + 1;
            
            nRxn = numel(kb.reactions) + numel(molIds) + 2;
            mc = sim.getState('MoleculeCounts');
            [~, cIdxs] = ind2sub([numel(mc.ids) numel(mc.compartments)], this.metabolite.mapping);
            this.rxnIds = [
                this.metabolism.reactionIds
                cellfun(@(x, y) ['ex_' x '_' y], this.metabolite.ids, {mc.compartments(cIdxs).id}', 'UniformOutput', false)
                'growth' 
                'ex_biomass'
                ];
            this.rxnNames = [
                this.metabolism.reactionNames
                cellfun(@(x, y) [x ' exchange (' y ')'], this.metabolite.names, {mc.compartments(cIdxs).name}', 'UniformOutput', false)
                'growth'
                'biomass exchange'
                ];
            this.rxnIdx = struct();
            this.rxnIdx.real = (1:numel(kb.reactions))';
            this.rxnIdx.exchange = this.rxnIdx.real(end) + (1:numel(molIds))';
            this.rxnIdx.growth = this.rxnIdx.exchange(end) + 1;
            this.rxnIdx.biomassExchange = this.rxnIdx.growth(end) + 1;
            
            mc = sim.getState('MoleculeCounts');
            [~, ~, cIdxs] = mc.getIndex(molIds);
            iExtracellular = find(strcmp({mc.compartments.id}, 'e'), 1, 'first');
            this.rxnIdx.internalExchange = this.rxnIdx.exchange(cIdxs ~= iExtracellular);
            this.rxnIdx.externalExchange = this.rxnIdx.exchange(cIdxs == iExtracellular);
            
            nEnz = numel(enzIds);
            
            %objective
            this.objective = zeros(nRxn, 1);
            this.objective(this.rxnIdx.growth) = 1;
            
            %stoichiometry matrix, enzymes, kinetics
            this.sMat = zeros(nMet, nRxn);
            this.sMat(this.metIdx.real, this.rxnIdx.exchange) = eye(numel(molIds));
            this.sMat(this.metIdx.biomass, this.rxnIdx.growth) = 1;
            this.sMat(this.metIdx.biomass, this.rxnIdx.biomassExchange) = -1;
            
            this.eMat = zeros(nRxn, nEnz);
            this.bounds = struct(...
                'kinetic',       struct('lo', -inf(nRxn, 1), 'up', inf(nRxn, 1)), ...
                'thermodynamic', struct('lo', -inf(nRxn, 1), 'up', inf(nRxn, 1)), ...
                'exchange',      struct('lo', -inf(nRxn, 1), 'up', inf(nRxn, 1)) ...
                );
            this.bounds.thermodynamic.lo(this.rxnIdx.growth) = 0;
            this.bounds.thermodynamic.lo(this.rxnIdx.biomassExchange) = 0;
            
            tmpSMat = cell(0, 3);
            tmpEMat = cell(0, 2);
            for iReaction = 1:numel(kb.reactions)
                r = kb.reactions(iReaction);
                rIdx = this.rxnIdx.real(iReaction);
                
                %stoichiometry
                for iMolecule = 1:numel(kb.reactions(iReaction).stoichiometry)
                    s = r.stoichiometry(iMolecule);                    
                    tmpSMat = [tmpSMat; {sprintf('%s:%s[%s]', s.molecule, s.form, s.compartment) rIdx s.coeff}]; %#ok<AGROW>
                end
                
                %enzyme
                if ~isempty(r.enzyme)
                    %catalysis
                    e = r.enzyme;
                    tmpEMat = [tmpEMat; {sprintf('%s:%s[%s]', e.id, e.form, e.compartment) rIdx}]; %#ok<AGROW>                   
                    
                    %kinetics
                    if ~isnan(e.kCatRev)
                        this.bounds.kinetic.lo(rIdx) = -e.kCatRev;
                    end
                    if ~isnan(e.kCatFor)
                        this.bounds.kinetic.up(rIdx) = e.kCatFor;
                    end
                end
                
                %thermodynamics
                if r.dir == 1
                    this.bounds.thermodynamic.lo(rIdx) = 0;
                elseif r.dir == -1
                    this.bounds.thermodynamic.up(rIdx) = 0;
                end
            end
            mIdx = this.metIdx.real(this.metabolite.getIndex(tmpSMat(:, 1)));
            this.sMat(sub2ind(size(this.sMat), mIdx, cell2mat(tmpSMat(:, 2)))) = cell2mat(tmpSMat(:, 3));
            
            eIdx = this.enzyme.getIndex(tmpEMat(:, 1));
            this.eMat(sub2ind(size(this.eMat), cell2mat(tmpEMat(:, 2)), eIdx)) = 1;
            
            %exchange
            metIds = arrayfun(@(x) [x.id ':mature[e]'], kb.metabolites, 'UniformOutput', false);
            metExs = [kb.metabolites.maxExchangeRate]';
            
            tfs = ismember(metIds, molIds);
            metIds = metIds(tfs);
            metExs = metExs(tfs);
            
            this.bounds.exchange.lo(:) = 0;
            this.bounds.exchange.up(:) = 0;
            
            metIdxs = this.metabolite.getIndex(metIds);
            this.bounds.exchange.lo(this.rxnIdx.exchange(metIdxs)) = -metExs;
            this.bounds.exchange.up(this.rxnIdx.exchange(metIdxs)) =  metExs;
            
            %objective
            objMets = kb.metabolites([kb.metabolites.metabolismFlux] ~= 0);
            
            metComps = cell(numel(objMets), 1);
            metComps( [objMets.hydrophobic]) = {'m'};
            metComps(~[objMets.hydrophobic]) = {'c'};
            
            realMetIds = cellfun(@(x, y) [x '[' y ']'], {objMets.id}', metComps, 'UniformOutput', false);
            realMetIdxs = this.metabolite.getIndex(realMetIds);
            
            this.sMat(this.metIdx.real(realMetIdxs), this.rxnIdx.growth) = -[objMets(realMetIdxs ~= 0).metabolismFlux]';
            
            %dConc/dt
            this.dConc_dt = zeros(nMet, 1);
            
            %constraint, variable types
            this.metConstTypes = repmat('S', nMet, 1);
            this.rxnVarTypes = repmat('C', nRxn, 1);
            
            %more indices
            this.rxnIdx.catalyzed = find(any(this.eMat, 2));
        end
        
        %calculate needed metabolites
        function val = calcReqMetabolites(this)
            val = ones(size(this.metabolite.fullCounts));
            val(this.metabolite.idx.ntps) = 0;
            val(this.metabolite.idx.h2o) = 0;
        end
        
        %calculate needed proteins
        function val = calcReqEnzyme(this)
            val = ones(size(this.enzyme.fullCounts));
        end
        
        %calculate temporal evolution
        function this = evolveState(this)
            %calculate flux bounds
            bounds = this.calcFluxBounds(this.metabolite.counts, this.enzyme.counts); %#ok<*PROP>
            
            %calculate growth rate
            [this.metabolism.growth, this.metabolism.fluxes, exRates] = this.calcGrowthRate(bounds);
            
            %update metabolite copy numbers
            this.metabolite.counts = round(...
                + this.metabolite.counts ...
                - this.sMat(this.metIdx.real, this.rxnIdx.growth) * (this.metabolism.growth / (3600 * this.avgCellInitMass)) * this.timeStepSec ...
                + exRates * this.timeStepSec ...
                );
            
%             %NTP-->NDP
%             tmp = min(0, this.metabolite.counts(this.metabolite.idx.ndps));
%             this.metabolite.counts(this.metabolite.idx.ntps) = ...
%                 + this.metabolite.counts(this.metabolite.idx.ntps) ...
%                 + tmp;
%             this.metabolite.counts(this.metabolite.idx.h2o) = ...
%                 + this.metabolite.counts(this.metabolite.idx.h2o) ...
%                 + sum(tmp);
%             this.metabolite.counts(this.metabolite.idx.ndps) = ...
%                 + this.metabolite.counts(this.metabolite.idx.ndps) ...
%                 - tmp;
%             this.metabolite.counts(this.metabolite.idx.pi) = ...
%                 + this.metabolite.counts(this.metabolite.idx.pi) ...
%                 - sum(tmp);
%             this.metabolite.counts(this.metabolite.idx.h) = ...
%                 + this.metabolite.counts(this.metabolite.idx.h) ...
%                 - sum(tmp);
%             
%             %NTP-->NMP
%             tmp = min(0, this.metabolite.counts(this.metabolite.idx.nmps));
%             this.metabolite.counts(this.metabolite.idx.ntps) = ...
%                 + this.metabolite.counts(this.metabolite.idx.ntps) ...
%                 + tmp;
%             this.metabolite.counts(this.metabolite.idx.h2o) = ...
%                 + this.metabolite.counts(this.metabolite.idx.h2o) ...
%                 + sum(tmp);
%             this.metabolite.counts(this.metabolite.idx.nmps) = ...
%                 + this.metabolite.counts(this.metabolite.idx.nmps) ...
%                 - tmp;
%             this.metabolite.counts(this.metabolite.idx.ppi) = ...
%                 + this.metabolite.counts(this.metabolite.idx.ppi) ...
%                 - sum(tmp);
%             this.metabolite.counts(this.metabolite.idx.h) = ...
%                 + this.metabolite.counts(this.metabolite.idx.h) ...
%                 - sum(tmp);
        end
        
        function [growth, fluxes, exchangeRates] = calcGrowthRate(this, bounds)
            import wholecell.util.linearProgramming;
            
            %cap bounds
            bounds.lo = min(0, max(bounds.lo * this.timeStepSec, -this.realMax));
            bounds.up = max(0, min(bounds.up * this.timeStepSec,  this.realMax));
            
            %flux-balance analysis 
            %TODO: make flexible similar to full whole-cell model
            [x, ~, ~, errorFlag, errorMsg] = linearProgramming(...
                'maximize', this.objective, ...
                this.sMat, this.dConc_dt, ...
                bounds.lo, bounds.up, ...
                this.metConstTypes, ...
                this.rxnVarTypes, ...
                struct('solver', this.lpSolver, 'solverOptions', struct('glpk', struct())));
            if errorFlag
                throw(MException('Metabolism:linearProgramming:error', 'Linear programming error: %s', errorMsg));
            end
            
            %extract growth, real fluxes
            x = x  / this.timeStepSec;
            growth = x(this.rxnIdx.growth) * 3600 * this.avgCellInitMass;
            fluxes = x(this.rxnIdx.real);
            exchangeRates = x(this.rxnIdx.exchange);
        end
        
        function val = calcFluxBounds(this, metCnts, enzCnts, ...
                applyThermoBounds, applyKineticBounds, applyExchangeBounds, ...
                applyMetAvailBounds)
            if nargin < 4, applyThermoBounds = true; end;
            if nargin < 5, applyKineticBounds = true; end;
            if nargin < 6, applyExchangeBounds = true; end;
            if nargin < 7, applyMetAvailBounds = true; end;

            %initialize
            nRxn = numel(this.rxnIds);
            lo = -inf(nRxn, 1);
            up =  inf(nRxn, 1);
            
            %thermodynamics
            if applyThermoBounds
                lo = max(lo, this.bounds.thermodynamic.lo);
                up = min(up, this.bounds.thermodynamic.up);
            end
            
            %kinetics
            if applyKineticBounds
                rxnEnzCnts = this.eMat(this.rxnIdx.catalyzed, :) * enzCnts;
                
                kLo = rxnEnzCnts .* this.bounds.kinetic.lo(this.rxnIdx.catalyzed);
                kUp = rxnEnzCnts .* this.bounds.kinetic.up(this.rxnIdx.catalyzed);
                kLo(rxnEnzCnts == 0) = 0;
                kUp(rxnEnzCnts == 0) = 0;
                
                lo(this.rxnIdx.catalyzed) = max(lo(this.rxnIdx.catalyzed), kLo);
                up(this.rxnIdx.catalyzed) = min(up(this.rxnIdx.catalyzed), kUp);
            end
            
            %exchange
            if applyExchangeBounds
                lo(this.rxnIdx.internalExchange) = 0;
                up(this.rxnIdx.internalExchange) = 0;
                
                lo(this.rxnIdx.externalExchange) = max(lo(this.rxnIdx.externalExchange), ...
                    this.bounds.exchange.lo(this.rxnIdx.externalExchange) * 6.022e23 * 1e-3 * 3600 * sum(this.mass.cellDry) * 1e-15);
                up(this.rxnIdx.externalExchange) = min(up(this.rxnIdx.externalExchange), ...
                    this.bounds.exchange.up(this.rxnIdx.externalExchange) * 6.022e23 * 1e-3 * 3600 * sum(this.mass.cellDry) * 1e-15);
            end
            
            %metabolite availability
            if applyMetAvailBounds
                idxs = find(ismember(this.rxnIdx.exchange, this.rxnIdx.externalExchange));
                idxs = idxs(setdiff(1:end, [2 54 95 126 137 195 197 200 208 211 227])); %Todo
                
                lo(this.rxnIdx.exchange(idxs)) = max(lo(this.rxnIdx.exchange(idxs)), -metCnts(idxs) / this.timeStepSec);
            end
            
            %protein %TODO
            
            %return value
            val = struct('lo', lo, 'up', up);
        end
    end
end