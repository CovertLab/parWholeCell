%MoleculeCounts
% State which represents the copy numbers of a class of molecules as an
% array.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
classdef MoleculeCounts < wholecell.sim.state.State
    %metadata
    properties
        meta = wholecell.util.struct(...
            'id', 'MoleculeCounts', ...
            'name', 'Molecule counts', ...
            'dynamics', {'counts'}, ...
                'units', struct(...
                'counts', 'molecules' ...
                ) ...
            )
    end
    
    %references to processes
    properties
        complexation
    end
    
    %constants
    properties
        ids              %molecule ids        
        names            %molecule names
        forms            %molecule forms (eg. nascent, mature, etc.)
        types            %molecule type (e.g. metabolite, RNA, protein)
        mws              %molecular weights (Da)
        localizations    %prefered intracellular localizations
        idx = struct()   %indices over molecules
        
        chamberVolume = 4.3667e-12 %L
        
        metMediaConc     %metabolite media concentration (mM)
        metBiomassConc   %metabolite biomass concentration (molecules/cell)
        
        rnaLens          %RNA lengths
        rnaExp           %mature RNA expression
        
        monLens          %protein monomer lengths
        monExp           %mature protein monomer expression
        
        compartments = [ %compartment ids, names
            struct('id', 'c', 'name', 'Cytosol')
            struct('id', 'e', 'name', 'Extracellular space')
            struct('id', 'm', 'name', 'Membrane')
            ];
        cIdx = struct(...
            'c', 1, ...
            'e', 2, ...
            'm', 3 ...
            )
    end
    
    %dynamical properties
    properties
        counts %molecule counts (molecules x compartments)
    end
    
    %form values
    properties
        formVals = struct(...
            'nascent', 2, ...
            'mature', 1 ...
            );
        typeVals = struct(...
            'metabolite', 1, ...
            'rna', 2, ...
            'protein', 3 ...
            );
    end
    
    %partitioning
    properties
        %used by parent
        partitionedCounts   %molecules partitioned at each time step (molecules x compartments x partitions)
        unpartitionedCounts %molecules not partitioned at each time step (molecules x compartments)
        
        %used by children
        fullCounts %full count in parent
        mapping    %index mapping between parent, partition
        reqFunc    %request function handle
        isReqAbs   %requesting absolute copy number (true/false)
    end
    
    %initialization
    methods
        %constructor
        function this = MoleculeCounts(varargin)
            this = this@wholecell.sim.state.State(varargin{:});
        end
        
        %calculate constants
        function this = initialize(this, sim, kb)
            this.initialize@wholecell.sim.state.State(sim, kb);
            
            this.complexation = sim.getProcess('Complexation');
            
            %molecule identities
            this.ids = [
                {kb.metabolites.id}'
                {kb.rnas.id}'
                {kb.rnas.id}'
                {kb.proteins.id}'
                {kb.proteins.id}'
                ];
            this.forms = [
                repmat(this.formVals.mature, numel(kb.metabolites), 1);
                repmat(this.formVals.nascent, numel(kb.rnas), 1)
                repmat(this.formVals.mature, numel(kb.rnas), 1)
                repmat(this.formVals.nascent, numel(kb.proteins), 1)
                repmat(this.formVals.mature, numel(kb.proteins), 1)
                ];
            this.types = [
                repmat(this.typeVals.metabolite, numel(kb.metabolites), 1);
                repmat(this.typeVals.rna, numel(kb.rnas), 1)
                repmat(this.typeVals.rna, numel(kb.rnas), 1)
                repmat(this.typeVals.protein, numel(kb.proteins), 1)
                repmat(this.typeVals.protein, numel(kb.proteins), 1)
                ];
            this.names = [
                {kb.metabolites.name}'
                {kb.rnas.name}'
                {kb.rnas.name}'
                {kb.proteins.name}'
                {kb.proteins.name}'
                ];
            this.mws = [
                [kb.metabolites.mw]';
                [kb.rnas.mw]'
                [kb.rnas.mw]'
                [kb.proteins.mw]'
                [kb.proteins.mw]'
                ];
            
            [~, this.idx.nmps] = this.getIndex({'AMP[c]'; 'CMP[c]'; 'GMP[c]'; 'UMP[c]'});
            [~, this.idx.aas] = this.getIndex({
                'ALA[c]'; 'ARG[c]'; 'ASN[c]'; 'ASP[c]'; 'CYS[c]'; 'GLU[c]'; 'GLN[c]'; 'GLY[c]'; 'HIS[c]'; 'ILE[c]';  'LEU[c]';
                'LYS[c]'; 'MET[c]'; 'PHE[c]'; 'PRO[c]'; 'SER[c]'; 'THR[c]'; 'TRP[c]'; 'TYR[c]'; 'VAL[c]'
                });
            [~, this.idx.h2o] = this.getIndex('H2O[c]');
            
            %localizations
            metLocs = zeros(numel(kb.metabolites), 1);
            metLocs( [kb.metabolites.hydrophobic]) = this.cIdx.m;
            metLocs(~[kb.metabolites.hydrophobic]) = this.cIdx.c;
            [~, protLocs] = ismember({kb.proteins.compartment}', {this.compartments.id}');
            this.localizations = [
                metLocs
                repmat(this.cIdx.c, numel(kb.rnas), 1)
                repmat(this.cIdx.c, numel(kb.rnas), 1)
                repmat(this.cIdx.c, numel(kb.proteins), 1)
                protLocs
                ];
            
            %composition
            this.rnaLens = sum([kb.rnas.ntCount], 1)';
            this.rnaExp = [kb.rnas.exp]';
            this.rnaExp = this.rnaExp / sum(this.rnaExp);
            this.idx.nascentRna = find(this.types == this.typeVals.rna & this.forms == this.formVals.nascent);
            this.idx.matureRna = find(this.types == this.typeVals.rna & this.forms == this.formVals.mature);
            this.idx.nascentMrna = find(this.types == this.typeVals.rna & this.forms == this.formVals.nascent & ismember(this.ids, {kb.rnas(strcmp({kb.rnas.type}, 'mRNA')).id}));
            this.idx.matureMrna = find(this.types == this.typeVals.rna & this.forms == this.formVals.mature & ismember(this.ids, {kb.rnas(strcmp({kb.rnas.type}, 'mRNA')).id}));
            
            mons = kb.proteins([kb.proteins.monomer]);
            this.monLens = sum([mons.aaCount], 1)';
            this.monExp = [kb.rnas(strcmp({kb.rnas.type}, 'mRNA')).exp]';
            this.monExp = this.monExp / sum(this.monExp);
            [~, this.idx.matureMonomers] = this.getIndex(arrayfun(@(x) [x.id ':mature[' x.compartment ']'], mons, 'UniformOutput', false));
            
            cpxs = kb.proteins(~[kb.proteins.monomer]);
            [~, this.idx.matureComplexes] = this.getIndex(arrayfun(@(x) [x.id ':mature[' x.compartment ']'], cpxs, 'UniformOutput', false));
            
            %media biomass
            this.metMediaConc = [kb.metabolites.mediaConc]';   %mM
            this.metBiomassConc = [kb.metabolites.biomassConc]'; %molecules/cell
        end
        
        %allocate memory
        function this = allocate(this)
            this.allocate@wholecell.sim.state.State();
            
            this.counts = zeros(numel(this.ids), numel(this.compartments));
            if isempty(this.parentState)
                this.partitionedCounts = zeros(numel(this.ids), numel(this.compartments), numel(this.partitions));
                this.unpartitionedCounts = zeros(numel(this.ids), numel(this.compartments));
            else
                this.fullCounts = zeros(numel(this.ids), numel(this.compartments));
            end
        end
        
        %calculate initial conditions
        function this = calcInitialConditions(this)
            import wholecell.util.Constants;
            
            this.counts(:) = 0;
            
            %media metabolites            
            this.counts(this.types == this.typeVals.metabolite, this.cIdx.e) = ...
                round(this.metMediaConc * this.chamberVolume * Constants.nAvogadro * 1e-3);
            
            %biomass metabolites
            metIdx = find(this.types == this.typeVals.metabolite);
            metCompIdxs = sub2ind(size(this.counts), metIdx, this.localizations(metIdx));
            this.counts(metCompIdxs) = round(this.metBiomassConc);
            
            %RNA
            rnaCnts = this.randStream.mnrnd(round(sum(this.counts(this.idx.nmps, this.cIdx.c)) / (this.rnaExp' * this.rnaLens)), this.rnaExp);
            this.counts(this.idx.nmps, this.cIdx.c) = 0;
            rnaCompIdxs = sub2ind(size(this.counts), this.idx.matureRna, this.localizations(this.idx.matureRna));
            this.counts(rnaCompIdxs) = rnaCnts;
            
            %protein monomers
            monCnts = this.randStream.mnrnd(round(sum(this.counts(this.idx.aas, this.cIdx.c)) / (this.monExp' * this.monLens)), this.monExp);
            this.counts(this.idx.aas, this.cIdx.c) = 0;
            protCompIdxs = sub2ind(size(this.counts), this.idx.matureMonomers, this.localizations(this.idx.matureMonomers));
            this.counts(protCompIdxs) = monCnts;
            
            %macromolecular complexation
            c = this.complexation;
            
            c.subunit.counts = this.counts(c.subunit.mapping);
            c.complex.counts = this.counts(c.complex.mapping);
            
            [c.subunit.counts, c.complex.counts] = c.calcNewComplexes(c.subunit.counts, c.complex.counts, 100);
            [c.subunit.counts, c.complex.counts] = c.calcNewComplexes(c.subunit.counts, c.complex.counts, 10);
            [c.subunit.counts, c.complex.counts] = c.calcNewComplexes(c.subunit.counts, c.complex.counts, 1);
            
            this.counts(c.subunit.mapping) = c.subunit.counts;
            this.counts(c.complex.mapping) = c.complex.counts;
        end
    end
    
    %partitioning into substates
    methods
        function partition = addPartition(this, process, reqMols, reqFunc, isReqAbs)
            if nargin < 5
                isReqAbs = false;
            end
            
            %super class
            partition = this.addPartition@wholecell.sim.state.State(process);
            
            %clear inherited properties only valid on parent
            partition.compartments = struct('id', 'merged__', 'name', 'merged');
            partition.partitionedCounts = [];
            partition.unpartitionedCounts = [];
            partition.idx = struct();
            
            %set child properties for mapping to parent
            [iMolFormComp, iMolForm] = this.getIndex(reqMols);
            if numel(unique(iMolFormComp)) < numel(iMolFormComp)
                throw(MException('MoleculeCounts:duplicateIds', 'Partition request cannot contain duplicate ids'));
            end
            
            partition.ids = this.ids(iMolForm);
            partition.names = this.names(iMolForm);
            partition.forms = this.forms(iMolForm);
            partition.types = this.types(iMolForm);            
            partition.mws = this.mws(iMolForm);
            
            partition.mapping = iMolFormComp;
            partition.reqFunc = reqFunc;
            partition.isReqAbs = isReqAbs;
        end
        
        %prepare to partition state among processes
        function prepartition(this)
            for iPartition = 1:numel(this.partitions)
                p = this.partitions{iPartition};
                p.fullCounts = this.counts(p.mapping);
            end
        end
        
        %partition state among processes
        function this = partition(this)
            %calculate requests
            touchs = zeros([size(this.counts) numel(this.partitions)]);
            reqs = zeros([size(this.counts) numel(this.partitions)]);
            for iPartition = 1:numel(this.partitions)
                p = this.partitions{iPartition};
                
                if ~p.isReqAbs
                    touch = zeros(size(this.counts));
                    touch(p.mapping) = 1;
                    touchs(:, :, iPartition) = touch;
                end
                
                req = zeros(size(this.counts));
                req(p.mapping) = max(0, p.reqFunc());
                reqs(:, :, iPartition) = req;
            end
            
            %partition
            tmp = cellfun(@(x) x.isReqAbs, this.partitions);
            absReqs = sum(reqs(:, :,  tmp), 3);
            relReqs = sum(reqs(:, :, ~tmp), 3);
            
            absScale = max(0, min(min(this.counts, absReqs) ./ absReqs, 1));
            relScale = max(0, max(0, this.counts - absReqs) ./ relReqs);
            relScale(relReqs == 0) = 0;
            
            unReqs = max(0, this.counts - absReqs) ./ sum(touchs, 3) .* (relReqs == 0);
            unReqs(sum(touchs, 3) == 0) = 0;
            
            for iPartition = 1:numel(this.partitions)
                p = this.partitions{iPartition};
                
                if p.isReqAbs
                    scale = absScale;
                else
                    scale = relScale;
                end
                
                alloc = floor(reqs(:, :, iPartition) .* scale + unReqs .* touchs(:, :, iPartition));
                this.partitionedCounts(:, :, iPartition) = alloc;
                p.counts = alloc(p.mapping);
            end
            
            %TODO: allocate unpartitioned molecules
            this.unpartitionedCounts = this.counts - sum(this.partitionedCounts, 3);
        end
        
        %merge sub-states partitioned to processes
        function this = merge(this)
            this.counts = this.unpartitionedCounts;
            for iPartition = 1:numel(this.partitions)
                p = this.partitions{iPartition};
                cnt = zeros(size(this.counts));
                cnt(p.mapping) = p.counts;
                this.counts = this.counts + cnt;
            end
        end
    end
    
    methods
        %get index of molecule by id (id, form, compartment)
        function [idxs, idFormIdxs, compIdxs] = getIndex(this, ids)
            if isempty(this.parentState)
                [idxs, idFormIdxs, compIdxs] = this.getIndex_parent(ids);
            else
                [~, idxs] = ismember(this.parentState.getIndex(ids), this.mapping);
                idFormIdxs = idxs;
                compIdxs = ones(size(idxs));
                if ~all(idxs)
                    throw(MException('MoleculeCounts:invalidId', 'Invalid ids:\n- %s', ...
                        strjoin(sprintf('\n- '), ids{idxs == 0})));
                end
            end
        end
        
        function [idxs, idFormIdxs, compIdxs] = getIndex_parent(this, ids)
            if ischar(ids)
                ids = {ids};
            end
            
            idForms = cell(numel(ids), 2);
            comps = cell(numel(ids), 1);
            for i = 1:numel(ids)
                match = regexp(ids{i}, '^(?<molecule>[^:\[\]]+)(?<form>:[^:\[\]]+)*(?<compartment>\[[^:\[\]]+\])*$', 'names');
                if isempty(match)
                    throw(MException('MoleculeCounts:invalidId', 'Invalid id: %s', ids{i}))
                end
                idForms{i, 1} = match.molecule;
                if isempty(match.form)
                    idForms{i, 2} = 1;
                else
                    idForms{i, 2} = this.formVals.(match.form(2:end));
                end
                if isempty(match.compartment)
                    comps{i, 1} = this.compartments(1).id;
                else
                    comps{i, 1} = match.compartment(2:end-1);
                end
            end
            [~, compIdxs] = ismember(comps, {this.compartments.id});
            if ~all(compIdxs)
                throw(MException('MoleculeCounts:invalidCompartment', 'Invalid compartment'))
            end
            
            [~, idFormIdxs] = ismember(...
                cellfun(@(x,y) [x ':' num2str(y)], idForms(:, 1), idForms(:, 2), 'UniformOutput', false), ...
                cellfun(@(x,y) [x ':' num2str(y)], this.ids, num2cell(this.forms), 'UniformOutput', false));
            if ~all(idFormIdxs)
                tmp = cellfun(@(x,y) [x ':' num2str(y)], ...
                    idForms(idFormIdxs == 0, 1), {this.compartments(cell2mat(idForms(idFormIdxs == 0, 2))).id}', ...
                    'UniformOutput', false);
                throw(MException('MoleculeCounts:invalidIdForm', 'Invalid id/form:\n- %s', ...
                    strjoin(sprintf('\n- '), tmp{:})));
            end
            
            idxs = sub2ind([numel(this.ids) numel(this.compartments)], idFormIdxs, compIdxs);
            idxs = reshape(idxs, size(ids));
        end
    end
end