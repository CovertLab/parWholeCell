%KnowledgeBase
% Whole-cell knowledge base.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/4/2013
classdef KnowledgeBase < handle
    properties (Access = protected)
        dataFile = 'data/KnowledgeBase.xlsx'
        seqFile = 'data/KnowledgeBase.fna'
    end
    
    properties
        metabolites
        translationTable
        genomeSeq
        genes
        rnas
        proteins
        reactions
    end
    
    methods
        function this = KnowledgeBase(varargin)
            %% parse input arguments
            ip = inputParser;
            ip.addParamValue('dataFile', 'data/KnowledgeBase.xlsx', @(x) exist(x, 'file'));
            ip.addParamValue('seqFile', 'data/KnowledgeBase.fna', @(x) exist(x, 'file'));
            ip.parse(varargin{:});
            
            this.dataFile = ip.Results.dataFile;
            this.seqFile = ip.Results.seqFile;
        
            %% parse data
            this.loadMetabolites();
            this.loadGenome();
            this.loadGenes();
            this.loadComplexes();
            this.loadReactions();
        end
    end
    
    methods (Access = protected)
        function this = loadMetabolites(this)
            [~, ~, raw] = xlsread(this.dataFile, 'Metabolites');
            this.metabolites = [];
            for i = 2:size(raw, 1)
                if all(cellfun(@(x) isnumeric(x) && isnan(x), raw(i, :)))
                    continue;
                end
                
                m = struct(...
                    'id', raw{i, 1}, ...
                    'name', raw{i, 2}, ...
                    'formula', raw{i, 3}, ...
                    'smiles', raw{i, 4}, ...
                    'charge', raw{i, 5}, ...
                    'mw', raw{i, 6}, ...
                    'hydrophobic', strcmp(raw{i, 7}, 'Y'), ...
                    'mediaConc', 0, ...
                    'biomassConc', 0, ...
                    'metabolismFlux', 0, ...
                    'maxExchangeRate', 0 ... %mmol/gDCW/h
                    );
                if isnumeric(m.name) && isnan(m.name)
                    m.name = '';
                end
                if isnumeric(m.smiles) && isnan(m.smiles)
                    m.smiles = '';
                end
                if ~isnan(raw{i, 8})
                    m.mediaConc = raw{i, 8};
                end
                if ~isnan(raw{i, 9})
                    m.biomassConc = raw{i, 9};
                end
                if ~isnan(raw{i, 10})
                    m.metabolismFlux = raw{i, 10};
                end
                if ~isnan(raw{i, 11})
                    m.maxExchangeRate = raw{i, 11};
                end
                this.metabolites = [this.metabolites; m];
            end
        end
        
        function this = loadGenome(this)
            this.translationTable = 4;
            [~, this.genomeSeq] = fastaread(this.seqFile);
        end
        
        function this = loadGenes(this)
            [~, ~, raw] = xlsread(this.dataFile, 'Genes');
            this.genes = [];
            this.rnas = [];
            this.proteins = [];
            for i = 2:size(raw, 1)
                if all(cellfun(@(x) isnumeric(x) && isnan(x), raw(i, :)))
                    continue;
                end
                
                %gene
                g = struct(...
                    'id', raw{i, 1}, ...
                    'name', raw{i, 2}, ...
                    'symbol', raw{i, 3}, ...
                    'type', raw{i, 4}, ...
                    'start', raw{i, 5}, ...
                    'len', raw{i, 6}, ...
                    'dir', strcmp(raw{i, 7}, 'forward'), ...
                    'seq', [], ...
                    'rnaId', raw{i, 1} ...
                    );
                if isnumeric(g.name) && isnan(g.name)
                    g.name = '';
                end
                if isnumeric(g.symbol) && isnan(g.symbol)
                    g.symbol = '';
                end
                g.seq = this.genomeSeq(1, g.start + (1:g.len) - 1);
                if ~g.dir
                    g.seq = seqrcomplement(g.seq);
                end
                this.genes = [this.genes; g];
                
                %RNA
                r = struct(...
                    'id', g.id, ...
                    'name', g.name, ...
                    'type', g.type, ...
                    'exp', raw{i, 8}, ...
                    'halfLife', raw{i, 9}, ... %s
                    'seq', [], ...
                    'ntCount', [], ...
                    'mw', [], ...
                    'geneId', g.id, ...
                    'monomerId', [g.id '_MONOMER'] ...
                    );
                r.seq = dna2rna(seqcomplement(g.seq));
                tmp = basecount(r.seq);
                r.ntCount = [tmp.A; tmp.C; tmp.G; tmp.T];
                r.mw = ...
                    + r.ntCount' * [345.20; 321.18; 361.20; 322.17] ...
                    - (numel(r.seq) - 1) * 17.01;
                this.rnas = [this.rnas; r];
                
                %protein monomer
                if strcmp(g.type, 'mRNA')
                    p = struct(...
                        'id', [g.id '_MONOMER'], ...
                        'name', g.name, ...
                        'monomer', true, ...
                        'composition', [], ...
                        'compartment', raw{i,10}, ...
                        'formationProcess', '', ...
                        'seq', '', ...
                        'aaCount', zeros(20, 1), ...
                        'ntCount', zeros(4, 1), ...
                        'mw', 0, ...
                        'geneId', g.id, ...
                        'rnaId', g.id ...
                        );
                    p.seq = nt2aa(g.seq, 'GeneticCode', this.translationTable);
                    if any(p.seq == '*')                
                        p.seq = p.seq(1:strfind(p.seq, '*') - 1);
                    end
                    
                    tmp = aacount(p.seq);
                    p.aaCount = [tmp.A; tmp.R; tmp.N; tmp.D; tmp.C; tmp.E; tmp.Q; tmp.G; tmp.H; tmp.I; tmp.L; tmp.K; tmp.M; tmp.F; tmp.P; tmp.S; tmp.T; tmp.W; tmp.Y; tmp.V];
                    
                    p.mw = molweight(p.seq);
                    this.proteins = [this.proteins; p];
                end
            end
        end
        
        function this = loadComplexes(this)
            [~, ~, raw] = xlsread(this.dataFile, 'Complexes');            
            for i = 2:size(raw, 1)
                if all(cellfun(@(x) isnumeric(x) && isnan(x), raw(i, :)))
                    continue;
                end
                
                p = struct(...
                    'id', raw{i, 1}, ...
                    'name', raw{i, 2}, ...
                    'monomer', false, ...
                    'composition', [], ...
                    'compartment', '', ...
                    'formationProcess', raw{i, 4}, ...
                    'seq', '', ...                    
                    'aaCount', zeros(20, 1), ...
                    'ntCount', zeros(4, 1), ...
                    'mw', 0, ...
                    'geneId', [], ...
                    'rnaId', [] ...
                    );
                if isnumeric(p.name) && isnan(p.name)
                    p.name = '';
                end
                
                this.proteins = [this.proteins; p];
            end
            
            metIds = {this.metabolites.id};
            rnaIds = {this.rnas.id};
            protIds = {this.proteins.id};
            for i = 2:size(raw, 1)
                if all(cellfun(@(x) isnumeric(x) && isnan(x), raw(i, :)))
                    continue;
                end
                
                iCpx = find(strcmp(protIds, raw{i, 1}), 1, 'first');
                p = this.proteins(iCpx);
                
                p.composition = this.parseReaction(raw{i, 3});
                for j = 1:numel(p.composition)
                    if strcmp(p.id, p.composition(j).molecule)
                        p.compartment = p.composition(j).compartment;
                    else
                        metIdx = find(strcmp(p.composition(j).molecule, metIds), 1, 'first');
                        rnaIdx = find(strcmp(p.composition(j).molecule, rnaIds), 1, 'first');
                        protIdx = find(strcmp(p.composition(j).molecule, protIds(1:iCpx)), 1, 'first');
                        if ~isempty(metIdx)
                            subunitMw = this.metabolites(metIdx).mw;
                            %TODO: p.metCnt
                        elseif ~isempty(rnaIdx)
                            subunitMw = this.rnas(rnaIdx).mw;
                            p.ntCount = p.ntCount - p.composition(j).coeff * this.rnas(rnaIdx).ntCount;
                        elseif ~isempty(protIdx)
                            subunitMw = this.proteins(protIdx).mw;
                            p.aaCount = p.aaCount - p.composition(j).coeff * this.proteins(protIdx).aaCount;
                        else
                            throw(MException('KnowledgeBase:undefinedSubunit', 'Undefined subunit: %s', p.composition(j).molecule));
                        end
                        p.mw = p.mw - p.composition(j).coeff * subunitMw;
                    end
                end
                
                this.proteins(iCpx) = p;
            end
        end
        
        function this = loadReactions(this)
            [~, ~, raw] = xlsread(this.dataFile, 'Reactions');
            this.reactions = [];
            for i = 2:size(raw, 1)
                if all(cellfun(@(x) isnumeric(x) && isnan(x), raw(i, :)))
                    continue;
                end
                
                r = struct(...
                    'id', raw{i, 1}, ...
                    'name', raw{i, 2}, ...
                    'process', raw{i, 3}, ...
                    'ec', raw{i, 4}, ...
                    'dir', [], ...
                    'stoichiometry', [], ...
                    'enzyme', [] ...
                    );
                if isnumeric(r.name) && isnan(r.name)
                    r.name = '';
                end
                if isnumeric(r.ec) && isnan(r.ec)
                    r.ec = '';
                end
                [r.stoichiometry, r.dir] = this.parseReaction(raw{i, 5});
                if ~(isnumeric(raw{i, 6}) && isnan(raw{i, 6}))
                    r.enzyme = struct('id', [], 'compartment', [], 'kCatFor', [], 'kCatRev', []);
                    match = regexp(raw{i, 6}, '^(?<id>[^:\[\]]+)(?<form>:[^:\[\]]+)*(?<compartment>\[[^:\[\]]+\])*', 'names');
                    r.enzyme.id = match.id;
                    if isempty(match.form)
                        r.enzyme.form = 'mature';
                    else
                        r.enzyme.form = match.form;
                    end
                    r.enzyme.compartment = match.compartment(2:end-1);
                    r.enzyme.kCatFor = this.calcKCat(r.enzyme.id, raw{i, 7}, raw{i, 8});
                    r.enzyme.kCatRev = this.calcKCat(r.enzyme.id, raw{i, 9}, raw{i, 10});
                end
                this.reactions = [this.reactions; r];
            end
        end
        
        function [stoich, dir] = parseReaction(this, str)
            %detect if global compartment
            match = regexp(str, '^\[(?<comp>.*?)\]: (?<stoich>.*)$', 'names');
            if ~isempty(match)
                globalComp = match.comp;
                stoich = match.stoich;
            else
                globalComp = [];
                stoich = str;
            end
            
            %divide equation into left, right, direction
            match = regexp(stoich, '^(?<lefts>.*) (?<dir><*==>*) (?<rights>.*)$', 'names');
            if isempty(match)
                throw(MException('KnowledgeBase:invalidStoichiometry', 'Invalid stoichiometry: %s', stoich))
            end
            
            %direction
            switch match.dir
                case '==>', dir = 1;
                case '<==', dir = -1;
                case '<==>', dir = 0;
            end
            
            %stoichiometry
            stoich = [];
            
            lefts = strsplit(' + ', match.lefts);
            for i = 1:numel(lefts)
                [coeff, mol, form, comp, type] = this.parseReactionComponent(lefts{i}, globalComp);
                stoich = [stoich; struct('coeff', -coeff, 'compartment', comp, 'molecule', mol, 'form', form, 'type', type)]; %#ok<AGROW>
            end
            
            rights = strsplit(' + ', match.rights);
            for i = 1:numel(rights)
                [coeff, mol, form, comp, type] = this.parseReactionComponent(rights{i}, globalComp);
                stoich = [stoich; struct('coeff', coeff, 'compartment', comp, 'molecule', mol, 'form', form, 'type', type)]; %#ok<AGROW>
            end
        end
        
        function [coeff, mol, form, comp, type] = parseReactionComponent(this, str, globalComp)
            if isempty(globalComp)
                tmp = regexp(str, '^(?<coeff>\(\d*\.*\d*\) )*(?<mol>.+?)(?<form>:.+)*\[(?<comp>.+)\]$', 'names');
                if isempty(tmp)
                    throw(MException('KnowledgeBase:invalidStoichiometry', 'Invalid stoichiometry: %s', str))
                end
                
                if isempty(tmp.coeff)
                    coeff = 1;
                else
                    coeff = str2double(tmp.coeff(2:end-2));
                end
                
                mol = tmp.mol;
                
                if isempty(tmp.form)
                    form = 'mature';
                else
                    form = tmp.form(2:end);
                end
                
                comp = tmp.comp;
            else
                tmp = regexp(str, '^(?<coeff>\(\d*\.*\d*\) )*(?<mol>.+?)(?<form>:.+)*$', 'names');
                if isempty(tmp)
                    throw(MException('KnowledgeBase:invalidStoichiometry', 'Invalid stoichiometry: %s', str))
                end
                if isempty(tmp.coeff)
                    coeff = 1;
                else
                    coeff = str2double(tmp.coeff(2:end-2));
                end
                
                mol = tmp.mol;
                
                if isempty(tmp.form)
                    form = 'mature';
                else
                    form = tmp.form(2:end);
                end
                
                comp = globalComp;
            end
            
            if any(strcmp(mol, {this.metabolites.id}))
                type = 'metabolite';
            elseif any(strcmp(mol, {this.rnas.id}))
                type = 'rna';
            elseif any(strcmp(mol, {this.proteins.id}))
                type = 'protein';
            else
                throw(MException('KnowledgeBase:undefinedMolecule', 'Undefined molecule: %s', mol));
            end
        end
        
        function val = calcKCat(this, enzId, vMax, units)
            if (isnumeric(enzId) && isnan(enzId)) || isnan(vMax)
                val = NaN;
                return;
            end
            
            switch units
                case 'U/mg'
                    idx = find(strcmp(enzId, {this.proteins.id}), 1, 'first');
                    if isempty(idx)
                        throw(MException('KnowledgeBase:undefinedEnzyme', 'Undefined enzyme: %s', enzId));
                    end
                    val = vMax / 60 * 1e-3 * this.proteins(idx).mw;
                case '1/min'
                    val = vMax / 60;
                otherwise
                    throw(MException('KnowledgeBase:invalidKCatUnits', 'Invalid kCat units: %s', units));
            end
        end
    end
end