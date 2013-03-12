%SimulationTest
% Tests whole-cell simulation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
classdef SimulationTest < TestCase
    %fixture
    properties (Access = protected)
        sim
    end
    
    %constructor
    methods
        function this = SimulationTest(name)
            this = this@TestCase(name);
        end
    end
    
    %setup, teardown
    methods
        function setUp(this) %#ok<*MANU>
            tmp = load('data/fixtures/Simulation.mat');
            this.sim = tmp.sim;
        end
        
        function tearDown(this)
            this.sim = [];
        end
    end
    
    %tests for run-time errors
    methods
        function testConstruction(this)
            import wholecell.sim.Simulation;
            
            %load cached KB
            tmp = load('data/fixtures/KnowledgeBase.mat');
            
            %construct simulation
            sim = Simulation(tmp.kb);
            sim.setOptions(struct('seed', 1));
            sim.calcInitialConditions();
        end
        
        function testRun(this)
            import wholecell.sim.logger.Shell;
            
            %simulate
            sim = this.sim; %#ok<*PROP>
            sim.setOptions(struct('lengthSec', 10));
            sim.run();
            
            %test time progressed
            assertEqual(10, sim.getState('Time').value);
        end
        
        %test logging with timeStepSec = 1 s
        function testLogging(this)
            import wholecell.sim.logger.Disk;
            import wholecell.sim.logger.Shell;
            
            %output directory
            outDir = fullfile(pwd, 'out/test/SimulationTest_testLogging');
            
            %run simulation
            sim = this.sim; %#ok<*PROP>
            sim.setOptions(struct('lengthSec', 25));
            sim.run({
                Shell()
                Disk('outDir', outDir, 'segmentLen', 5)
                });
            
            %retrieve data & check sizes and values
            time = Disk.load(outDir, 'Time', 'value');
            assertEqual(permute(0:25, [1 3 2]), time);
            
            assertEqual([1 3 26], size(Disk.load(outDir, 'Mass', 'cell')));
            assertEqual([size(sim.getState('MoleculeCounts').counts) 26], size(Disk.load(outDir, 'MoleculeCounts', 'counts')));
            assertEqual([1 1 26], size(Disk.load(outDir, 'Metabolism', 'growth')));
        end
        
        %test logging with timeStepSec = 5 s
        function testLogging2(this)
            import wholecell.sim.logger.Disk;
            import wholecell.sim.logger.Shell;
            
            %output directory
            outDir = fullfile(pwd, 'out/test/SimulationTest_testLogging2');
            
            %run simulation
            sim = this.sim; %#ok<*PROP>
            sim.setOptions(struct('lengthSec', 50, 'timeStepSec', 5));
            sim.run({
                Shell()
                Disk('outDir', outDir, 'segmentLen', 5)
                });
            
            %retrieve data & check sizes and values
            time = Disk.load(outDir, 'Time', 'value');
            assertEqual(permute(0:5:50, [1 3 2]), time);
            
            assertEqual([1 3 11], size(Disk.load(outDir, 'Mass', 'cell')));
            assertEqual([size(sim.getState('MoleculeCounts').counts) 11], size(Disk.load(outDir, 'MoleculeCounts', 'counts')));
            assertEqual([1 1 11], size(Disk.load(outDir, 'Metabolism', 'growth')));
        end
        
        function testRunSimulation(~)
            %knowledge base options
            kbOpts = {
                'dataFile', 'data/KnowledgeBase.xlsx'
                'seqFile', 'data/KnowledgeBase.fna'
                }';
            
            %simulation options
            simOpts = struct(...
                'lengthSec', 100 ...
                );
            
            %disk logger options
            diskOpts = {
                'outDir', fullfile(pwd, 'out/test/SimulationTest_testRunSimulation')
                'metadata', struct(...
                    'name', 'Test simulation', ...
                    'description', 'Test simulation', ...
                    'investigator', struct(...
                        'first', 'Jonathan', ...
                        'last', 'Karr', ...
                        'email', 'jkarr@stanford.edu', ...
                        'affiliation', 'Stanford University' ...
                        ), ...
                    'ip', '171.65.92.146' ...
                    )
                'segmentLen', 10
                }';
            
            %run simulation
            runSimulation('kbOpts', kbOpts(:), 'simOpts', simOpts, 'diskOpts', diskOpts(:));
        end
    end
    
    %test biology
    methods
        function testMetabolicNetwork(this)
            import wholecell.util.Constants;
            
            sim = this.sim; %#ok<*PROP>
            met = sim.getProcess('Metabolism');
            
            %% reaction stoichiometry
            %Ex 1
            rIdx = find(strcmp('Cls3', met.rxnIds));
            assertEqual(3, nnz(met.sMat(:, rIdx)))
            assertEqual(-2, met.sMat(met.metabolite.getIndex('PG181[m]'), rIdx));
            assertEqual( 1, met.sMat(met.metabolite.getIndex('CL181[m]'), rIdx));
            assertEqual( 1, met.sMat(met.metabolite.getIndex('GL[c]'), rIdx));
            
            %Ex 2
            rIdx = find(strcmp('HinT_GMP_Mor_MG132', met.rxnIds));
            assertEqual(5, nnz(met.sMat(:, rIdx)))
            assertEqual(-1, met.sMat(met.metabolite.getIndex('GMP_Mor[e]'), rIdx));
            assertEqual(-1, met.sMat(met.metabolite.getIndex('H2O[e]'), rIdx));
            assertEqual( 1, met.sMat(met.metabolite.getIndex('GMP[e]'), rIdx));
            assertEqual( 2, met.sMat(met.metabolite.getIndex('H[e]'), rIdx));
            assertEqual( 1, met.sMat(met.metabolite.getIndex('MOR[e]'), rIdx));
            
            %Ex 3
            rIdx = find(strcmp('MsrA', met.rxnIds));
            assertEqual(5, nnz(met.sMat(:, rIdx)))
            assertEqual(-1, met.sMat(met.metabolite.getIndex('H2O[c]'), rIdx));
            assertEqual(-1, met.sMat(met.metabolite.getIndex('MET[c]'), rIdx));
            assertEqual(-1, met.sMat(met.metabolite.getIndex('MG_124_MONOMER_ox[c]'), rIdx));
            assertEqual( 1, met.sMat(met.metabolite.getIndex('METSOXSL[c]'), rIdx));
            assertEqual( 1, met.sMat(met.metabolite.getIndex('MG_124_MONOMER[c]'), rIdx));
            
            %% exchange reactions
            assertEqual(eye(numel(met.metIdx.real)), met.sMat(met.metIdx.real, met.rxnIdx.exchange))
            assertFalse(any(any(met.sMat(setdiff(1:end, met.metIdx.real), met.rxnIdx.exchange))))
            
            %% objective
            %biomass composition
            assertElementsAlmostEqual(13.1, ...
                -met.sMat(met.metIdx.real, met.rxnIdx.growth)' * met.metabolite.mws / Constants.nAvogadro * 1e15, ...
                'relative', 1e-3);
            
            assertEqual(786393.089973227, -met.sMat(met.metabolite.getIndex('ALA[c]'), met.rxnIdx.growth));
            assertEqual(-23513832.1575775, -met.sMat(met.metabolite.getIndex('AMP[c]'), met.rxnIdx.growth));
            
            %objective
            assertEqual(1, nnz(met.objective))
            assertEqual(1, met.objective(met.rxnIdx.growth));
            
            assertEqual( 1, met.sMat(met.metIdx.biomass, met.rxnIdx.growth))
            assertEqual(-1, met.sMat(met.metIdx.biomass, met.rxnIdx.biomassExchange))
            assertEqual(1, nnz(met.sMat(:, met.rxnIdx.biomassExchange)));
            
            %% enzymes & kinetic bounds
            assertTrue(all(met.eMat(:) >=0))
            assertTrue(all(sum(met.eMat, 2) <= 1))
            assertTrue(all(isinf(met.bounds.kinetic.lo(~any(met.eMat, 2)))))
            assertTrue(all(isinf(met.bounds.kinetic.up(~any(met.eMat, 2)))))
            assertTrue(all(met.bounds.kinetic.lo <= 0))
            assertTrue(all(met.bounds.kinetic.up >= 0))
            
            %Ex 1
            iRxn = find(strcmp('LdhA', met.rxnIds));
            iEnz = met.enzyme.getIndex('MG_460_TETRAMER[c]');
            assertEqual(2000 / 60 * 1e-3 * met.enzyme.mws(iEnz), met.bounds.kinetic.up(iRxn));
            assertEqual(-inf, met.bounds.kinetic.lo(iRxn));
            
            %Ex 2
            iRxn = find(strcmp('MetF', met.rxnIds));
            iEnz = met.enzyme.getIndex('MG_228_TETRAMER[c]');
            assertEqual(1, met.eMat(iRxn, iEnz));
            assertEqual(inf, met.bounds.kinetic.up(iRxn));
            assertEqual(-9.7 / 60 * 1e-3 * met.enzyme.mws(iEnz), met.bounds.kinetic.lo(iRxn));
            
            %% exchange bounds
            assertTrue(all(met.bounds.exchange.lo <= 0))
            assertTrue(all(met.bounds.exchange.up >= 0))
            
            assertFalse(any(met.bounds.exchange.lo(met.rxnIdx.internalExchange)));
            assertFalse(any(met.bounds.exchange.up(met.rxnIdx.internalExchange)));
            
            iRxn = met.rxnIdx.exchange(met.metabolite.getIndex('AD[c]'));
            assertEqual(0, met.bounds.exchange.lo(iRxn))
            assertEqual(0, met.bounds.exchange.up(iRxn))
            
            iRxn = met.rxnIdx.exchange(met.metabolite.getIndex('AD[e]'));
            assertEqual(-12, met.bounds.exchange.lo(iRxn))
            assertEqual( 12, met.bounds.exchange.up(iRxn))
            
            iRxn = met.rxnIdx.exchange(met.metabolite.getIndex('H2O[e]'));
            assertEqual(-20, met.bounds.exchange.lo(iRxn))
            assertEqual( 20, met.bounds.exchange.up(iRxn))
            
            %% thermodynamics bounds
            assertTrue(all(met.bounds.thermodynamic.lo <= 0))
            assertTrue(all(met.bounds.thermodynamic.up >= 0))
            
            iRxn = find(strcmp(met.rxnIds, 'AckA'));
            assertEqual(-inf, met.bounds.thermodynamic.lo(iRxn))
            assertEqual( inf, met.bounds.thermodynamic.up(iRxn))
            
            iRxn = find(strcmp(met.rxnIds, 'AcpS'));
            assertEqual(0, met.bounds.thermodynamic.lo(iRxn))
            assertEqual( inf, met.bounds.thermodynamic.up(iRxn))
            
            %% growth rate
            state_met = sim.getState('Metabolism');
            growth = zeros(10, 1);
            for i = 1:numel(growth)
                sim.setOptions(struct('seed', i));
                sim.calcInitialConditions();
                growth(i) = state_met.growth;
            end
            assertElementsAlmostEqual(log(2) / (9 * 3600) * 3600 * 13.1, mean(growth), 'relative', 10e-2); %fg/h
        end
        
        function testInitialConditions(this)
            sim = this.sim; %#ok<*PROP>
            mass = sim.getState('Mass');
            mc = sim.getState('MoleculeCounts');
            met = sim.getState('Metabolism');
            
            %calculate initial conditions
            sim.setOptions(struct('seed', 1));
            sim.calcInitialConditions();
            
            %metabolite counts
            assertTrue(all(isfinite(mc.counts(:))));
            assertTrue(all(mc.counts(:) >= 0));
            
            %mass
            cellCompIdxs = [mass.cIdx.c; mass.cIdx.m];
            assertElementsAlmostEqual(13.1, sum(mass.cell(cellCompIdxs)), 'relative', 1e-2); %fg
            assertElementsAlmostEqual(13.1 * (0.7  + 0.3 * (1 - 0.0929 - 0.6197)), sum(mass.metabolite(cellCompIdxs)), 'relative', 1e-2); %fg
            assertElementsAlmostEqual(13.1 * 0.3 * 0.0929, sum(mass.rna(cellCompIdxs)), 'relative', 6e-1); %fg
            assertElementsAlmostEqual(13.1 * 0.3 * 0.6197, sum(mass.protein(cellCompIdxs)), 'relative', 1e-1); %fg
            
            %growth
            assertElementsAlmostEqual(log(2) / (9 * 3600) * 3600 * 13.1, met.growth, 'relative', 1e-1); %fg/h
        end
        
        function testGrowth(this)
            import wholecell.util.Constants;
            
            sim = this.sim;
            met = sim.getState('Metabolism');
            mc = sim.getState('MoleculeCounts');
            mass = sim.getState('Mass');
            
            %% simulate
            %set options
            sim.setOptions(struct('seed', 1, 'lengthSec', 15000, 'timeStepSec', 100));
            
            %calculate initial conditions
            sim.calcInitialConditions();
            
            %record initial state
            init = sim.getDynamics();
            init.Mass.matureRna       = mc.mws(mc.idx.matureRna)'       * sum(mc.counts(mc.idx.matureRna,       :), 2) / Constants.nAvogadro * 1e15;
            init.Mass.matureMonomers  = mc.mws(mc.idx.matureMonomers)'  * sum(mc.counts(mc.idx.matureMonomers,  :), 2) / Constants.nAvogadro * 1e15;
            init.Mass.matureComplexes = mc.mws(mc.idx.matureComplexes)' * sum(mc.counts(mc.idx.matureComplexes, :), 2) / Constants.nAvogadro * 1e15;
            
            %simulate dynamics
            time = sim.getState('Time');
            mc = sim.getState('MoleculeCounts');
            trl = sim.getProcess('Translation');
            pm = sim.getProcess('ProteinMaturation');
            rm = sim.getProcess('RnaMaturation');
            cpx = sim.getProcess('Complexation');
            
            matureRna = zeros(1, sim.lengthSec / sim.timeStepSec);
            matureMrna = zeros(1, sim.lengthSec / sim.timeStepSec);
            matureMonomer = zeros(1, sim.lengthSec / sim.timeStepSec);
            matureComplex = zeros(1, sim.lengthSec / sim.timeStepSec);
            for iSec = sim.timeStepSec:sim.timeStepSec:sim.lengthSec
                fprintf('Time = %d s\n', iSec);
                time.value = iSec;
                sim.evolveState();
                
                matureRna(iSec / sim.timeStepSec) = sum(mc.counts(rm.matureRna.mapping));
                matureMrna(iSec / sim.timeStepSec) = sum(mc.counts(trl.mrna.mapping));
                matureMonomer(iSec / sim.timeStepSec) = sum(mc.counts(pm.matureProteinMonomer.mapping));
                matureComplex(iSec / sim.timeStepSec) = sum(mc.counts(cpx.complex.mapping));
            end
            matureNoncodingRna = matureRna - matureMrna;
            
            %% plot
            time = sim.timeStepSec:sim.timeStepSec:sim.lengthSec;
            plot(time, matureNoncodingRna);
            plot(time, matureMrna);
            plot(time, matureMonomer);
            plot(time, matureComplex);
            
            %% assert
            expCumGrowth = exp(log(2) * iSec / 30000);
            
            %growth
            assertElementsAlmostEqual(expCumGrowth * init.Metabolism.growth, met.growth, 'relative', 2e-1);
            
            %mass
            assertElementsAlmostEqual(expCumGrowth * init.Mass.cell,       mass.cell,    'relative', 2e-1);
            %assertElementsAlmostEqual(expCumGrowth * init.Mass.cellDry,    mass.cellDry, 'relative', 2e-1);
            %assertElementsAlmostEqual(expCumGrowth * init.Mass.rna,        mass.rna,     'relative', 2e-1);
            assertElementsAlmostEqual(expCumGrowth * init.Mass.protein,    mass.protein, 'relative', 2e-1);
            
            cellCompIdxs = [mass.cIdx.c; mass.cIdx.m];
            assertElementsAlmostEqual(expCumGrowth * init.Mass.metabolite(cellCompIdxs), mass.metabolite(cellCompIdxs), 'relative', 2e-1);
            
            %physical counts
            %assertTrue(all(isfinite(mc.counts(:)) & mc.counts(:) >= 0)); %TODO
            
            %RNA, protein matured and complexed
            assertElementsAlmostEqual(expCumGrowth * init.Mass.matureRna, ...
                mc.mws(mc.idx.matureRna)' * sum(mc.counts(mc.idx.matureRna, :), 2) / Constants.nAvogadro * 1e15, ...
                'relative', 2e-1);
            assertElementsAlmostEqual(expCumGrowth * init.Mass.matureMonomers, ...
                mc.mws(mc.idx.matureMonomers)' * sum(mc.counts(mc.idx.matureMonomers, :), 2) / Constants.nAvogadro * 1e15, ...
                'relative', 2e-1);
            assertElementsAlmostEqual(expCumGrowth * init.Mass.matureComplexes, ...
                mc.mws(mc.idx.matureComplexes)' * sum(mc.counts(mc.idx.matureComplexes, :), 2) / Constants.nAvogadro * 1e15, ...
                'relative', 2e-1);
            
            %no mature RNA, protein in wrong compartments
            matIdxs = find(mc.types == mc.typeVals.rna & mc.forms == mc.formVals.mature);
            matCompIdxs = sub2ind(size(mc.counts), matIdxs(:, ones(3, 1)), repmat(1:3, size(matIdxs)));
            assertFalse(any(mc.counts(setdiff(matCompIdxs(:), mc.idx.matureRna))));
            
            matIdxs = find(mc.types == mc.typeVals.protein & mc.forms == mc.formVals.mature);
            matAllCompIdxs = sub2ind(size(mc.counts), matIdxs(:, ones(3, 1)), repmat(1:3, size(matIdxs)));
            matCompIdxs = sub2ind(size(mc.counts), matIdxs, mc.localizations(matIdxs));
            assertFalse(any(mc.counts(setdiff(matAllCompIdxs(:), matCompIdxs))));
        end
    end
end