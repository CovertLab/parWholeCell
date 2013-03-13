%KnowledgeBaseTest
% Tests whole-cell knowledge base.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/8/2013
classdef KnowledgeBaseTest < TestCase
    %constructor
    properties
        kb
    end
    
    methods
        function this = KnowledgeBaseTest(name)
            this = this@TestCase(name);
        end
        
        function setUp(this)
            tmp = load('data/fixtures/KnowledgeBase.mat');
            this.kb = tmp.kb;
        end
        
        function tearDown(this)
            this.kb = [];
        end
    end
    
    %tests
    methods
        function testConstruction(~)
            import wholecell.kb.KnowledgeBase;
            
            KnowledgeBase('dataFile', 'data/KnowledgeBase.xlsx', 'seqFile', 'data/KnowledgeBase.fna');
        end
        
        function testMetabolites(this)
            kb = this.kb; %#ok<*PROP>
            
            assertEqual(722, numel(kb.metabolites));
            
            assertEqual(82, sum([kb.metabolites.hydrophobic]));
            assertEqual(640, sum(~[kb.metabolites.hydrophobic]));
            
            met = kb.metabolites(strcmp({kb.metabolites.id}, 'ACCOA'));
            assertTrue(isstruct(met));
            assertEqual('ACCOA', met.id);
            assertEqual('Acetyl-CoA', met.name);
            assertEqual('C23H34N7O17P3S1', met.formula);
            assertEqual('CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)N1C=NC2=C1N=CN=C2N', met.smiles);
            assertEqual(-4, met.charge);
            assertEqual(805.538, met.mw);
            assertEqual(false, met.hydrophobic);
            assertEqual(0, met.mediaConc);
            assertElementsAlmostEqual(3.5812e+03, met.biomassConc, 'relative', 1e-4);
            assertElementsAlmostEqual(3.5812e+03, met.metabolismNewFlux, 'relative', 1e-4);
            assertElementsAlmostEqual(0, met.metabolismRecyclingFlux, 'relative', 1e-4);
            
            met = kb.metabolites(strcmp({kb.metabolites.id}, 'AC'));
            assertElementsAlmostEqual(0.304753, met.mediaConc, 'relative', 1e-4);
        end
        
        function testGeneticCode(this)
            kb = this.kb;
            
            assertEqual(4, kb.translationTable);
        end
        
        function testGenome(this)
            kb = this.kb;
            
            assertEqual([1 580076], size(kb.genomeSeq));
            assertEqual('ACGT', unique(kb.genomeSeq));
        end
        
        function testGenes(this)
            kb = this.kb;
            
            assertEqual(525, numel(kb.genes));
            assertEqual(482, sum(strcmp({kb.genes.type}, 'mRNA')));
            assertEqual(3, sum(strcmp({kb.genes.type}, 'rRNA')));
            assertEqual(4, sum(strcmp({kb.genes.type}, 'sRNA')));
            assertEqual(36, sum(strcmp({kb.genes.type}, 'tRNA')));
            
            gene = kb.genes(strcmp({kb.genes.id}, 'MG_001'));
            assertTrue(isstruct(gene));
            assertEqual('MG_001', gene.id);
            assertEqual('DNA polymerase III, beta subunit', gene.name);
            assertEqual('dnaN', gene.symbol);
            assertEqual('mRNA', gene.type);
            assertEqual(686, gene.start);
            assertEqual(1143, gene.len);
            assertEqual(true, gene.dir);
            assertEqual([
                'ATGAAAATATTAATTAATAAAAGTGAATTGAATAAAATTTTGAAAAAAAT' ...
                'GAATAACGTTATTATTTCCAATAACAAAATAAAACCACATCATTCATATT' ...
                'TTTTAATAGAGGCAAAAGAAAAAGAAATAAACTTTTATGCTAACAATGAA' ...
                'TACTTTTCTGTCAAATGTAATTTAAATAAAAATATTGATATTCTTGAACA' ...
                'AGGCTCCTTAATTGTTAAAGGAAAAATTTTTAACGATCTTATTAATGGCA' ...
                'TAAAAGAAGAGATTATTACTATTCAAGAAAAAGATCAAACACTTTTGGTT' ...
                'AAAACAAAAAAAACAAGTATTAATTTAAACACAATTAATGTGAATGAATT' ...
                'TCCAAGAATAAGGTTTAATGAAAAAAACGATTTAAGTGAATTTAATCAAT' ...
                'TCAAAATAAATTATTCACTTTTAGTAAAAGGCATTAAAAAAATTTTTCAC' ...
                'TCAGTTTCAAATAATCGTGAAATATCTTCTAAATTTAATGGAGTAAATTT' ...
                'CAATGGATCCAATGGAAAAGAAATATTTTTAGAAGCTTCTGACACTTATA' ...
                'AACTATCTGTTTTTGAGATAAAGCAAGAAACAGAACCATTTGATTTCATT' ...
                'TTGGAGAGTAATTTACTTAGTTTCATTAATTCTTTTAATCCTGAAGAAGA' ...
                'TAAATCTATTGTTTTTTATTACAGAAAAGATAATAAAGATAGCTTTAGTA' ...
                'CAGAAATGTTGATTTCAATGGATAACTTTATGATTAGTTACACATCGGTT' ...
                'AATGAAAAATTTCCAGAGGTAAACTACTTTTTTGAATTTGAACCTGAAAC' ...
                'TAAAATAGTTGTTCAAAAAAATGAATTAAAAGATGCACTTCAAAGAATTC' ...
                'AAACTTTGGCTCAAAATGAAAGAACTTTTTTATGCGATATGCAAATTAAC' ...
                'AGTTCTGAATTAAAAATAAGAGCTATTGTTAATAATATCGGAAATTCTCT' ...
                'TGAGGAAATTTCTTGTCTTAAATTTGAAGGTTATAAACTTAATATTTCTT' ...
                'TTAACCCAAGTTCTCTATTAGATCACATAGAGTCTTTTGAATCAAATGAA' ...
                'ATAAATTTTGATTTCCAAGGAAATAGTAAGTATTTTTTGATAACCTCTAA' ...
                'AAGTGAACCTGAACTTAAGCAAATATTGGTTCCTTCAAGATAA'
                ], gene.seq);
            assertEqual('MG_001', gene.rnaId);
        end
        
        function testRnas(this)
            kb = this.kb;
            
            assertEqual(525, numel(kb.rnas));
            
            rna = kb.rnas(strcmp({kb.rnas.id}, 'MG_001'));
            assertTrue(isstruct(rna));
            assertEqual('MG_001', rna.id);
            assertEqual('DNA polymerase III, beta subunit', rna.name);
            assertEqual('mRNA', rna.type);
            assertElementsAlmostEqual(8.8983e-05, rna.exp, 'relative', 1e-4);
            assertElementsAlmostEqual(146.9388, rna.halfLife, 'relative', 1e-4);
            assertEqual([
                'UACUUUUAUAAUUAAUUAUUUUCACUUAACUUAUUUUAAAACUUUUUUUA' ...
                'CUUAUUGCAAUAAUAAAGGUUAUUGUUUUAUUUUGGUGUAGUAAGUAUAA' ...
                'AAAAUUAUCUCCGUUUUCUUUUUCUUUAUUUGAAAAUACGAUUGUUACUU' ...
                'AUGAAAAGACAGUUUACAUUAAAUUUAUUUUUAUAACUAUAAGAACUUGU' ...
                'UCCGAGGAAUUAACAAUUUCCUUUUUAAAAAUUGCUAGAAUAAUUACCGU' ...
                'AUUUUCUUCUCUAAUAAUGAUAAGUUCUUUUUCUAGUUUGUGAAAACCAA' ...
                'UUUUGUUUUUUUUGUUCAUAAUUAAAUUUGUGUUAAUUACACUUACUUAA' ...
                'AGGUUCUUAUUCCAAAUUACUUUUUUUGCUAAAUUCACUUAAAUUAGUUA' ...
                'AGUUUUAUUUAAUAAGUGAAAAUCAUUUUCCGUAAUUUUUUUAAAAAGUG' ...
                'AGUCAAAGUUUAUUAGCACUUUAUAGAAGAUUUAAAUUACCUCAUUUAAA' ...
                'GUUACCUAGGUUACCUUUUCUUUAUAAAAAUCUUCGAAGACUGUGAAUAU' ...
                'UUGAUAGACAAAAACUCUAUUUCGUUCUUUGUCUUGGUAAACUAAAGUAA' ...
                'AACCUCUCAUUAAAUGAAUCAAAGUAAUUAAGAAAAUUAGGACUUCUUCU' ...
                'AUUUAGAUAACAAAAAAUAAUGUCUUUUCUAUUAUUUCUAUCGAAAUCAU' ...
                'GUCUUUACAACUAAAGUUACCUAUUGAAAUACUAAUCAAUGUGUAGCCAA' ...
                'UUACUUUUUAAAGGUCUCCAUUUGAUGAAAAAACUUAAACUUGGACUUUG' ...
                'AUUUUAUCAACAAGUUUUUUUACUUAAUUUUCUACGUGAAGUUUCUUAAG' ...
                'UUUGAAACCGAGUUUUACUUUCUUGAAAAAAUACGCUAUACGUUUAAUUG' ...
                'UCAAGACUUAAUUUUUAUUCUCGAUAACAAUUAUUAUAGCCUUUAAGAGA' ...
                'ACUCCUUUAAAGAACAGAAUUUAAACUUCCAAUAUUUGAAUUAUAAAGAA' ...
                'AAUUGGGUUCAAGAGAUAAUCUAGUGUAUCUCAGAAAACUUAGUUUACUU' ...
                'UAUUUAAAACUAAAGGUUCCUUUAUCAUUCAUAAAAAACUAUUGGAGAUU' ...
                'UUCACUUGGACUUGAAUUCGUUUAUAACCAAGGAAGUUCUAUU'
                ], rna.seq);
            assertEqual([sum(rna.seq == 'A'); sum(rna.seq == 'C'); sum(rna.seq == 'G'); sum(rna.seq == 'U')], rna.ntCount);
            assertElementsAlmostEqual(362605.040400, rna.mw, 'relative', 1e-5);
            assertEqual('MG_001', rna.geneId);
            assertEqual('MG_001_MONOMER', rna.monomerId);
        end
        
        function testProteins(this)
            kb = this.kb;
            
            assertEqual(482 + 164, numel(kb.proteins));
            assertEqual(482, sum([kb.proteins.monomer]));
            assertEqual(164, sum(~[kb.proteins.monomer]));
            
            %monomers
            mon = kb.proteins(strcmp({kb.proteins.id}, 'MG_001_MONOMER'));
            assertTrue(isstruct(mon));
            assertEqual('DNA polymerase III, beta subunit', mon.name);
            assertTrue(mon.monomer);
            assertTrue(isempty(mon.composition));
            assertEqual('c', mon.compartment);
            assertTrue(isempty(mon.formationProcess))
            assertEqual([
                'MKILINKSELNKILKKMNNVIISNNKIKPHHSYFLIEAKEKEINFYANNE' ...
                'YFSVKCNLNKNIDILEQGSLIVKGKIFNDLINGIKEEIITIQEKDQTLLV' ...
                'KTKKTSINLNTINVNEFPRIRFNEKNDLSEFNQFKINYSLLVKGIKKIFH' ...
                'SVSNNREISSKFNGVNFNGSNGKEIFLEASDTYKLSVFEIKQETEPFDFI' ...
                'LESNLLSFINSFNPEEDKSIVFYYRKDNKDSFSTEMLISMDNFMISYTSV' ...
                'NEKFPEVNYFFEFEPETKIVVQKNELKDALQRIQTLAQNERTFLCDMQIN' ...
                'SSELKIRAIVNNIGNSLEEISCLKFEGYKLNISFNPSSLLDHIESFESNE' ...
                'INFDFQGNSKYFLITSKSEPELKQILVPSR'
                ], mon.seq);
            assertEqual([
                sum(mon.seq == 'A'); sum(mon.seq == 'R'); sum(mon.seq == 'N'); sum(mon.seq == 'D'); sum(mon.seq == 'C');
                sum(mon.seq == 'E'); sum(mon.seq == 'Q'); sum(mon.seq == 'G'); sum(mon.seq == 'H'); sum(mon.seq == 'I');
                sum(mon.seq == 'L'); sum(mon.seq == 'K'); sum(mon.seq == 'M'); sum(mon.seq == 'F'); sum(mon.seq == 'P');
                sum(mon.seq == 'S'); sum(mon.seq == 'T'); sum(mon.seq == 'W'); sum(mon.seq == 'Y'); sum(mon.seq == 'V');
                ], mon.aaCount);
            assertEqual(zeros(4, 1), mon.ntCount);
            assertElementsAlmostEqual(44317.8348 - (177.22 - 149.21), mon.mw, 'relative', 1e-4);
            assertEqual('MG_001', mon.geneId);
            assertEqual('MG_001', mon.geneId);
            
            %complexes
            cpx = kb.proteins(strcmp({kb.proteins.id}, 'DNA_GYRASE'));
            assertTrue(isstruct(cpx));
            assertEqual('DNA_GYRASE', cpx.id);
            assertEqual('DNA gyrase', cpx.name);
            assertFalse(cpx.monomer);
            assertEqual([
                struct('molecule', 'MG_003_MONOMER', 'compartment', 'c', 'coeff', -2, 'form', 'mature', 'type', 'protein')
                struct('molecule', 'MG_004_MONOMER', 'compartment', 'c', 'coeff', -2, 'form', 'mature', 'type', 'protein')
                struct('molecule', 'DNA_GYRASE', 'compartment', 'c', 'coeff', 1, 'form', 'mature', 'type', 'protein')
                ], cpx.composition);
            assertEqual('c', cpx.compartment);
            assertEqual('Complexation', cpx.formationProcess);
            assertTrue(isempty(cpx.seq));
            assertEqual(...
                + 2 * kb.proteins(strcmp({kb.proteins.id}, 'MG_003_MONOMER')).aaCount ...
                + 2 * kb.proteins(strcmp({kb.proteins.id}, 'MG_004_MONOMER')).aaCount ...
                , cpx.aaCount);
            assertEqual(zeros(4, 1), cpx.ntCount);
            assertElementsAlmostEqual(334028.2216, cpx.mw, 'relative', 1e-2);
            assertTrue(isempty(cpx.geneId));
            assertTrue(isempty(cpx.rnaId));
            
            cpx = kb.proteins(strcmp({kb.proteins.id}, 'MG_014_015_DIMER'));
            assertEqual('m', cpx.compartment);
        end
        
        function testReactions(this)
            kb = this.kb;
            
            assertEqual(643, numel(kb.reactions));
            
            rxn = kb.reactions(strcmp({kb.reactions.id}, 'AckA'));
            assertTrue(isstruct(rxn));
            assertEqual('AckA', rxn.id);
            assertEqual('acetate kinase', rxn.name);
            assertEqual('Metabolism', rxn.process);
            assertEqual('2.7.2.1', rxn.ec);
            assertEqual(0, rxn.dir);
            assertEqual([
                struct('molecule', 'ACTP', 'form', 'mature', 'compartment', 'c', 'coeff', -1, 'type', 'metabolite')
                struct('molecule', 'ADP', 'form', 'mature', 'compartment', 'c', 'coeff', -1, 'type', 'metabolite')
                struct('molecule', 'AC', 'form', 'mature', 'compartment', 'c', 'coeff', 1, 'type', 'metabolite')
                struct('molecule', 'ATP', 'form', 'mature', 'compartment', 'c', 'coeff', 1, 'type', 'metabolite')
                ], rxn.stoichiometry)
            assertEqual(struct(...
                'id', 'MG_357_DIMER', ...
                'form', 'mature', ...
                'compartment', 'c', ...
                'kCatFor', 68 / 60 * 1e-3 * kb.proteins(strcmp({kb.proteins.id}, 'MG_357_DIMER')).mw, ...
                'kCatRev', 70 / 60 * 1e-3 * kb.proteins(strcmp({kb.proteins.id}, 'MG_357_DIMER')).mw ...
                ), rxn.enzyme)
            
            
            rxn = kb.reactions(strcmp({kb.reactions.id}, 'Aas1'));
            assertEqual(1, rxn.dir)
            
            rxn = kb.reactions(strcmp({kb.reactions.id}, 'Cls1'));
            assertEqual([
                struct('molecule', 'PG160', 'form', 'mature', 'compartment', 'm', 'coeff', -2, 'type', 'metabolite')
                struct('molecule', 'CL160', 'form', 'mature', 'compartment', 'm', 'coeff', 1, 'type', 'metabolite')
                struct('molecule', 'GL', 'form', 'mature', 'compartment', 'c', 'coeff', 1, 'type', 'metabolite')
                ], rxn.stoichiometry)
        end
    end
end