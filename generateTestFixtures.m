%generateTestFixtures
% Generates fixtures for whole-cell tests.
%
% Example:
% >> setWarnings();
% >> setPath();
% >> generateTestFixtures();
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/8/2013
function generateTestFixtures()
import wholecell.kb.KnowledgeBase;
import wholecell.sim.Simulation;

%create output directory
outDir = 'data/fixtures';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%construct KB
kb = KnowledgeBase('dataFile', 'data/KnowledgeBase.xlsx', 'seqFile', 'data/KnowledgeBase.fna'); %construct
save(fullfile(outDir, 'KnowledgeBase.mat'), 'kb'); %save

%construct simulation
sim = Simulation(kb); %construct
sim.setOptions(struct('seed', 1)); %seed
sim.calcInitialConditions(); %calculate initial conditions
save(fullfile(outDir, 'Simulation.mat'), 'sim'); %save