%runSimulation
% Runs and logs whole-cell simulation.
%
% Example:
% >> setWarnings();
% >> setPath();
% >> runSimulation()
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
function runSimulation(varargin)
%% import classes
import wholecell.sim.logger.Disk;
import wholecell.sim.logger.Shell;
import wholecell.sim.Simulation;

%% parse inputs
ip = inputParser;
ip.addParamValue('simOpts', struct());
ip.addParamValue('diskOpts', {}, @(x) iscell(x) && iseven(numel(x)));
ip.addParamValue('kbOpts', {}, @(x) iscell(x) && iseven(numel(x)));
ip.parse(varargin{:});

diskOpts = ip.Results.diskOpts;
simOpts = ip.Results.simOpts;
kbOpts = ip.Results.kbOpts;

%% instantiate knowledge base
kb = wholecell.kb.KnowledgeBase(kbOpts{:});

%% instantiate loggers
diskLogger = Disk(diskOpts{:});
shellLogger = Shell();
loggers = {
    diskLogger
    shellLogger
    };

%% run simulation
sim = Simulation(kb);
sim.setOptions(simOpts);
sim.run(loggers);

%% plot result
%load data
time = permute(Disk.load(diskLogger.outDir, 'Time', 'value') / 3600, [1 3 2]);
mass = permute(sum(Disk.load(diskLogger.outDir, 'Mass', 'cell'), 2), [1 3 2]);

%plot
plot(time, mass);
xlabel('Time (h)');
ylabel('Mass (fg)');
xlim(time([1 end]));
line(time([1 end]), mass(1) * [1 1], 'Color', 0.5 * [1 1 1], 'LineStyle', ':');
line(time([1 end]), mass(1) * [1 1], 'Color', 0.5 * [1 1 1], 'LineStyle', ':');
ylim([min(mass) max(mass)] + 0.05 * range(mass) * [-1 1]);