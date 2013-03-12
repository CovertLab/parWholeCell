%Simulation
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
classdef Simulation < handle
    %metadata
    properties
        meta = wholecell.util.struct(...
            'options', {'lengthSec'; 'timeStepSec'; 'seed'}, ...
            'units', struct(...
                'lengthSec', 's', ...
                'timeStepSec' , 's' ...
                ) ...
            );
    end
    
    %options
    properties
        lengthSec   = 50000; %simulation length (s)
        timeStepSec = 1;     %simulation time step (s)        
    end
    
    properties (Dependent = true)
        seed
    end
    
    %rand stream, state, process handles
    properties
        randStream %rand stream
        
        states     %cell array of states
        time       %time state
        
        processes  %cell array of processes
    end
    
    %construct simulation
    methods
        %constructor
        function this = Simulation(kb)
            this.constructRandStream();
            this.constructStates();
            this.constructProcesses();
            this.initialize(kb);
            this.allocateMemory();
        end
    end
    
    methods (Access = protected)
        %construct rand stream
        function this = constructRandStream(this)
            import wholecell.util.RandStream;
            
            this.randStream = RandStream('mcg16807');
        end
        
        %construct states
        function this = constructStates(this)
            this.states = {
                wholecell.sim.state.Mass()
                wholecell.sim.state.Metabolism()
                wholecell.sim.state.MoleculeCounts()
                wholecell.sim.state.Time()
                };
            
            this.time = this.getState('Time');
        end
        
        %construct processes
        function this = constructProcesses(this)
            this.processes = {
                wholecell.sim.process.Complexation()
                wholecell.sim.process.Metabolism()
                wholecell.sim.process.ProteinMaturation()
                wholecell.sim.process.RnaDegradation()
                wholecell.sim.process.RnaMaturation()
                wholecell.sim.process.Transcription()
                wholecell.sim.process.Translation()
                };
        end
        
        %link states and processes
        function this = initialize(this, kb)
            for iState = 1:numel(this.states)
                this.states{iState}.initialize(this, kb);
            end
            
            for iProcess = 1:numel(this.processes)
                this.processes{iProcess}.initialize(this, kb);
            end
        end        
        
        %allocate memory
        function this = allocateMemory(this)
            for iState = 1:numel(this.states)
                this.states{iState}.allocate();
            end
        end
    end
    
    %run simulation
    methods
        %run simulation
        function this = run(this, loggers)
            %process arguments
            if nargin < 2
                loggers = cell(0, 1);
            elseif ~iscell(loggers)
                loggers = {loggers};
            end
            
            %calculate initial conditions
            this.calcInitialConditions();
            
            %initialize logs
            for iLogger = 1:numel(loggers)
                loggers{iLogger}.initialize(this);
            end
            
            %calculate temporal evolution
            for iSec = this.timeStepSec:this.timeStepSec:this.lengthSec
                this.time.value = iSec;
                this.evolveState();
                
                %append logs
                for iLogger = 1:numel(loggers)
                    loggers{iLogger}.append(this);
                end
            end
            
            %finalize logs
            for iLogger = 1:numel(loggers)
                loggers{iLogger}.finalize(this);
            end
        end
        
        %calculate initialial conditions
        function this = calcInitialConditions(this)
            %calculate initial conditions
            this.getState('Time').calcInitialConditions();
            this.getState('MoleculeCounts').calcInitialConditions();
            this.getState('Mass').calculate();
            this.getState('Mass').partition();
            this.getState('Metabolism').calcInitialConditions();
            
            %calculate dependent state
            for iState = 1:numel(this.states)
                s = this.states{iState};
                s.calculate();
            end
        end
        
        %calculate temporal evolution
        function this = evolveState(this)
            %prepare to partition states among processes
            for iState = 1:numel(this.states)
                s = this.states{iState};
                s.prepartition();
            end
            
            %partition states among processes
            for iState = 1:numel(this.states)
                s = this.states{iState};
                s.partition();
            end
            
            %simulate submodels
            for iProcess = 1:numel(this.processes)
                p = this.processes{iProcess};
                p.evolveState();
            end
            
            %partition state among processes
            for iState = 1:numel(this.states)
                s = this.states{iState};
                s.merge();
            end
            
            %recalculate dependent state
            for iState = 1:numel(this.states)
                s = this.states{iState};
                s.calculate();
            end
        end
    end
    
    %helper methods
    methods
        function value = getStateIndex(this, stateId)
            value = [];
            for iState = 1:numel(this.states)
                if strcmp(this.states{iState}.meta.id, stateId)
                    value = iState;
                    return;
                end
            end
        end
        
        function value = getProcessIndex(this, processId)
            value = [];
            for iProcess = 1:numel(this.processes)
                if strcmp(this.processes{iProcess}.meta.id, processId)
                    value = iProcess;
                    return;
                end
            end
        end
        
        function value = getState(this, stateId)
            iState = this.getStateIndex(stateId);
            if ~isempty(iState)
                value = this.states{iState};
            else
                value = [];
            end
        end
        
        function value = getProcess(this, processId)
            iProcess = this.getProcessIndex(processId);
            if ~isempty(iProcess)
                value = this.processes{iProcess};
            else
                value = [];
            end
        end
    end
    
    %get, set options, parameters
    methods
        %returns options as struct
        function val = getOptions(this)
            %initialize output
            val = struct('states', struct(), 'processes', struct());
            
            %top-level
            if isfield(this.meta, 'options')
                for iOpt = 1:numel(this.meta.options)
                    opt = this.meta.options{iOpt};
                    val.(opt) = this.(opt);
                end
            end
            
            %states
            for iState = 1:numel(this.states)
                s = this.states{iState};
                val.states.(s.meta.id) = s.getOptions();
            end
            
            %processes
            for iProcess = 1:numel(this.processes)
                p = this.processes{iProcess};
                val.processes.(p.meta.id) = p.getOptions();
            end
        end
        
        %sets options values based on passed in struct
        function this = setOptions(this, val)
            %top-level
            opts = setdiff(fieldnames(val), {'states'; 'processes'});
            if ~isempty(opts) && (~isfield(this.meta, 'options') || ~all(ismember(opts, this.meta.options)))
                tfs = ismember(opts, this.meta.options);
                throw(MException('Simulation:invalidOption', 'Invalid options:\n- %s', ...
                    strjoin(sprintf('\n- '), opts{~tfs})));
            end
            for iOpt = 1:numel(opts)
                opt = opts{iOpt};
                this.(opt) = val.(opt);
            end
            
            %states
            if isfield(val, 'states')
                fields = fieldnames(val);
                for iField = 1:numel(fields)
                    s = this.getState(fields{iField});
                    s.setOptions(val.(fields{iField}));
                end
            end
            
            %processes
            if isfield(val, 'processes')
                fields = fieldnames(val);
                for iField = 1:numel(fields)
                    p = this.getProcess(fields{iField});
                    p.setOptions(val.(fields{iField}));
                end
            end
        end
        
        %returns parameters as struct
        function val = getParameters(this)
            %initialize output
            val = struct('states', struct(), 'processes', struct());
            
            %states
            for iState = 1:numel(this.states)
                s = this.states{iState};
                val.states.(s.meta.id) = s.getParameters();
            end
            
            %processes
            for iProcess = 1:numel(this.processes)
                p = this.processes{iProcess};
                val.processes.(p.meta.id) = p.getParameters();
            end
        end
        
        %sets parameters values based on passed in struct
        function this = setParameters(this, val)
            %states
            if isfield(val, 'states')
                fields = fieldnames(val);
                for iField = 1:numel(fields)
                    s = this.getState(fields{iField});
                    s.setParameters(val.(fields{iField}));
                end
            end
            
            %processes
            if isfield(val, 'processes')
                fields = fieldnames(val);
                for iField = 1:numel(fields)
                    p = this.getProcess(fields{iField});
                    p.setParameters(val.(fields{iField}));
                end
            end
        end
        
        function val = getDynamics(this)
            val = struct();
            for iState = 1:numel(this.states)
                s = this.states{iState};
                val.(s.meta.id) = s.getDynamics();
            end
        end
        
        function this = setDynamics(this, val)
            fields = fieldnames(val);
            for iField = 1:numel(fields)
                s = this.getState(fields{iField});
                s.setDynamics(val.(fields{iField}));
            end
        end
        
        function set.timeStepSec(this, val)
            this.timeStepSec = val;
            for iProcess = 1:numel(this.processes) %#ok<MCSUP>
                this.processes{iProcess}.timeStepSec = val; %#ok<MCSUP>
            end
        end
        
        function val = get.seed(this)
            val = this.randStream.seed;
        end
        
        function set.seed(this, val)
            this.randStream.seed = val;
        end
    end
end