%Disk
% Logs whole-cell simulations and metadata to disk. Also provides a
% function (load) for reading stored simulation data.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
classdef Disk < wholecell.sim.logger.Logger
    properties (SetAccess = protected)
        metadata   %metadata
        outDir     %output directory
        segmentLen %number of time steps per segment
    end
    
    properties (Access = protected)
        iSegment
        iStep
        stateLog
        randStreamLog
    end
    
    methods
        function this = Disk(varargin)
            import wholecell.util.mkdir;
            
            %super class constructor
            this = this@wholecell.sim.logger.Logger();
            
            %process input arguments
            ip = inputParser;
            ip.addParamValue('metadata', struct(), @(x) isstruct(x) && this.isMetadaValid(x));
            ip.addParamValue('outDir', '', @(x) ischar(x) && ~isempty(x));
            ip.addParamValue('segmentLen', 1000, @(x) isnumeric(x) && mod(x, 1) == 0);
            ip.parse(varargin{:});
            
            %set properties
            fields = fieldnames(ip.Results);
            for iField = 1:numel(fields);
                this.(fields{iField}) = ip.Results.(fields{iField});
            end
            
            %initialization
            mkdir(this.outDir, true);
            delete(fullfile(this.outDir, '*.mat'));
        end
        
        function this = initialize(this, sim)
            %% metadata
            this.metadata.startTime = datestr(clock, 'yyyy-mm-dd HH:MM:SS');
            this.metadata.endTime = [];
            this.metadata.lengthSec = [];
            this.metadata.timeStepSec = [];
            this.metadata.segmentLen = this.segmentLen;
            
            %% save initial state
            %setup segment, step counters
            this.iSegment = 0;
            this.iStep = 1;
            
            %allocate memory
            this.allocateMemory(sim, 1);
            
            %copy data from states
            this.copyDataFromStates(sim);
            
            %save initial state to disk
            this.saveSegmentToDisk();
            
            %% setup to save dynamics
            %setup segment, step counters
            this.iSegment = 1;
            this.iStep = 0;
            
            %allocate memory
            this.allocateMemory(sim, this.segmentLen);
        end
        
        function this = append(this, sim)
            %increment step counter
            this.iStep = this.iStep + 1;
            
            %copy data from states
            this.copyDataFromStates(sim);
            
            %rotate segment
            if this.iStep == this.segmentLen
                %save segment to disk
                this.saveSegmentToDisk();
                
                %increment segment counter
                this.iStep = 0;
                this.iSegment = this.iSegment + 1;
            end
        end
        
        function this = finalize(this, sim)
            %% save final segment
            if this.iStep > 0
                %contract segment
                this.contractSegment(sim);
                
                %save last segment to disk
                this.saveSegmentToDisk();
            end
            
            %% metadata
            %record
            this.metadata.endTime = datestr(clock, 'yyyy-mm-dd HH:MM:SS');
            this.metadata.lengthSec = sim.getState('Time').value;
            this.metadata.timeStepSec = sim.timeStepSec;
            
            %save
            this.saveMetadata(sim.getOptions(), sim.getParameters());            
            
            %% clear log
            this.iSegment = [];
            this.iStep = [];
            this.stateLog = [];
            this.randStreamLog = [];
        end
    end
    
    methods (Access = protected)
        %TODO: implement. check for
        %- name, description
        %- investigator
        %- revision
        %- username, hostname, ip address
        function val = isMetadaValid(this, metadata) %#ok<INUSD>
            val = true;
        end
        
        function allocateMemory(this, sim, nSteps)
            %state
            this.stateLog = struct();
            
            for iState = 1:numel(sim.states)
                s = sim.states{iState};
                this.stateLog.(s.meta.id) = struct();
                for iProp = 1:numel(s.meta.dynamics)
                    p = s.meta.dynamics{iProp};
                    this.stateLog.(s.meta.id).(p) = zeros([size(s.(p)) nSteps]);
                end
            end
            
            %rand stream
            this.randStreamLog = zeros(numel(sim.randStream.state), nSteps);
        end
        
        function contractSegment(this, sim)
            %state
            for iState = 1:numel(sim.states)
                s = sim.states{iState};
                for iProp = 1:numel(s.meta.dynamics)
                    p = s.meta.dynamics{iProp};
                    this.stateLog.(s.meta.id).(p) = this.stateLog.(s.meta.id).(p)(:, :, 1:this.iStep);
                end
            end
            
            %rand stream
            this.randStreamLog = this.randStreamLog(:, 1:this.iStep);
        end
        
        function copyDataFromStates(this, sim)
            %state
            for iState = 1:numel(sim.states)
                s = sim.states{iState};
                for iProp = 1:numel(s.meta.dynamics)
                    p = s.meta.dynamics{iProp};
                    this.stateLog.(s.meta.id).(p)(:, :, this.iStep) = s.(p);
                end
            end
            
            %rand stream
            this.randStreamLog(:, this.iStep) = sim.randStream.state;
        end
        
        function saveSegmentToDisk(this)
            %dynamics
            tmp = this.stateLog; %#ok<NASGU>
            save(fullfile(this.outDir, sprintf('state-%d.mat', this.iSegment)), '-struct', 'tmp');
            
            %rand stream
            tmp = struct('randStream', this.randStreamLog); %#ok<NASGU>
            save(fullfile(this.outDir, sprintf('randStream-%d.mat', this.iSegment)), '-struct', 'tmp');
        end
        
        %save metadata, options, parameters
        function saveMetadata(this, options, parameters) %#ok<INUSD>
            metadata = this.metadata; %#ok<PROP,NASGU>
            save(fullfile(this.outDir, 'metadata.mat'),   '-struct', 'metadata');
            save(fullfile(this.outDir, 'options.mat'),    '-struct', 'options');
            save(fullfile(this.outDir, 'parameters.mat'), '-struct', 'parameters');
        end
    end
    
    methods (Static = true)
        function value = load(dir, state, prop)
            %load metadata
            md = load(fullfile(dir, 'metadata.mat'));
            
            %load first time point
            tmp = load(fullfile(dir, 'state-0.mat'), state);
            
            %allocate memory for result
            value = zeros([size(tmp.(state).(prop))  md.lengthSec / md.timeStepSec + 1]);
            
            %record first segment
            value(:, :, 1) = tmp.(state).(prop);
            
            %load subsequent segments
            for iTime = md.timeStepSec:md.segmentLen*md.timeStepSec:md.lengthSec
                iSegment = (iTime - md.timeStepSec) / (md.timeStepSec * md.segmentLen) + 1; %#ok<PROP>
                tmp = load(fullfile(dir, sprintf('state-%d.mat', iSegment)), state); %#ok<PROP>
                
                minIdx = (iTime - md.timeStepSec) / md.timeStepSec + 1 + 1;
                maxIdx = (iTime - md.timeStepSec) / md.timeStepSec + 1 + md.segmentLen;
                
                value(:, :, minIdx:min(maxIdx, end)) = tmp.(state).(prop);
            end
        end
    end
end