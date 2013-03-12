%Process
% Process submodel base class. Defines interface that processes expose to
% the simulation and to the states.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
classdef Process < handle
    properties (Abstract = true)
        meta
    end
    
    %constants
    properties
        timeStepSec
    end
    
    %handles
    properties
        randStream %handle to simulation random stream
    end
    
    methods
        %constructor
        function this = Process()
        end
        
        %construct object graph, calculate constants
        function this = initialize(this, sim, kb) %#ok<INUSD>
            this.timeStepSec = sim.timeStepSec;
            this.randStream = sim.randStream;
        end
        
        %calculate submodel contribution to temporal evolution of cell
        %state
        function evolveState(this) %#ok<MANU>
        end
    end
    
    %get, set options, parameters
    methods
        function val = getOptions(this)
            val = struct();
            
            if ~isfield(this.meta, 'options')
                return;
            end
            
            for iOpt = 1:numel(this.meta.options)
                opt = this.meta.options{iOpt};
                val.(opt) = this.(opt);
            end
        end
        
        function this = setOptions(this, val)
            fields = fieldnames(val);
            
            if isempty(fields)
                return;
            end
            
            if ~isfield(this.meta, 'options') || ~all(ismember(fields, this.meta.options))
                throw(MException('State:invalidOptions', 'Invalid option'))
            end
            
            for iOpt = 1:numel(fields)
                opt = fields{iOpt};
                this.(opt) = val.(opt);
            end
        end
        
        function val = getParameters(this)
            val = struct();
            
            if ~isfield(this.meta, 'parameters')
                return;
            end
            
            for iParam = 1:numel(this.meta.parameters)
                param = this.meta.parameters{iParam};
                val.(param) = this.(param);
            end
        end
        
        function this = setParameters(this, val)
            fields = fieldnames(val);
            
            if isempty(fields)
                return;
            end
            
            if ~isfield(this.meta, 'parameters') || ~all(ismember(fields, this.meta.parameters))
                throw(MException('State:invalidParameters', 'Invalid parameter'))
            end
            
            for iParam = 1:numel(fields)
                param = fields{iOpt};
                this.(param) = val.(param);
            end
        end
    end
end