%State
% State variable base class. Defines the interface states expose to the
% simulation and processes.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
classdef State < handle
    properties (Abstract = true)
        meta %metadata: id, name, list of dynamic properties, units
    end
    
    properties
        randStream %handle to simulation random stream
    end
    
    %partitioning
    properties
        %used by parent
        partitions = cell(0, 1) %cell array of interacting processes and handles to partition calculation functions
        
        %used by children
        parentState
        parentProcess
    end
    
    %initialization
    methods
        %constructor
        function this = State(propValues)
            if nargin == 0
                return;
            end
            
            fields = fieldnames(propValues);
            for iField = 1:numel(fields)
                this.(fields{iField}) = propValues.(fields{iField});
            end
        end
        
        %construct state-process graph, calculate constants
        function this = initialize(this, sim, kb) %#ok<INUSD>
            this.randStream = sim.randStream;        
        end
        
        %allocate memory
        function this = allocate(this)
            for iPartition = 1:numel(this.partitions)
                this.partitions{iPartition}.allocate();
            end
        end
        
        %calculate initial conditions
        function this = calcInitialConditions(this)
        end
    end
    
    %partitioning into substates
    methods
        function partition = addPartition(this, process, varargin)
            partition = this.constructPartition(process, varargin{:});
            this.partitions = [this.partitions; {partition}];
        end
        
        function partition = constructPartition(this, process)
            %get property values of this instance
            metaClass = metaclass(this);
            propVals = struct();
            for iProp = 1:numel(metaClass.PropertyList)
                prop = metaClass.PropertyList(iProp);
                if ~(prop.Abstract || prop.Constant || prop.Dependent)
                    propVals.(prop.Name) = this.(prop.Name);
                end
            end
            
            propVals.meta.id = sprintf('%s_%s', propVals.meta.id, process.meta.id);
            propVals.meta.name = sprintf('%s - %s', propVals.meta.name, process.meta.name);
            
            propVals.partitions = cell(0, 1);
            
            propVals.parentState = this;
            propVals.parentProcess = process;
            
            %create partition
            partition = feval(class(this), propVals);
        end
        
        %prepare to partion state
        function this = prepartition(this)
        end
        
        %partition state among processes
        function this = partition(this)
            for iPartition = 1:numel(this.partitions)
                p = this.partitions{iPartition};
                for iDynamic = 1:numel(this.meta.dynamics)
                    d = this.meta.dynamics{iDynamic};
                    p.(d) = this.(d);
                end
            end
        end
        
        %merge sub-states partitioned to processes
        %- Copy value from substate with write access
        function this = merge(this)
            for iDynamic = 1:numel(this.meta.dynamics)
                d = this.meta.dynamics{iDynamic};
                oldVal = this.(d);
                newVal = oldVal;
                nNewVal = 0;
                
                for iPartition = 1:numel(this.partitions)
                    p = this.partitions{iPartition};
                    
                    if ~isequal(oldVal, p.(d))
                        newVal = p.(d);
                        nNewVal = nNewVal + 1;
                    end
                end
                
                if nNewVal > 1
                    throw(MException('State:multiProcessWriting', 'Multiple processes cannot simultaneously edit state properties'))
                end
                this.(d) = newVal;
            end
        end
    end
    
    %calculations
    methods
        %calculate (and cache) any dependent properties
        function this = calculate(this)
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
        
        function val = getDynamics(this)
            val = struct();
            
            if ~isfield(this.meta, 'dynamics')
                return;
            end
            
            for iProp = 1:numel(this.meta.dynamics)
                p = this.meta.dynamics{iProp};
                val.(p) = this.(p);
            end
        end
        
        function this = setDynamics(this, val)
            fields = fieldnames(val);
            
            if isempty(fields)
                return;
            end
            
            if ~isfield(this.meta, 'dynamics') || ~all(ismember(fields, this.meta.dynamics))
                throw(MException('State:invalidDynamics', 'Invalid dynamics'))
            end
            
            for iProp = 1:numel(fields)
                p = fields{iProp};
                this.(p) = val.(p);
            end
        end
    end
end