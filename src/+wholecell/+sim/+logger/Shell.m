%Shell
% Prints a very breif summary of a whole-cell simulation to standard output
% (command window).
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Stanford University
% Created: 3/3/2013
classdef Shell < wholecell.sim.logger.Logger
    properties
        iterFreq = 1;
    end
    
    properties (Access = protected)
        nLines
        startTime
    end
    
    properties (Access = protected)
        columns
    end
    
    methods
        function this = Shell()
            this = this@wholecell.sim.logger.Logger();
        end
        
        function this = initialize(this, sim)
            %array of columns
            this.columns = [
                struct('header', 'Time (s)', 'state', 'Time', 'property', 'value', 'length', 8, 'format', 'd', 'sum', false)
                struct('header', 'Mass (fg)', 'state', 'Mass', 'property', 'cell', 'length', 9, 'format', '.2f', 'sum', true)
                struct('header', 'Growth (fg/h)', 'state', 'Metabolism', 'property', 'growth', 'length', 13, 'format', '.2f', 'sum', false)
                ];
            
            %collect metadata
            this.nLines = -1;
            this.startTime = now;
            
            %print headers
            for iColumn = 1:numel(this.columns)
                this.columns(iColumn).stateIdx = sim.getStateIndex(this.columns(iColumn).state);
                if iColumn > 1
                    fprintf('  ');
                end
                fprintf(['%' num2str(this.columns(iColumn).length)  's'], this.columns(iColumn).header);
            end
            fprintf('\n');
            
            for iColumn = 1:numel(this.columns)
                if iColumn > 1
                    fprintf('  ');
                end
                fprintf(['%' num2str(this.columns(iColumn).length)  's'], repmat('=', 1, this.columns(iColumn).length));
            end
            fprintf('\n');
            
            %print initial state
            this.append(sim);
        end
        
        function this = append(this, sim)
            this.nLines = this.nLines + 1;
            if mod(this.nLines, this.iterFreq) ~= 0
                return
            end
            
            for iColumn = 1:numel(this.columns)
                c = this.columns(iColumn);
                if iColumn > 1
                    fprintf('  ');
                end
                val = sim.states{c.stateIdx}.(c.property);
                if c.sum
                    val = sum(val);
                end
                fprintf(['%' num2str(c.length)  c.format], val);
            end
            fprintf('\n');
        end
        
        function this = finalize(this, sim)
            %print summary
            fprintf('\n');
            fprintf('Simulation finished:\n');
            
            %length
            h = floor(sim.getState('Time').value / 3600);
            m = floor(mod(sim.getState('Time').value, 3600) / 60);
            s = mod(sim.getState('Time').value, 60);
            fprintf('- Length: %d:%02d:%02.0f\n', h, m, s);
            
            %runtime
            [~, ~, d, h, m, s] = datevec(now - this.startTime);
            fprintf('- Runtime: %d:%02d:%02.0f\n', 24 * d + h, m, s);
            
            %new line
            fprintf('\n');
        end
    end
end