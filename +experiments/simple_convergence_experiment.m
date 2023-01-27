classdef simple_convergence_experiment < experiments.experiment
    %SIMPLE_CONVERGENCE_EXPERIMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
        
    methods
        function obj = simple_convergence_experiment()
            obj@experiments.experiment();
        end
        
        function input = get_input(obj, ~)
            input = obj.input;
        end
        
        function input = set_up_input(obj, input)
            if nargin<2; input = 0; end
            obj.input = input;
        end
    end
end

