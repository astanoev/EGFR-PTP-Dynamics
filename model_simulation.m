classdef model_simulation < handle
    properties
        input_vec;
        state;
        time;
        experiment;
        model;
        sol;
    end
    
    methods
        function obj = model_simulation(model, experiment)
            if nargin>=1; obj.set_model(model); end
            if nargin>=2; obj.set_experiment(experiment); end
        end
        
        function set_experiment(obj, experiment)
            obj.experiment = experiment;
            obj.experiment.set_up_model(obj.model);
            %obj.input_vec = obj.experiment.set_up_input(); 
        end
        
        function set_model(obj, model)
            obj.model = model;
        end
        
        function simulate(obj)
            [obj.time, obj.state] = obj.time_profile(obj.model.init_conds);
        end
        
        function [T, Y] = time_profile(obj, init_cond)
            %options = odeset('RelTol',1e-8,'AbsTol',1e-10);%,'Events',@obj.eventfun);
            options = odeset('RelTol',1e-3,'AbsTol',1e-6);
            obj.sol = ode45(@(tt,y) obj.model.df_model(tt, y, obj.experiment), [obj.experiment.time(1) obj.experiment.time(end)], init_cond, options);
            T = obj.sol.x;
            Y = obj.sol.y;
        end
    end
end