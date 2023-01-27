classdef experiment < handle
    properties
        model
        %model_simulation
        input = [];
        t_total_min
        t_total_sec
        t_total_frames
        fr_per_sec = 10;
        time
        time_min
        pre_stimulus_min
        post_stimulus_min
        pre_stimulus_frames
        post_stimulus_frames
        minute_to_frames
    end
    
    methods
        function obj = experiment(model)
            if nargin>=1; obj.model = model; end
            %obj.model_simulation = model_simulation;
            %obj.set_up_doses();
        end
        
        function time_min = get.time_min(obj)
            time_min = obj.correct_time(obj.time);
        end
        
        function time_min = correct_time(obj, time)
            time_min = (time+1.5)./60 -obj.pre_stimulus_min; % convert to mins and shift+correct
        end
        
        function time = reverse_correct_time(obj, time_min)
            time = 60.*(time_min + obj.pre_stimulus_min) - 1.5;
        end
        
        function minute_to_frames = get.minute_to_frames(obj)
            minute_to_frames = 60*obj.fr_per_sec;
        end
        
        function pre_stimulus_frames = get.pre_stimulus_frames(obj)
            pre_stimulus_frames = obj.pre_stimulus_min*obj.minute_to_frames;
        end
        
        function post_stimulus_frames = get.post_stimulus_frames(obj)
            post_stimulus_frames = obj.post_stimulus_min*obj.minute_to_frames;
        end
        
        function t_total_sec = get.t_total_sec(obj)
            t_total_sec = obj.t_total_min*60;
        end
        
        function t_total_frames = get.t_total_frames(obj)
            t_total_frames = obj.t_total_min*obj.minute_to_frames;
        end
                
        function set_up_model(obj, model)
            obj.model = model;
        end
        
        function set_up_time(obj)
            obj.time = 0:1/obj.fr_per_sec:(obj.t_total_sec-1/obj.fr_per_sec);
        end
        
        function input = get_input(obj, t)
            input = 0;
            if ~isempty(obj.input)
                index = round(t*obj.fr_per_sec+1);
                %[~, index] = min(abs(obj.time-t));
                input = obj.input(index, :);
            end
        end
        
        function set_up_input(obj)
        end
    end
end