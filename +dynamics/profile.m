classdef profile < handle
    properties
        vars_cont_set = {}
        LP_inxs_set = {}
        vars_cont = []
        LP_inxs = []
        state_stable = []
        ret_status = 0;
    end

    properties (SetAccess = private)
        parent_cont
    end

    methods (Access = public)
        function vec = get_stable_branches(obj, stable)
            vec = nan;
            if isempty(obj.state_stable); return; end
            inxs = find(obj.state_stable == stable);
            ind_jumps = find(diff(inxs)>1)+1;
            vec = nan(size(obj.vars_cont,1),numel(inxs)+numel(ind_jumps));
            if isempty(inxs); return; end
            ind_jumps = [1,ind_jumps,numel(inxs)+1];
            for i = 1:(numel(ind_jumps)-1)
                vec(:,(ind_jumps(i)+(i-1)):(ind_jumps(i+1)-1+(i-1))) = obj.vars_cont(:,inxs(ind_jumps(i)):inxs(ind_jumps(i+1)-1));
            end
        end
    end

    methods (Access = ?dynamics.continuation)  %Only dynamics.continuation is allowed to construct a child
        function obj = profile(parent_cont)
            obj.parent_cont = parent_cont;            
        end

        function append(obj, vars_cont_single, LP_inxs_single, ret_status_single)
            obj.vars_cont_set{1, end+1} = vars_cont_single;
            obj.LP_inxs_set{1, end+1} = LP_inxs_single;
            obj.ret_status = obj.ret_status | ret_status_single;
        end

        function calc_stability(obj, bif_par)
            obj.state_stable = zeros(1,size(obj.vars_cont,2));
            %eigvals = zeros(size(obj.vars_cont,1)-1,size(obj.vars_cont,2));
            for i = 1:size(obj.vars_cont,2)
                state = obj.vars_cont(:,i);
                if any(isnan(state)); continue; end
                [~, stable, ~] = obj.parent_cont.model.jacobian(state, bif_par);
                if ~stable; continue; end
                obj.state_stable(i) = stable;
            end
            % let stable branches take the LPts (if there is a stable branch right next to LP)
            obj.state_stable(obj.LP_inxs) = any(obj.state_stable([obj.LP_inxs'-1, obj.LP_inxs', obj.LP_inxs'+1])',1);
            obj.state_stable(isnan(obj.vars_cont(1,:))) = nan;
            %eigvals(:,isnan(obj.vars_cont(1,:))) = nan;
        end

        function calc_stability_v1(obj)
            % determine stability for each branch of the bifurcation profile, using the jacobian eigenvalues
            % find limit and breaking (nan) points first
            pts_lb = sort([find(isnan(obj.vars_cont(1,:))),obj.LP_inxs]);
            pts_rep = [1, 1+pts_lb]; % select first points of the branches as representative
            pts_end = [1, pts_lb, size(obj.vars_cont,2)];
            obj.state_stable = zeros(1,size(obj.vars_cont,2));
            for i = 1:numel(pts_rep)
                state_rep = obj.vars_cont(:,pts_rep(i)); % use bifpar+state, cause of the fp_cont that is used
                [~, stable] = obj.parent_cont.model.jacobian(state_rep);
                if ~stable; continue; end
                % let stable branches take the LPts
                obj.state_stable(pts_end(i):pts_end(i+1)) = stable;
            end
            obj.state_stable(isnan(obj.vars_cont(1,:))) = nan;
        end

        function sort_link_profiles(obj)
            end_points = nan(2*size(obj.vars_cont_set, 2), 2);
            graph = zeros(size(end_points, 1), size(end_points, 1));
            % get data for end-points, and link them into a graph
            for i = 1:size(obj.vars_cont_set, 2)
                end_points(2*i-1, :) = [obj.vars_cont_set{i}(1, 1), obj.parent_cont.readout(obj.vars_cont_set{i}(:, 1))];
                end_points(2*i, :) = [obj.vars_cont_set{i}(1, end), obj.parent_cont.readout(obj.vars_cont_set{i}(:, end))];
                graph(2*i-1, 2*i) = 1;
                graph(2*i, 2*i-1) = -1;
            end
            s_min = min(end_points(:, 1));
            s_max = max(end_points(:, 1));
            s_min_inxs = find(abs(end_points(:, 1) -s_min)<1e-5);
            s_max_inxs = find(abs(end_points(:, 1) -s_max)<1e-5);
            [~, inxs_1_order] = sort(end_points(s_min_inxs, 2), obj.parent_cont.order_mode);
            %[~, inxs_2_order] = sort(end_points(s_min_inxs, 2), obj.parent_cont.order_mode);
            
            for i_o = 1:numel(inxs_1_order)
                % denote starting and ending nodes in the graph
                ind_start = s_min_inxs(inxs_1_order(i_o));
                %ind_end = s_max_inxs(ind2);
                
                q = {[ind_start]}; %#ok<NBRAK>
                solutions = [];
                
                while ~isempty(q)
                    curr_node_list = q{1};
                    q = q(1, 2:end);
                    if isempty(curr_node_list)
                        break;
                    end
                    if (length(curr_node_list) == size(graph, 1))
                        % check if viable - paths do not cross
                        no_sol = false;
                        pos_pairs = end_points([curr_node_list(1:2:end),curr_node_list(end)],1);
                        virt_pairs = reshape([-inf,end_points(curr_node_list,2)',inf],2,[])';
                        for i = 2:size(virt_pairs,1)
                            for j = 1:(i-1)
                                if pos_pairs(i)~=pos_pairs(j); continue; end
                                if mod(sum(sum(virt_pairs(i,:)'>virt_pairs(j,:))),2)==1
                                    % if the paths are crossing (virt_pairs positions are intercrossing) - not viable
                                    no_sol = true;
                                    break;
                                end
                            end
                            if no_sol; break; end
                        end
                        if ~no_sol
                            solutions = [solutions; curr_node_list]; %#ok<AGROW>
                        end
                        continue;
                    end
                    next_id = setdiff(find(graph(curr_node_list(end), :)), curr_node_list);
                    if isempty(next_id)
                        if ismember(curr_node_list(end), s_min_inxs)
                            next_ids = setdiff(s_min_inxs(sum(abs(graph(s_min_inxs, :)),2)<2), curr_node_list);
                        else
                            next_ids = setdiff(s_max_inxs(sum(abs(graph(s_max_inxs, :)),2)<2), curr_node_list);
                        end
                        if numel(next_ids)>1
                            [~,ind_sort] = sort(pdist2(end_points(curr_node_list(end),:),end_points(next_ids,:)));
                            next_ids = next_ids(ind_sort);
                        end
                        for i = 1:length(next_ids)
                            next_node_list = [curr_node_list, next_ids(i)];
                            q{1, end+1} = next_node_list; %#ok<AGROW>
                        end
                    else
                        next_node_list = [curr_node_list, next_id];
                        q{1, end+1} = next_node_list; %#ok<AGROW>
                    end
                end
                if ~isempty(solutions)
                    break;
                end
            end
            if isempty(solutions)
                disp('problem linking');
                return;
            end
            order = ceil(solutions(1, 1:2:end)/2); % ceil(indices/2) gives the vars_cont index
            flip = mod(solutions(1, 1:2:end),2)==0; % odd positions should have odd indices, otherwise flip direction within vars_cont
            for i = 1:length(flip)
                if flip(i)
                    obj.vars_cont_set{order(i)} = fliplr(obj.vars_cont_set{order(i)});
                    obj.LP_inxs_set{order(i)} = size(obj.vars_cont_set{order(i)}, 2) - obj.LP_inxs_set{order(i)} + 1;
                end
            end
            obj.vars_cont_set = obj.vars_cont_set(order);
            obj.LP_inxs_set = obj.LP_inxs_set(order);
            for i = 1:length(obj.vars_cont_set)
                if isempty(obj.vars_cont)
                    obj.LP_inxs = obj.LP_inxs_set{i};
                    obj.vars_cont = obj.vars_cont_set{i};
                else
                    obj.LP_inxs = [obj.LP_inxs, size(obj.vars_cont, 2) + 1 + obj.LP_inxs_set{i}];
                    obj.vars_cont = [obj.vars_cont, nan(size(obj.vars_cont,1),1), obj.vars_cont_set{i}];
                end
            end
        end
    end
end