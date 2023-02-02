classdef continuation
    %CONTINUATION class for execuing numerical continuation procedure
    % This class is used with an input model, where the continuation derivatives are defined, 
    % in terms of the bifurcation parameter and variables of the system
    % It works with arbitrary number of variables

    % For the EGFR-PTPRG model, the output is defined via readout_wts, as a linear combination of the variables (EGFRp+EGF_EGFRp)
    % The input is passed through the EGF_EGFRt parameter
    % Bifurcation parameter G1 is also used, as all the rest can also be used
    
    properties
        model;
        readout_wts;
        tol_n = 1.0E-09;
        marker = '.';
        order_mode = 'asc'; % or 'desc'
    end

    properties (Transient)
        out
    end
    
    methods
        function obj = continuation(model, readout_wts)
            if nargin < 1
                % if called directly from here, plot DR continuation of the best parameter set
                obj.model = models.egfr_ptprg_model;
                try
                    load(fullfile('data','model_v10_fit_v12.5.mat'),'results');
                    par_new = results.converged.par_bestfit;
                    fn = fieldnames(obj.model.par);
                    subset_cond = 1;
                    for i=1:length(fn)
                        obj.model.par.(fn{i}) = par_new(i +(subset_cond-1)*length(fn));
                    end
                    % parameter values can also be altered before bifurcation
                    % obj.model.par.EGF_EGFRtt = 0.01;
                catch
                    disp('problem');
                end
            else
                obj.model = model;
            end
            if nargin < 2
                obj.readout_wts = zeros(numel(obj.model.labels),1);
                obj.readout_wts([1,3]) = 1;
            else
                obj.readout_wts = readout_wts;
            end
            warning('off','MATLAB:rankDeficientMatrix');
            warning('off','MATLAB:singularMatrix');
            if nargin < 1
                bif_par = 'g1'; % or 'EGF_EGFRtt';
                bif_max_val = 3*obj.model.par.(bif_par); 
                % bif_max_val = 1;
                obj.out = obj.calc_profile(bif_par, [0, bif_max_val], [0, 1], true);
            end
        end

        function sss = find_steady_states(obj, bif_par, x_minmax, y_range, sss_master, varargin)
            p = inputParser;
            addParameter(p,'order_mode',obj.order_mode);
            addParameter(p,'k_steps',2);
            parse(p,varargin{:});
            k = p.Results.k_steps; % how many subsections per cube side - to explore unstable solutions
            if strcmp(bif_par,'EGF_EGFRtt')
                init_conds = [x_minmax(1),0,1,x_minmax(1); x_minmax(1),1-x_minmax(1),0,x_minmax(1);...
                                x_minmax(2),0,1,x_minmax(2); x_minmax(2),1-x_minmax(2),0,x_minmax(2)];
            else
                init_conds = [x_minmax(1),0,1,0; x_minmax(1),1-obj.model.par.EGF_EGFRtt,0,obj.model.par.EGF_EGFRtt;...
                                x_minmax(2),0,1,0; x_minmax(2),1-obj.model.par.EGF_EGFRtt,0,obj.model.par.EGF_EGFRtt];
            end
            m = size(init_conds,2)-1; % num variables
            n = k^m; % use k-base system to fetch all combinations
            sss = nan(2*n, numel(obj.model.labels)+1);
            for i_i = 1:2
                % loop through x_limits
                iconds = init_conds([2*i_i-1,2*i_i],:);
                for i=1:n
                    % loop through combinations of fractions within the cube
                    % implementation is through looping through numbers in base k format
                    wgts = str2num(char(num2cell(dec2base(i-1,k,m))))./(k-1); %#ok<ST2NM> 
                    wgts = [[.5,.5];[wgts, 1-wgts]];
                    ss0 = sum(wgts'.*iconds,1)';
                    [status, ss0] = obj.newton(ss0, bif_par, 1);
                    if status == 0; sss(i+n*(i_i-1),:) = ss0; end
                end
            end
            sss = sss(~any(isnan(sss),2), :);
            sss = uniquetol(sss, 1e-4, 'ByRows', true, 'DataScale', 1);
            if ~isempty(sss_master)
                sss(ismember(round(sss, 4), round(sss_master, 4), 'rows'), :) = []; % remove already calculated ss
            end
            sss((sss(:,1)<x_minmax(1)-1e-5)|(sss(:,1)>x_minmax(2)+1e-5), :) = [];
            sss((obj.readout(sss')<y_range(1)-1e-5)|(obj.readout(sss')>y_range(end)+1e-5), :) = [];
            % remove steady states that go out of bounds for the variables
            sss(any(sss(:,2:end)<min(init_conds(:,2:end),[],1)-1e-5,2)|any(sss(:,2:end)>max(init_conds(:,2:end),[],1)+1e-5,2), :) = [];

            [~, ind_sort] = sortrows([sss(:,1), obj.readout(sss')'], [1, 2]);
            if strncmpi(p.Results.order_mode, 'asc', 3)                
                sss = sss(ind_sort,:);
            else
                sss = sss(flipud(ind_sort),:);
            end
            sss(abs(sss(:,1)-x_minmax(1))<1e-5,1) = x_minmax(1); % snap to min-max
            sss(abs(sss(:,1)-x_minmax(2))<1e-5,1) = x_minmax(2);
        end
        
        function [vars_cont, LP_inxs, ret_status] = calc_profile(obj, bif_par, x_minmax, y_minmax, animate)
            if nargin < 5; animate = false; end
            vars_cont = {};
            LP_inxs = {};
            ret_status = 0;
            sss_master = zeros(0, numel(obj.model.labels)+1);
            k_steps = 3;
            while k_steps < 10
                % detect steady states for a given k_steps value 
                % (k_steps value per variable - k_steps x k_steps x k_steps variable sets to use as initial conditions for a newton step, 
                % to find a steady state solution)
                sss = obj.find_steady_states(bif_par, x_minmax, y_minmax, sss_master, 'k_steps', k_steps);
                processed = zeros(size(sss,1), 1); % vector od processed sss
                problematic = zeros(size(sss,1), 1); % vector of problematic sss (detected ss, but not able to continue from them)
                ind_ss = 1; % current ss that is being processed
                while ~all(processed)
                    if processed(ind_ss)
                        % if ss(ind_ss) is already processed, find an unprocessed one
                        for j = mod(ind_ss:(size(sss, 1)+ind_ss-2), size(sss,1))+1
                            if ~processed(j)
                                ind_ss = j;
                                break;
                            end
                        end
                    end
                    ss = sss(ind_ss, :);
                    processed(ind_ss) = 1; % mark as processed ss
                    % set of unprocessed sss - potential endpoints of the continuation algorithm
                    sss_finals = sss(~processed, :);
                    if ss(1) == x_minmax(1)
                        % if starting from left x-limit, move forward, otherwise move backwards
                        t0 = [1; 0; 0; 0];
                    else
                        t0 = [-1; 0; 0; 0];
                    end
                    x0_minmax = [x_minmax', repmat([0;1],1,3)];
                    % main continuation execution from a given starting point ss' and initial direction t0
                    [vars_cont_single, LP_inxs_single, ret_status_single, animate] = obj.continuation_single(bif_par, ss', sss_finals, x0_minmax, t0, animate);
                    if isempty(vars_cont_single) || ret_status_single>0
                        % if there is no possible continuation from the given starting ss
                        % (rare case, but possible when there is limit stability)
                        problematic(ind_ss) = 1;
                        continue;
                    end
                    ss_final = vars_cont_single(:, end);
                    exists_as_ss = false;
                    for j = 1:size(sss, 1)
                        if processed(j); continue; end
                        if exists_as_ss
                            % select the next ss for processing
                            ind_ss = j;
                            break;
                        elseif norm(ss_final - sss(j, :)') < 2e-4
                            % if the final state exists as ss in the sss dataset
                            processed(j) = 1;
                            exists_as_ss = true;
                            continue; % continue to select the next candidate
                        end
                    end
                    if ~exists_as_ss
                        % if not detected previously, add as a new ss in the sss set
                        sss = [sss; ss_final']; %#ok<AGROW>
                        processed = [processed; 1]; %#ok<AGROW>
                        problematic = [problematic; 0]; %#ok<AGROW>
                        if ~all(processed)
                            ind_ss = find(~processed, 1);
                        end
                    end
                    if ss(1) == x_minmax(2) && abs(ss_final(1)-x_minmax(2))>1e-5
                        % if started from the end, flip everything
                        vars_cont_single = fliplr(vars_cont_single);
                        LP_inxs_single = fliplr(size(vars_cont_single,2)-LP_inxs_single+1);
                    end
                    vars_cont{1, end+1} = vars_cont_single; %#ok<AGROW>
                    LP_inxs{1, end+1} = LP_inxs_single; %#ok<AGROW>
                    ret_status = ret_status | ret_status_single;
                end
                sss_master = [sss_master; sss(~problematic, :)]; %#ok<AGROW>
                if ((mod(length(find(abs(sss_master(:, 1) -x_minmax(1))<1e-5)), 2) ~= 1) || (mod(length(find(abs(sss_master(:, 1) -x_minmax(2))<1e-5)), 2) ~= 1))
                    % if there are undetected steady states, decrease the step size
                    k_steps = k_steps+1;
                else
                    % sort and link lists
                    try
                        [vars_cont, LP_inxs] = obj.sort_link_profiles(vars_cont, LP_inxs);
                    catch
                        disp('problem linking');
                        %model_problem = obj.model;
                        %save('data\temp\problem_model.mat','model_problem');
                    end
                    return;
                end
            end
            disp('warning: not all border steady-states have been identified');
            %model_problem = obj.model;
            %save('data\temp\problem_model.mat','model_problem');
        end

        function [vars_cont, LP_inxs, ret_status, ax] = continuation_single(obj, bif_par, x0, xfinals, x0_minmax, t0, animate)
            if nargin < 2
                x0 = [0; 2];
            end
            if nargin < 3
                xfinals = [];
            end
            if nargin < 4
                % bounding box
                x0_minmax = [[0, 0]; [5, 5]];
            end
            if nargin < 5
                % starting direction - towards right (if you start at zero/from the left)
                if x0(1) == x0_minmax(1,1)
                    t0 = [1; 0];
                else
                    t0 = [-1; 0];
                end
            end
            if nargin < 6
                % animate is for debugging purposes
                animate = false;
            end
            
            p0 = 1; % parameterize variable p0 initially
            % tolerance has been inferred from the data, seems like if the
            % step size is adapted to yield error smaller than ~1e-4,
            % it is likely that solution is 'nearby', i.e. can be found
            tol = 6*1e-5;%1e-04;
            tol1 = 0.7*tol; % if error lower than tol1, increase step size
            tol2 = 0.8*tol; % if error higher than tol, decrease step size
            h = 0.015; % starting step size
            hmax = 0.02; % maximal step size
            hmin = 0.0005; % minimal step size
            step_max = 30000; % maximal number of steps (can change to inf)
            step_num = 1; 
            vars_cont = zeros(numel(x0), step_max); % save solutions
            LP_inxs = []; % save limit points
            
            % just for safety - if already at a final stop
            for i = 1:size(xfinals, 1)
                if norm(x0-xfinals(i, :)') < obj.tol_n
                    ret_status = 0;
                    vars_cont = zeros(2, 0);
                    return;
                end
            end

            vars_cont(:, step_num) = x0(:, 1);
            
            if ~ (animate == false)
                % if called from within continuation - animate solutions
                if ~ishandle(animate)
                    fig = figure();
                    ax = axes('Parent',fig);
                    hold(ax,'on'); box(ax,'on');
                    if strcmp(bif_par,'EGF_EGFRtt')
                        axis(ax,'equal');
                    end
                    xlim(ax, x0_minmax(:, 1)');
                    ylim(ax, x0_minmax(:, 2)');
                else
                    ax = animate;
                end
                bif_profile = animatedline(ax,'LineWidth',1.5,'Color',[1,0,0],'Marker',obj.marker);
                tang = plot(ax, [nan], [nan],'b-','LineWidth',1); %#ok<NBRAK>
                addpoints(bif_profile, x0(1), obj.readout(x0));
            else
                ax = false; % for returning false, if called from outside continuation
            end
            
            hh = h; % starting step size
            while step_num < step_max
                step_num = step_num + 1;
                if step_num > 2
                    t0 = vars_cont(:, step_num-1)-vars_cont(:, step_num-2);
                    t0 = t0./norm(t0);
                end
                found_solution = false;
                % loop and halve step size every iteration if no solution
                % is found, but continue reducing step size if at first
                % time step
                while (hh >= hmin) || ((step_num==2) && (hh >= 1e-10))
                    % loop while parameterizing different variable
                    j = 1;
                    while j <= numel(x0)
                        [status, x2, t2, p2] = obj.step(x0, t0, bif_par, p0, hh, x0_minmax);
                        if (status ~= 0) % if not solved
                            if ~ (animate == false) && ~any(imag(x2) ~= 0)
                                set(tang,'XData',[x0(1), x2(1)]);
                                set(tang,'YData',[obj.readout(x0), obj.readout(x2)]);
                                drawnow;
                            end
                            p0 = mod(p0,numel(x0))+1;
                            j = j + 1;
                        else % if solution is found
                            if step_num == 2 && obj.out_of_bounds(x2, x0_minmax)
                                % if initial move is out of bounds, reduce step
                                if hh >= hmin
                                    hh = hh/2;
                                else
                                    hh = h;
                                    p0 = mod(p0,numel(x0))+1;
                                    j = j+1;
                                end
                                status = 1;
                                continue;
                            end
                            tp2 = obj.tangent(x2, t2, bif_par, p2);
                            if sign(tp2(1))~=sign(t2(1))
                                % finely search between the two points
                                % for the limit point (saddle node)
                                x3 = obj.binary_search(0, x0, x2, t2, tp2, bif_par, x0_minmax);
                                if isnan(x3)
                                    p0 = mod(p0,numel(x0))+1;
                                    j = j+1;
                                    continue;
                                end
                                vars_cont(:, step_num) = x3(:, 1); % add to solutions
                                LP_inxs = [LP_inxs, step_num]; %#ok<AGROW>
                                step_num = step_num + 1;
                                if ~ (animate == false)
                                    plot(ax, x3(1), obj.readout(x3), 'ko');
                                    text_offset = abs(x0_minmax(2,1)-x0_minmax(1,1))*(-0.055*(x0(1)>x3(1))+0.02*(x3(1)>x0(1)));
                                    text(ax, x3(1)+text_offset, obj.readout(x3), 'LP');
                                end
                            end
                            % adapt step size for next iteration,
                            % relative to current error
                            err = norm(x2 - (x0 + hh * t2));
                            if err < tol1
                                hh = min([hh * tol1/err, hh * 10, hmax]);
                            elseif err > tol
                                if hh ~= hmin
                                    hh = max(hh * tol2/err, hmin);
                                end
                            end
                            found_solution = true;
                            break;
                        end
                    end
                    if found_solution
                        break;
                    end
                    if step_num == 2
                        % if no solution is found from the start - abort
                        ret_status = 1;
                        return;
                    end
                    hh = hh/2; % halve step size
                end
                hh = max(hh, hmin);
                
                % if still no solution is found (very sharp cusp and
                % unreliable tangent direction) do brute-force direction
                % search
                if status ~= 0
                    %return;
                    % disp(strcat('exhaustive search mode at: ',num2str(x0(1)),', ',num2str(x0(2))));
                    n_angles = 35; % 35 angles (10deg)
                    hh = hmin; % very small step size
                    while hh <= hmax
                        x2s = zeros(2*n_angles, 2);
                        t02s = zeros(2*n_angles, 2);
                        statuses = zeros(2*n_angles, 1);
                        dot_prods = zeros(2*n_angles, 1);
                        for p0 = 1:numel(x0)
                            for i=1:n_angles
                                ind_i = (p0-1)*n_angles +i;
                                % exclude the backwards angle (0 from the offset)
                                angle_offset = atan2(-t0(2),-t0(1));
                                t01 = [cos(angle_offset + i/(n_angles+1)*2*pi); sin(angle_offset + i/(n_angles+1)*2*pi)];
                                x1 = x0 + hh*t01;
                                [statuses(ind_i), x2s(ind_i, :)] = obj.newton(x1, bif_par, p0);
                                t02 = x2s(ind_i, :)' - x0;
                                t02 = t02./norm(t02);
                                t02s(ind_i,:) = t02';
                                dot_prods(ind_i) = t0'*t02;
                            end
                        end
                        % take distinct solution that is in the most forward direction
                        inxs = find((statuses == 0) & (vecnorm((x2s-repmat(x0',size(x2s,1),1)),2,2)<h) & (vecnorm((x2s-repmat(x0',size(x2s,1),1)),2,2)>1e-04) & ~any(x2s<0,2));
                        if ~isempty(inxs)
                            [~, ind] = max(dot_prods(inxs));
                            x2 = x2s(inxs(ind), :)';
                            p2 = ceil(inxs(ind)/n_angles);
                            if sign(t0(2))~=sign(t02s(inxs(ind),2)) % limit point in-between
                                x3 = obj.binary_search(0, x0, x2, t0, t02s(inxs(ind),:)', bif_par, x0_minmax);
                                if isnan(x3) % if some it failed choose one of x0 or x2 as limit points
                                    x02 = [x0, x2];
                                    if sign(t0(2))==-1
                                        [~, ind_lp] = max(x02(1, :));
                                    else
                                        [~, ind_lp] = min(x02(1, :));
                                    end
                                    LP_inxs = [LP_inxs, step_num-(2-ind_lp)]; %#ok<AGROW>
                                    if ~ (animate == false)
                                        plot(x02(1, ind_lp), obj.readout(x02(:, ind_lp)),'ko');
                                        text_offset = -0.16*(x0(1)>x02(1, ind_lp))+0.05*(x02(1, ind_lp)>x0(1));
                                        text(x02(1, ind_lp)+text_offset,obj.readout(x02(:, ind_lp)),'LP');
                                    end
                                else
                                    vars_cont(:, step_num) = x3(:, 1); % add to solutions
                                    LP_inxs = [LP_inxs, step_num]; %#ok<AGROW>
                                    step_num = step_num + 1;
                                    if ~ (animate == false)
                                        plot(x3(1),obj.readout(x3),'ko');
                                        text_offset = -0.16*(x0(1)>x3(1))+0.05*(x3(1)>x0(1));
                                        text(x3(1)+text_offset,obj.readout(x3),'LP');
                                    end
                                end
                            end
                            found_solution = true;
                            break;
                        end
                        hh = hh*2;
                    end
                    hh = min(hh, hmax);
                    if ~found_solution
                        % if no solution is found at the end - abort
                        ret_status = 1;
                        vars_cont = vars_cont(:, 1:step_num-1);
                        return;
                    end
                end
                
                % check if it converged into a previously known end-state - stop if so
                for i = 1:size(xfinals, 1)
                    if ((x2(1) >= xfinals(i, 1)) && (x0(1) < xfinals(i, 1))) || ((x2(1) <= xfinals(i, 1)) && (x0(1) > xfinals(i, 1)))
                        % if passed the Fextfinal point
                        x1 = x0 + (x2-x0)*((xfinals(i, 1)-x0(1))/(x2(1)-x0(1)));
                        if norm(x1-xfinals(i,:)') < 1e-4
                            status = 1;
                            for jj=1:numel(x0)
                                p0 = mod(p0+jj-2,numel(x0))+1;
                                [status, x3] = obj.newton(x1, bif_par, p0);
                                if (status == 0) && ~obj.out_of_bounds(x3, x0_minmax)
                                    break;
                                end
                            end
                            if (status == 0) 
                                if ~ (animate == false)
                                    addpoints(bif_profile, x3(1), obj.readout(x3));
                                    set(tang,'XData',[nan, nan]);
                                    set(tang,'YData',[nan, nan]);
                                    drawnow;
                                end
                                ret_status = 0;
                                vars_cont(:, step_num) = xfinals(i, :);
                                vars_cont = vars_cont(:, 1:step_num);
                                return;
                            end
                        end
                    elseif norm(x2-xfinals(i, :)')<h && (t0'*((xfinals(i, :)'-x2)/norm(xfinals(i, :)'-x2)) > 0)
                        % if close enough to a final point try to make the jump
                        x1 = x2; x1(p0) = xfinals(p0);
                        [status, x3] = obj.newton(x1, bif_par, p0);
                        if status ~= 0 || obj.out_of_bounds(x3, x0_minmax)
                            p0 = mod(p0,numel(x0))+1;
                            x1 = x2; x1(p0) = xfinals(p0);
                            [status, x3] = obj.newton(x1, bif_par, p0);
                            if obj.out_of_bounds(x3, x0_minmax); status = 1; end
                        end
                        if (status == 0) && (norm(x3 - xfinals(i, :)') < obj.tol_n)
                            if ~ (animate == false)
                                addpoints(bif_profile, x3(1), obj.readout(x3));
                                set(tang,'XData',[nan, nan]);
                                set(tang,'YData',[nan, nan]);
                                drawnow;
                            end
                            ret_status = 0;
                            vars_cont(:, step_num) = xfinals(i, :);
                            vars_cont = vars_cont(:, 1:step_num);
                            return;
                        end
                    end
                end
                
                % check if it went out of bounds, and previously it was not - stop if so
                if obj.out_of_bounds(x2, x0_minmax)
                    while obj.out_of_bounds(x2, x0_minmax) % iterate if needed
                        [r, c] = ind2sub([2,2], find([x2'; -x2'] < [x0_minmax(1,:); -x0_minmax(2,:)]));
                        scale = inf(1,length(r));
                        for i = 1:length(r) % in case it finished diagonally out of bounds
                            scale(i) = (x0_minmax(r(i),c(i))-x0(c(i)))/(x2(c(i))-x0(c(i)));
                        end
                        [scale_min, ind_min] = min(scale);
                        x1 = x0 + (x2-x0)*scale_min;
                        p0 = c(ind_min);
                        x1(p0) = x0_minmax(r(ind_min),c(ind_min));
                        [status, x2] = obj.newton(x1, bif_par, p0);
                        x2(p0) = x1(p0);
                        if (status == 1)
                            [status, x2] = obj.newton(x1, bif_par, mod(p0,numel(x0))+1);
                            x2(mod(p0,numel(x0))+1) = x1(mod(p0,numel(x0))+1);
                        end
                    end
                    if status == 0
                        vars_cont(:, step_num) = x2(:, 1);
                        if ~ (animate == false)
                            addpoints(bif_profile, x2(1), obj.readout(x2));
                            set(tang,'XData',[nan, nan]);
                            set(tang,'YData',[nan, nan]);
                        end
                        vars_cont = vars_cont(:, 1:step_num);
                    else
                        vars_cont = vars_cont(:, 1:step_num-1);
                    end
                    ret_status = 0;
                    return;
                end
                
                vars_cont(:, step_num) = x2(:, 1); % add to solutions
                
                if ~ (animate == false)
                    addpoints(bif_profile, x2(1), obj.readout(x2));
                    t20 = ((vars_cont(:, step_num)-vars_cont(:, step_num-1)))./norm(vars_cont(:, step_num)-vars_cont(:, step_num-1));
                    [t3, ~] = obj.tangent(x2, t20, bif_par, p0);
                    set(tang,'XData',[x2(1), x2(1)+h*t3(1)]);
                    set(tang,'YData',[obj.readout(x2), obj.readout(x2+h*t3)]);
                    if mod(step_num, 10) == 1
                        drawnow;
                    end
                end

                p0 = p2;
                x0 = x2;
            end
            ret_status = 0;
            
            function ret = vecnorm(A, p, dim)
                % finds the p-norm along the dimension DIM of A.
                ret = sum(abs(A).^p, dim).^(1/p);
            end
        end
        
        function [vars_cont_ret, LP_inxs_ret] = sort_link_profiles(obj, vars_cont, LP_inxs, order_mode)
            if nargin < 4
                order_mode = obj.order_mode;
            end
            end_points = nan(2*size(vars_cont, 2), 2);
            graph = zeros(size(end_points, 1), size(end_points, 1));
            % get data for end-points, and link them into a graph
            for i = 1:size(vars_cont, 2)
                end_points(2*i-1, :) = [vars_cont{i}(1, 1), obj.readout(vars_cont{i}(:, 1))];
                end_points(2*i, :) = [vars_cont{i}(1, end), obj.readout(vars_cont{i}(:, end))];
                graph(2*i-1, 2*i) = 1;
                graph(2*i, 2*i-1) = -1;
            end
            s_min = min(end_points(:, 1));
            s_max = max(end_points(:, 1));
            s_min_inxs = find(abs(end_points(:, 1) -s_min)<1e-5);
            s_max_inxs = find(abs(end_points(:, 1) -s_max)<1e-5);
            [~, inxs_1_order] = sort(end_points(s_min_inxs, 2), order_mode);
            %[~, inxs_2_order] = sort(end_points(s_min_inxs, 2), order_mode);
            
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
                                if mod(sum(virt_pairs(i,:)'>virt_pairs(j,:),'all'),2)==1
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
            vars_cont_ret = [];
            LP_inxs_ret = [];
            if isempty(solutions)
                disp('problem linking');
                return;
            end
            order = ceil(solutions(1, 1:2:end)/2); % ceil(indices/2) gives the vars_cont index
            flip = mod(solutions(1, 1:2:end),2)==0; % odd positions should have odd indices, otherwise flip direction within vars_cont
            for i = 1:length(flip)
                if flip(i)
                    vars_cont{order(i)} = fliplr(vars_cont{order(i)});
                    LP_inxs{order(i)} = size(vars_cont{order(i)}, 2) - LP_inxs{order(i)} + 1;
                end
            end
            vars_cont = vars_cont(order);
            LP_inxs = LP_inxs(order);
            for i = 1:length(vars_cont)
                if isempty(vars_cont_ret)
                    LP_inxs_ret = LP_inxs{i};
                    vars_cont_ret = vars_cont{i};
                else
                    LP_inxs_ret = [LP_inxs_ret, size(vars_cont_ret, 2) + 1 + LP_inxs{i}]; %#ok<AGROW>
                    vars_cont_ret = [vars_cont_ret, nan(size(vars_cont_ret,1),1), vars_cont{i}]; %#ok<AGROW>
                end
            end
        end

        function [x3] = binary_search(obj, n, x1, x2, t1, t2, bif_par, x0_minmax)
            ro = t2(1)/(t2(1)-t1(1)); % weigh closer to the point whose fp is closer to zero
            [status, x3] = obj.newton(x1*ro+x2*(1-ro), bif_par, 2);
            if (status ~= 0) 
                x3 = nan;
                return;
            end
            t3 = obj.tangent(x3, t1*ro+t2*(1-ro), bif_par, 2);
            if (abs(t3(1)) < obj.tol_n) || (n>=100) % if below tolerance threshold or recursion sufficiently deep
                if any(any([x3'; -x3'] < [x0_minmax(1,:); -x0_minmax(2,:)]))
                    % if the found LP solution is out of bounds do not add it
                    x3 = nan;
                end
                return;
            end
            if sign(t1(1))~=sign(t3(1))
                x3 = obj.binary_search(n+1, x1, x3, t1, t3, bif_par, x0_minmax);
            else
                x3 = obj.binary_search(n+1, x3, x2, t3, t2, bif_par, x0_minmax);
            end            
        end
        
        function [status, x] = newton(obj, x0, bif_par, p)
            alpha = x0(p);
            x = x0;

            it = 0;
            it_max = 100;

            while true
                if ( it_max < it )
                    status = 1;
                    return;
                end

                fx = obj.model.f_continuation(x, bif_par);
                fx(end+1,:) = x(p) - alpha; %#ok<AGROW>
                
                if any(any(isnan(fx))) || any(any(imag(fx)~=0))
                    status = 1;
                    return;
                end
                
                fx_norm = max(abs(fx));

                if(fx_norm <= obj.tol_n)
                    if any(imag(x) ~= 0)
                        status = 1;
                    else
                        status = 0;
                    end
                    return;
                end

                it = it + 1;

                fpx = obj.model.fp_continuation(x, bif_par);
                fpx(end+1, :) = 0.0; %#ok<AGROW>
                fpx(end, p) = 1.0;
                if (any(isnan(fpx(:)))) || any(imag(fpx(:))~=0) || rank(fpx) < size(fpx,1) % if singular matrix as fpx(1) approaches zero (limit point)
                    status = 2;
                    return;
                end

                dx = -fpx \ fx;
                x = x + dx;
                if any(imag(x) ~= 0)
                    status = 1;
                    return;
                end
            end
        end
        
        function [status, x2, t2, p2] = step(obj, x0, t0, bif_par, p0, h, x0_minmax)
            [t2, p2] = obj.tangent(x0, t0, bif_par, p0);

            if any(isnan(t2))
                status = 2; x2 = nan(size(x0));
                return;
            end

            x1 = x0 + h * t2;

            [status, x2] = obj.newton(x1, bif_par, p0);
            if status == 0
                if norm(x2-x0) > abs(2*h) || ((x2-x0)./norm(x2-x0))' * t0 < 0.0 %obj.out_of_bounds(x2, x0_minmax) || 
                    % if the jump is too big - bifurcation switch
                    % if out of bounds for some variable
                    % if new direction is backwards
                    status = 1;
                end
            end
        end
        
        function [t2, p2] = tangent(obj, x, t1, bif_par, p)
            fpx = obj.model.fp_continuation(x, bif_par);
            fpx(end+1, :) = 0.0;
            fpx(end, p) = 1.0;

            b = zeros(size(fpx,1), 1);
            b(end) = 1.0;
            t2 = fpx \ b;

            t2_norm = norm(t2);
            t2 = t2 / t2_norm;
            [~, p2] = max(abs(t2));
            %
            %  Especially when switching parameters, we need to make sure
            %  the sign of the tangent is chosen correctly.  In general,
            %  it should have a POSITIVE projection on the previous tangent.
            %
            if (t2' * t1 < 0.0)
                t2 = - t2;
            end 
        end

        function ret = readout(obj, x)
            % calculate readout - typically sum of some of the variables
            ret = sum(obj.readout_wts'*x(2:end,:),1);
        end    
    end

    methods (Static)
        function oob = out_of_bounds(x0, x0_minmax)
            oob = any(any([x0'; -x0'] - [x0_minmax(1,:); -x0_minmax(2,:)] < -1e-5));
        end

        function [inxs, vals] = combvec(varargin)
            size1 = length(varargin); % number of input arrays
            lens = zeros(size1,1);
            for i=1:size1
                lens(i) = length(varargin{i}); % length of each array
            end
            size2 = prod(lens); % number of combinations
            inxs = zeros(size1,size2);
            vals = zeros(size1,size2);
            ps = flipud([1;cumprod(flipud(lens(2:end)))]);
            for j=1:size2
                for i=1:size1
                    inxs(i,j) = mod(floor((j-1)/ps(i)),lens(i))+1;
                end
            end
            for i=1:size1
                vals(i,:) = varargin{i}(inxs(i,:));
            end
        end
    end
end