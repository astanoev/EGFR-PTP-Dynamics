classdef data_fitting < handle
    properties
        conditions
        model
        filename_fitting
        xx_master
        yy_master
        par_master
        par_min_master
        par_max_master
        par_shared
    end

    properties (Transient)
        n_pars_model
        inxs_pars
    end

    methods
        function obj = data_fitting(filename_fitting, conditions, model)
            if nargin<1
                filename_fitting = 'model_v10_fit_v12.5.mat';
            end
            if nargin<2
                s = load(fullfile('data',filename_fitting));
                conditions = s.conditions;
                model = s.model;
            end
            obj.filename_fitting = filename_fitting;
            obj.conditions = conditions;
            obj.model = model;
            obj.n_pars_model = numel(fieldnames(obj.model.par));
            obj.inxs_pars = containers.Map(fieldnames(obj.model.par),1:obj.n_pars_model); % dict - pars to ind
            obj.set_up_data();
            [obj.par_master, ~] = data_fitting.metropolis_hastings.load_results(fullfile('data',obj.filename_fitting));
            obj.set_up_par_bounds();
            obj.par_master(obj.par_min_master==obj.par_max_master) = obj.par_max_master(obj.par_min_master==obj.par_max_master);
            obj.set_par_master_shared_values();
        end

        function fit_data(obj, filename_fit_save, pars_gibbs)
            if nargin<3; pars_gibbs = []; end
            gibbs_burn_in = false;
            if ~isempty(pars_gibbs); gibbs_burn_in = true; end
            [obj.par_master, ss_fit, results] = data_fitting.metropolis_hastings(50).fit(@data_fitting.dose_response_est_linked, obj.par_master, obj.xx_master, ...
                    obj.yy_master, obj.par_min_master, obj.par_max_master, 'fun_args',{'model', obj.model},'par_shared',obj.par_shared,'parallel',false,...
                    'save_temp',true,'save_results_fname',filename_fit_save,'save_args',{'model', obj.model, 'conditions',obj.conditions},'track_progress',true,...
                    'gibbs_burn_in',gibbs_burn_in,'pars_gibbs',pars_gibbs); %#ok<ASGLU> 
        end

        function set_up_data(obj)
            if exist('readmatrix','file')>0; fn_read = @readmatrix; else; fn_read = @importdata; end %#ok<DLMRD> % for older versions
            folder_data = 'data'; % folder where data is stored
            % data used for fitting, compiled from all of the conditions
            obj.xx_master = nan(0,numel(obj.conditions));
            obj.yy_master = nan(0,numel(obj.conditions));
            % model used for estimating non-saturated cells, those that don't reach 100% binding, w/o ceiling for normalization
            model_satur = @(par,x) par(3).*(1-exp(-(x./par(1)).^par(2)));
            % loop through the conditions
            for i_cond = 1:numel(obj.conditions)
                condition = obj.conditions{i_cond};
                % load data for the condition
                filename_flb = strcat(condition,'_EGF_EGFR');
                filename_alpha = strcat(condition,'_alpha');
                filename_doses = strcat(condition,'_doses');
                folder_data_type = 'imaging';
                if contains(condition,'MCF7_WB'); folder_data_type = 'WB'; end
                if ~exist(fullfile(folder_data,folder_data_type,strcat(filename_flb,'.csv')), 'file') || ~exist(fullfile(folder_data,folder_data_type,strcat(filename_alpha,'.csv')), 'file') || ~exist(fullfile(folder_data,folder_data_type,strcat(filename_alpha,'.csv')), 'file')
                    fprintf('Data for condition %s not found!\n', condition);
                    continue;
                end
                EGF_EGFR = fn_read(fullfile(folder_data,folder_data_type,strcat(filename_flb,'.csv'))); % EGF/EGFR -> to convert to fraction of ligand-bound EGFR
                alpha = fn_read(fullfile(folder_data,folder_data_type,strcat(filename_alpha,'.csv'))); % fraction of phosphoryated EGFR
                doses = fn_read(fullfile(folder_data,folder_data_type,strcat(filename_doses,'.csv'))); % needed for estimating non-saturated cells
                
                if  contains(condition,'MCF7_WB')
                    % use WT doses-flb mapping to infer the EGF-EGFR fractions for WBs, especially since they might not be saturated
                    EGF_EGFR_wt = fn_read(fullfile(folder_data,'imaging','MCF7_WT_EGF_EGFR.csv'));
                    doses_wt = fn_read(fullfile(folder_data,'imaging','MCF7_WT_doses.csv'));
                    EGF_EGFR_wt_mean = mean(EGF_EGFR_wt,2);
                    EGF_EGFR_wt_mean = EGF_EGFR_wt_mean./EGF_EGFR_wt_mean(end);
                    xx = interp1(doses_wt,EGF_EGFR_wt_mean,doses);
                    xx = repmat(xx, 1, size(alpha,2))'; xx = xx(:);
                    yy = alpha'; yy = yy(:);
                    ind_nan = isnan(xx) | isnan(yy);
                    yy(ind_nan) = [];
                    xx(ind_nan) = [];
                else
                    exclude = [];
                    flb = nan(size(EGF_EGFR)); % fraction of ligand-bound
                    thr = 1;
                    % loop through cells
                    for i = 1:size(EGF_EGFR,2)
                        % fit simple model for saturation
                        par = nlinfit(doses,EGF_EGFR(:,i),model_satur,[50,1,1]);
                        % if 98% of binding is not reached, exclude cell
                        % the cells are pre-screened so there shouldn't be any non-saturated cell
                        if isempty(find(EGF_EGFR(:,i)>=0.98*par(3), 1))
                            disp('non-saturated cell excluded');
                            exclude = [exclude,i]; %#ok<AGROW>
                        end
                        % find first dose where saturation is reached - discard the rest of the data
                        inx = find(EGF_EGFR(:,i)./(max(EGF_EGFR(:,i)))>=thr,1);
                        % calculate ligand-bound fraction accordingly
                        flb(1:inx,i) = EGF_EGFR(1:inx,i)./EGF_EGFR(inx,i);
                    end
            
                    flb(:,exclude) = nan;
                    alpha(:,exclude) = nan;
                
                    % collect data points (flb-alpha pairs) from all of the cells
                    xx = flb(:);
                    yy = alpha(:);
                    ind_nan = isnan(xx) | isnan(yy);
                    yy(ind_nan) = [];
                    xx(ind_nan) = [];
                    % sort the points according to x-value
                    [xx, xx_inxs] = sort(xx);
                    yy = yy(xx_inxs);
                    xy_data = [xx'; yy'];
                    % calculate moving median
                    mm = data_fitting.moving_median(xy_data, 2, 20);
                    xx = mm(1,:);
                    yy = mm(2,:);
                    % remove offset of y-data
                    if ~strcmp(condition, 'MCF7_PTPRG_KO') && ~strcmp(condition, 'MCF7_PTPN2_KO')
                        % do not subtract initial point only for the KO cases with real offset
                        yy = yy -yy(1);
                    end
                end
                
                % if the xx_master array has grown, remove trailing nan's
                if size(obj.xx_master,1)<numel(xx)
                    obj.xx_master((size(obj.xx_master,1)+1):numel(xx),:) = nan;
                end
                if size(obj.yy_master,1)<numel(yy)
                    obj.yy_master((size(obj.yy_master,1)+1):numel(yy),:) = nan;
                end
                % store data for condition in master matrices
                obj.xx_master(1:length(xx),i_cond) = xx';
                obj.yy_master(1:length(yy),i_cond) = yy';
            end
        end

        function set_up_par_bounds(obj)
            % bounds of the parameters
            obj.par_min_master = [];
            obj.par_max_master = [];
            % loop through the conditions
            for j = 1:numel(obj.conditions)
                % initialize bounds to (0, infinity)
                par_min = zeros(1,obj.n_pars_model);
                par_max = inf(1,obj.n_pars_model);
                % fix k2/k1 parameter to 0.087, given the estimated measurement of alpha_ox=0.08 from the NOX-KO condition
                % k21 = alpha_ox/(1-alpha_ox) (the formula is derived from the ODE, taking alpha_ox=PTPRGi)
                par_min(obj.get_par_inxs('k21')) = 0.087;
                par_max(obj.get_par_inxs('k21')) = par_min(obj.get_par_inxs('k21'));
                if j==find(strcmp(obj.conditions,'MCF7_PTPRG_KO'),1) || j==find(strcmp(obj.conditions,'MCF7_p22_KO'),1)
                    par_max(obj.get_par_inxs('b1')) = 0;
                end
        
                % add parameter bounds to master matrix
                obj.par_min_master = [obj.par_min_master, par_min];
                obj.par_max_master = [obj.par_max_master, par_max];
            end
        end
        
        function set_par_master_shared_values(obj)
            ind_wt = find(strcmp(obj.conditions,'MCF7_WT'),1);
            if numel(obj.conditions) ~= numel(obj.par_master)/obj.n_pars_model
                fprintf('adding %d extra condition(s) to parameter set\n', numel(obj.conditions)-numel(obj.par_master)/obj.n_pars_model);
                % if new conditions are introduced, copy the values from WT condition
                for i_n = (numel(obj.par_master)/obj.n_pars_model):(numel(obj.conditions)-1)
                    obj.par_master(i_n*obj.n_pars_model +(1:obj.n_pars_model)) = obj.par_master(ind_wt-1 +1:obj.n_pars_model);
                end
            end
        
            % which parameters are shared across conditions - par_shared(a) = b -> means a-th value will be inherited from b-th value
            obj.par_shared = zeros(size(obj.par_master));
        
            for i=1:numel(obj.conditions)
                if i>1
                    obj.par_shared(obj.get_par_inxs('e1','e2','e3','e4','a1','a2','a3','a4','k21') +(i-1)*obj.n_pars_model) = obj.get_par_inxs('e1','e2','e3','e4','a1','a2','a3','a4','k21');
                end
            end
            
            inxs_n2 = find(strncmp(obj.conditions,'MCF7_PTPN2',9));
            for ind_n2 = inxs_n2
                obj.par_shared(obj.get_par_inxs('g1','g3','b1') +(ind_n2-1)*obj.n_pars_model) = obj.get_par_inxs('g1','g3','b1');
            end
        
            inxs_rg = find(strncmp(obj.conditions,'MCF7_PTPRG',9));
            for ind_rg = inxs_rg
                obj.par_shared(obj.get_par_inxs('g2','g4') +(ind_rg-1)*obj.n_pars_model) = obj.get_par_inxs('g2','g4');
            end
            
            ind_wb = find(strcmp(obj.conditions,'MCF7_WB'),1);
            if ~isempty(ind_wb)
                obj.par_shared(obj.get_par_inxs('g1','g2','g3','g4','b1') +(ind_wb-1)*obj.n_pars_model) = 0;
            end
        
            ind_nox = find(strcmp(obj.conditions,'MCF7_p22_KO'),1);
            if ~isempty(ind_nox)
                obj.par_shared(obj.get_par_inxs('g1','g3') +(ind_nox-1)*obj.n_pars_model) = obj.get_par_inxs('g1','g3');
                obj.par_shared(obj.get_par_inxs('b1') +(ind_nox-1)*obj.n_pars_model) = 0;
            end
        
            inxs_shared = find(obj.par_shared > 0);
            for i = 1:numel(inxs_shared)
                i_s = inxs_shared(i);
                obj.par_master(i_s) = obj.par_master(obj.par_shared(i_s));
            end
        end

        function ret = get_par_inxs(obj, varargin)
            ret = cell2mat(values(obj.inxs_pars,{varargin{:}})); %#ok<CCAT1> 
        end
    end

    methods (Static)
        function y_est = dose_response_est_linked(par, x_data, varargin)
            p = inputParser;
            addParameter(p,'model',[]);
            parse(p,varargin{:});
            model = p.Results.model;
            if isempty(model)
                model = models.egfr_ptprg_model;
            end
            fn = fieldnames(model.par);
            y_est = nan(size(x_data));
            for j=1:size(x_data,2)
                model.set_params(par((1:length(fn))+(j-1)*length(fn)));
                y_est_j = data_fitting.dose_response_est_single(model, x_data(:,j));
                y_est(1:numel(y_est_j),j) = y_est_j;
            end
        end

        function [y_ret, rmse] = dose_response_est_single(model, x_data, y_data)
            ct = dynamics.continuation(model);
            ct.calc_profile('EGF_EGFRtt', [0, 1], [0, 1]);

            x_data_nn = x_data(~isnan(x_data));
            [xx, ~, xx_ind] = unique(x_data_nn);

            %x = ct.prof_bif.vars_cont(1,:);
            %y = ct.readout(ct.prof_bif.vars_cont);
            branches_stable = ct.prof_bif.get_stable_branches(1);
            x = branches_stable(1,:);
            y = ct.readout(branches_stable);
            ind_nan = find(isnan(x),1,'first');
            if ~isempty(ind_nan)
                ind_21 = find(x((ind_nan+1):end)>x(ind_nan-1),1,'first') + ind_nan;
                y_21 = interp1(x((ind_nan+1):end), y((ind_nan+1):end), x(ind_nan-1)+1e-10, 'PCHIP', nan);
                yy = interp1([x(1:(ind_nan-1)),x(ind_nan-1)+1e-10,x(ind_21:end)], [y(1:(ind_nan-1)),y_21,y(ind_21:end)], xx, 'PCHIP', nan);
            else
                yy = interp1(x, y, xx, 'PCHIP', nan);
            end
            y_ret = nan(size(x_data));
            y_ret(~isnan(x_data)) = yy(xx_ind);
            rmse = nan;
            if nargin>2
                y_data_nn = y_data(~isnan(y_data));
                rmse = sqrt(sum((y_data_nn -yy(xx_ind)).^2)/numel(~isnan(y_data)));
            end
        end

        function xy_ret = moving_median(xy_data, flag, window_size)
            if nargin<2; flag = 1; end % flag=1 is mean+sd, flag=2 is median+mad
            if nargin<3; window_size = 10; end
            [~, ind_sort] = sort(xy_data(1, :));
            xy_data = xy_data(:, ind_sort);
            xy_ret = zeros(3, size(xy_data,2)-window_size+1);
            for i = 1:(size(xy_data,2)-window_size+1)
                xy_data_window = xy_data(:,i:(i+window_size-1));
                if flag == 2
                    xy_ret(1, i) = median(xy_data_window(1, :));
                    xy_ret(2, i) = median(xy_data_window(2, :));
                    xy_ret(3, i) = mad(xy_data_window(2, :),1);
                elseif flag == 1
                    xy_ret(1, i) = mean(xy_data_window(1, :));
                    xy_ret(2, i) = mean(xy_data_window(2, :));
                    xy_ret(3, i) = std(xy_data_window(2, :));
                end
            end
        end
    end
end