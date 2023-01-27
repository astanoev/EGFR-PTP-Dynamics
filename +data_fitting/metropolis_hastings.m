classdef metropolis_hastings < handle
    properties
        nr_random_walks = 1;
        n_iter = 300;
        n_bigiter = 100;
        continue_mcmc = false;
        iter_update_sigma = 50; % saving and correcting periodically
        sigma_init = 0.01; %initial sigma - later adjusted to match the target acceptance ratio
        sigma_min = 1e-5;
        slope_thr = -1e-5;
        k_burn_in = 15;
        steps_converge = 0;
        target_acceptance_ratio = .4; 
        error_function = @obj.sum_of_squares;
        rng_seed;
        save_temp = false;
        record_ydata_fun = false;
    end
    
    methods
        function obj = metropolis_hastings(n_iter)
            if nargin>0
                obj.n_iter = n_iter;
            end
        end
        
        function [par_fit, ss_fit, results] = fit(obj, fun, par_init, xdata, ydata, varargin)
           %% input is the same as the standard fitting methods:
            % fun - f'n for fitting; par_init - initial guess for the parameter set;
            % xdata - input data (MxN), each column is M-dimensional single input vector
            % ydata - output data (LxN), each column is L-dimensional
            % single output vector (in principle L can be different than M)
            p = inputParser;
            addOptional(p,'par_lb',[]); % lower bounds of the parameters (array)
            addOptional(p,'par_ub',[]); % upper bounds of the parameters
            addOptional(p,'par_shared',[]); % sharing of parameters
            addOptional(p,'data_weights',ones(1,size(xdata,2))); % weights of the data points (N-dim)
            addOptional(p,'error_function',@obj.sum_of_squares); % log-probability
            addParameter(p,'fun_args',{}); % run multiple parallel executions
            addParameter(p,'parallel',false); % run multiple parallel executions
            addParameter(p,'continue_fit',struct()); % continue a previous run - structure
            addParameter(p,'n_iter',obj.n_iter); % no. of iterations
            addParameter(p,'gibbs_burn_in',false); % sample individually the parameters
            addParameter(p,'pars_gibbs',[]); % loop through selected parameters only during gibbs sampling
            addParameter(p,'save_temp',obj.save_temp); % save execution file at output
            addParameter(p,'save_results_fname',''); % save execution file at output
            addParameter(p,'save_args',{}); % additional arguments to save in file
            addParameter(p,'track_progress',false); % plot ss, fit and residuals with the execution
            addParameter(p,'save_progress_png',false); % save progress iterations as images
            parse(p,varargin{:});
            par_lb = p.Results.par_lb;
            par_ub = p.Results.par_ub;
            par_shared = p.Results.par_shared;
            data_weights = p.Results.data_weights;
            obj.error_function = p.Results.error_function;
            fun_args = p.Results.fun_args;
            save_args = p.Results.save_args;
            obj.save_temp = p.Results.save_temp;
            if isempty(par_lb); par_lb = -inf(size(par_init)); end
            if isempty(par_ub); par_ub = inf(size(par_init)); end
            if isempty(par_shared); par_shared = zeros(size(par_init)); end
            if size(par_lb,2)<size(par_lb,1); par_lb = par_lb'; end
            if size(par_ub,2)<size(par_ub,1); par_ub = par_ub'; end
            if size(par_shared,2)<size(par_shared,1); par_shared = par_shared'; end
            if size(par_init,2) <size(par_init,1); par_init = par_init'; end
            obj.n_iter = p.Results.n_iter;
            continue_fit = p.Results.continue_fit;
            if p.Results.parallel
                pp1 = gcp('nocreate');
                if isempty(pp1)
                    try
                        pp1 = parpool;
                        obj.nr_random_walks = pp1.NumWorkers;
                    catch ex
                        disp(strcat('cannot open a parpool because: ', ex.message));
                        disp('operating in single mode instead');
                        obj.nr_random_walks = 1;
                    end
                else
                    obj.nr_random_walks = pp1.NumWorkers;
                end
            end
            
            par_fixed = par_lb==par_ub;
            par_free = find(~par_fixed & par_shared==0);
            if p.Results.gibbs_burn_in && ~isempty(p.Results.pars_gibbs)
                par_free = p.Results.pars_gibbs;
            end
            par_init(par_fixed) = par_lb(par_fixed); % set parameter if pre-fixed
            par_init(par_shared>0) = par_init(par_shared(par_shared>0)); % copy values for shared par's
            par_fit = par_init;
            
            ydata_est = fun(par_init, xdata, fun_args{:});
            ss_fits = obj.error_function(ydata_est, ydata, data_weights);
            
            if p.Results.track_progress
                [~, ~, pp1, pp1s, ~, pp2, pp2s, ~, pp3] = obj.progress_plot(ss_fits, xdata, ydata, fun, par_init, fun_args, p.Results.save_progress_png);
            end
            
            par_lb_master = par_lb;
            par_ub_master = par_ub;
            
            i_burn_in = inf;
            filename_temp = strcat(tempname('data\temp'),'.mat');
            results = struct();
            results.burn_in = cell(0);
            for i = 1:obj.n_bigiter
                if p.Results.track_progress; fprintf([num2str(i),' ']); end
                
                if isinf(i_burn_in) && p.Results.gibbs_burn_in
                    par_lb = par_fit; par_ub = par_fit;
                    par_lb(par_free(mod(i-1,numel(par_free))+1)) = par_lb_master(par_free(mod(i-1,numel(par_free))+1));
                    par_ub(par_free(mod(i-1,numel(par_free))+1)) = par_ub_master(par_free(mod(i-1,numel(par_free))+1));
                else
                    par_lb = par_lb_master; par_ub = par_ub_master;
                end
                
                if obj.nr_random_walks > 1
                    [par_fit, ss_fit, sigma_last, results_fit] = obj.fit_parallel(fun, par_fit, xdata, ydata, fun_args, par_lb, par_ub, par_shared, data_weights, continue_fit);
                else
                    [par_fit, ss_fit, sigma_last, results_fit] = obj.fit_single(fun, par_fit, xdata, ydata, fun_args, par_lb, par_ub, par_shared, data_weights, continue_fit);
                end
                ss_fits = [ss_fits, ss_fit]; %#ok<AGROW>
                if isinf(i_burn_in)
                    % if during burn-in period the parallel search yields
                    % the same results k-times as in the previous iterations 
                    % the burn-in period is done
                    if length(ss_fits)>=obj.k_burn_in && external_tools.weightedfit([1:obj.k_burn_in; ss_fits((end-(obj.k_burn_in-1)):end)]').slope > obj.slope_thr
                        i_burn_in = i;
                        results.converged = results_fit;
                        continue_fit = results_fit;
                        if p.Results.track_progress
                            fprintf('burn-in finished ');
                            pp1s.XData = [i-.5, i-.5]; pp1s.YData = mean(ss_fits(end-1:end))+.5*1e-4*[-1,1];
                        end
                        results.burn_in = [];
                    else
                        results.burn_in{length(results.burn_in)+1} = results_fit;
                        continue_fit.sigmas = sigma_last;
                        if length(results.burn_in)>1
                            results.burn_in{length(results.burn_in)-1} = [];
                        end
                    end
                else
                    results.converged = results_fit;
                    continue_fit = results_fit;
                end
                if obj.save_temp
                    struct_save = struct('results',results);
                    for i_sa=1:numel(save_args)/2
                        struct_save.(save_args{2*i_sa-1}) = save_args{2*i_sa};
                    end
                    save(filename_temp, '-struct', 'struct_save');
                end
                if p.Results.track_progress
                    [pp1, pp1s, pp2, pp2s, pp3] = obj.progress_plot_update(ss_fits, results_fit, pp1, pp1s, pp2, pp2s, pp3, i, i_burn_in, xdata, ydata, fun, par_fit, fun_args, p.Results.save_progress_png);
                end
                if i >= i_burn_in +obj.steps_converge
                    break;
                end
            end
            if p.Results.track_progress; fprintf('\n'); end
            if ~isempty(p.Results.save_results_fname)
                if obj.save_temp
                    movefile(filename_temp, fullfile('data\temp', strcat(p.Results.save_results_fname,'.mat')), 'f');
                else
                    var_results = whos('results');
                    if var_results.bytes/2^30>0.5
                        % if larger than 500mb, delete the burn-in part
                        if isfield(results,'converged')
                            results.burn_in = [];
                        end
                    end
                    save(fullfile('data\temp', strcat(p.Results.save_results_fname,'.mat')), 'results', '-v7.3');
                end
            elseif obj.save_temp
                delete(filename_temp);
            end
        end

        function [par_fit, ss_fit, results_fit] = fit2(obj, fun, par_init, xdata, ydata, varargin)
            %% input is the same as the standard fitting methods:
            % fun - f'n for fitting; par_init - initial guess for the parameter set;
            % xdata, ydata - data used for fitting
            p = inputParser;
            addOptional(p,'par_lb',[]); % lower bounds of the parameters (array)
            addOptional(p,'par_ub',[]); % upper bounds of the parameters
            addParameter(p,'parallel',false); % run multiple parallel executions
            addParameter(p,'continue_fit',struct()); % continue a previous run - structure
            addParameter(p,'n_iter',obj.n_iter); % no. of iterations
            parse(p,varargin{:});
            par_lb = p.Results.par_lb;
            par_ub = p.Results.par_ub;
            if isempty(par_lb); par_lb = -inf(size(par_init)); end
            if isempty(par_ub); par_ub = inf(size(par_init)); end
            if size(par_init,2) <size(par_init,1); par_init = par_init'; end
            obj.n_iter = p.Results.n_iter;
            if p.Results.parallel
                pp = gcp('nocreate');
                if isempty(pp)
                    try
                        pp = parpool;
                        obj.nr_random_walks = pp.NumWorkers;
                    catch ex %#ok<NASGU>
                        obj.nr_random_walks = 1;
                    end
                else
                    obj.nr_random_walks = pp.NumWorkers;
                end
            end
            if obj.nr_random_walks > 1
                [par_fit, ss_fit, results_fit] = obj.fit_parallel(fun, par_init, xdata, ydata, par_lb, par_ub, p.Results.continue_fit);
            else
                [par_fit, ss_fit, results_fit] = obj.fit_single(fun, par_init, xdata, ydata, par_lb, par_ub, p.Results.continue_fit);
            end
        end

        function [par_bestfit, ss_bestfit, sigma_last, results_fit] = fit_parallel(obj, fun, par_init, xdata, ydata, fun_args, par_lb, par_ub, par_shared, data_weights, continue_fit)
            if nargin<9 || (~iscell(continue_fit) && isempty(setdiff(fieldnames(continue_fit),{'sigmas'})))
                sigmas = [];
                if any(strcmp(fieldnames(continue_fit),'sigmas'))
                    sigmas = continue_fit.sigmas;
                end
                continue_fit = cell(obj.nr_random_walks,1);
                continue_fit(:) = {struct()};
                for cpu_id = 1:obj.nr_random_walks
                    seed = randi(floor(intmax/10)) + cpu_id;
                    continue_fit{cpu_id}.rng_init = rng(seed,'Threefry');
                    continue_fit{cpu_id}.rng_state = continue_fit{cpu_id}.rng_init;
                    if ~isempty(sigmas); continue_fit{cpu_id}.sigmas = sigmas; end
                end
            end 
            spmd
            %for labindex=1:obj.nr_random_walks
                cpu_id = labindex;
                cf = continue_fit{cpu_id};
                [~, ~, ~, results_fit_single] = obj.fit_single(fun, par_init, xdata, ydata, fun_args, par_lb, par_ub, par_shared, data_weights, cf);
            end
            results_fit = cell(obj.nr_random_walks,1);
            results_fit(:) = {results_fit_single{:}}; %#ok<CCAT1>
            ss_bestfit = inf;
            par_bestfit = [];
            for cpu_id=1:obj.nr_random_walks
                [ss_min_cpu_id, ind_min] = min(results_fit{cpu_id}.ss_mat);
                if ss_min_cpu_id < ss_bestfit
                    ss_bestfit = ss_min_cpu_id;
                    par_bestfit = results_fit{cpu_id}.par_mat(ind_min, :);
                    sigma_last = results_fit{cpu_id}.sigmas(end);
                end
            end
        end
        
        function [par_bestfit, ss_bestfit, sigma_last, results_fit] = fit_single(obj, fun, par_init, xdata, ydata, fun_args, par_lb, par_ub, par_shared, data_weights, continue_fit)
            if nargin<9; continue_fit = struct(); end 
            par_mat = zeros(obj.n_iter, length(par_init)); % record sampled parameter sets
            if obj.record_ydata_fun; ydata_fun = zeros(obj.n_iter, size(ydata,1), size(ydata,2)); end % record sampled outputs
            ss_mat = zeros(obj.n_iter,1); % records sum-of-squares for parameter sets
            beta_mat = zeros(obj.n_iter,1); % beta values - to update?
            if ~any(strcmp(fieldnames(continue_fit),'par_mat'))
                par = par_init;
            else
                par = continue_fit.par_mat(end,:);
            end
            if ~any(strcmp(fieldnames(continue_fit),'sigmas'))
                sigma = obj.sigma_init; % initial sigma - later adjusted to match the target acceptance ratio
            else
                sigma = continue_fit.sigmas(end);
            end
            if ~any(strcmp(fieldnames(continue_fit),'rng_state'))
                rng_init = rng('shuffle','combRecursive');
                rng(rng_init);
            else
                rng_init = continue_fit.rng_init;
                rng(continue_fit.rng_state);
            end
            par_mat(1,:) = par;
            yf = fun(par, xdata, fun_args{:});
            if obj.record_ydata_fun; ydata_fun(1,:,:) = yf; end
            ss = obj.error_function(yf, ydata, data_weights);
            ss_mat(1,:) = ss;
            beta = obj.update_beta(ss, size(xdata,2));
            beta_mat(1,:) = beta;
            sigmas = zeros(obj.n_iter/obj.iter_update_sigma,1);
            
            for iter=2:obj.n_iter
                par_new = obj.generate_new_params(sigma, par, par_lb, par_ub, par_shared);
                try
                    yf = fun(par_new, xdata, fun_args{:});
                catch ex
                    save('data\temp\par_problem.mat','par_new');
                    disp('parameter set problem');
                    %throw(ex);
                    par_new = par;
                end
                ss_new = obj.error_function(yf,ydata,data_weights);
                beta_new = obj.update_beta(ss_new, size(xdata,2));
                % acceptance criteria 
                if (ss_new < ss) || ((ss_new >= ss) && (log(rand) < obj.log_acceptance_ratio(ss_new, ss, par_new, par, beta_new, sigma, par_lb, par_ub, par_shared)))
                    ss = ss_new;
                    par = par_new;
                    beta = beta_new;
                end
                par_mat(iter,:) = par;
                ss_mat(iter,:) = ss;
                beta_mat(iter,:) = beta;
                if obj.record_ydata_fun; ydata_fun(iter,:,:) = yf; end
                
                if mod(iter,obj.iter_update_sigma)==0
                    % update sigma periodically to match the acceptance
                    % ratio of the underlying distribution
                    sigma = obj.update_sigma(sigma, par_mat, iter, obj.iter_update_sigma);
                    sigmas(iter/obj.iter_update_sigma) = sigma;
                end
            end
            if ~any(strcmp(fieldnames(continue_fit),'par_mat'))
                results_fit = struct();
                results_fit.par_mat = par_mat;
                if obj.record_ydata_fun; results_fit.ydata_fun = ydata_fun; end
                results_fit.ss_mat = ss_mat;
                results_fit.beta_mat = beta_mat;
                results_fit.sigmas = sigmas;
                results_fit.rng_init = rng_init;
            else
                results_fit = continue_fit;
                results_fit.par_mat = cat(1,results_fit.par_mat,par_mat);
                if obj.record_ydata_fun; results_fit.ydata_fun = cat(1,results_fit.ydata_fun,ydata_fun); end
                results_fit.ss_mat = cat(1,results_fit.ss_mat,ss_mat);
                results_fit.beta_mat = cat(1,results_fit.beta_mat,beta_mat);
                results_fit.sigmas = cat(1,results_fit.sigmas,sigmas);
            end
            results_fit.rng_state = rng();
            [ss_bestfit, ind_min] = min(results_fit.ss_mat);
            par_bestfit = results_fit.par_mat(ind_min, :);
            sigma_last = sigmas(end);
            results_fit.ss_bestfit = ss_bestfit;
            results_fit.par_bestfit = par_bestfit;
        end
        
        function results_fit = optimize_pars(obj, fun, inxs_pars, par_init, xdata, ydata, varargin)
            n_perms = 100;
            fun_args = varargin{find(strcmp(varargin,'fun_args'),1)+1};
            model = fun_args{find(strcmp(fun_args,'model'),1)+1};
            par_names = fieldnames(model.par);
            fig = figure; ax = axes('parent',fig); hold(ax,'on'); box(ax,'on');
            for i_p = 1:n_perms
                perm = randperm(numel(inxs_pars));
                inxs_perm = inxs_pars(perm);
                if i_p==1
                    results_fit = varargin{find(strcmp(varargin,'continue_fit'),1)+1};
                    ss_prev = results_fit.ss_mat(end);
                    results_fit.par_ind = zeros(size(results_fit.ss_mat));
                    ind_init = numel(results_fit.ss_mat)+1;
                else
                    ss_prev = ss_last;
                end
                for i = 1:numel(inxs_perm)
                    ind_par = inxs_perm(i);
                    varargin{find(strcmp(varargin,'continue_fit'),1)+1} = results_fit;
                    results_fit = obj.optimize_single_par(fun, ind_par, par_init, xdata, ydata, varargin{:});
                    ind_init_par = find(results_fit.par_ind~=results_fit.par_ind(end),1,'last')+1;
                    plot(ax, ind_init_par-ind_init+(0:(numel(results_fit.ss_mat)-ind_init_par)), results_fit.ss_mat(ind_init_par:end));
                    txt_ypos = results_fit.ss_mat(ind_init_par) +rand*(results_fit.ss_mat(max([ind_init_par-200,ind_init]))-results_fit.ss_mat(ind_init_par));
                    text(ax, ind_init_par-ind_init+(numel(results_fit.ss_mat)-ind_init_par)/2, txt_ypos, par_names(ind_par), 'VerticalAlignment','bottom','FontSize',12);
                    xlim(ax, [max([numel(results_fit.ss_mat)-ind_init-200,0]), numel(results_fit.ss_mat)-ind_init]);
                    ylim(ax, [results_fit.ss_mat(end), results_fit.ss_mat(max([ind_init_par-200,ind_init]))]);
                    drawnow;
                    par_init = results_fit.par_mat(end,:);
                end
                ss_last = results_fit.ss_mat(end);
                if ss_prev-ss_last < 1e-8
                    break;
                end
            end
        end

        function results_fit = optimize_single_par(obj, fun, ind_par, par_init, xdata, ydata, varargin)
            p = inputParser;
            addOptional(p,'par_lb',[]); % lower bounds of the parameters (array)
            addOptional(p,'par_ub',[]); % upper bounds of the parameters
            addOptional(p,'par_shared',[]); % sharing of parameters
            addOptional(p,'data_weights',ones(1,size(xdata,2))); % weights of the data points (N-dim)
            addOptional(p,'error_function',@obj.sum_of_squares); % log-probability
            addParameter(p,'fun_args',{}); % run multiple parallel executions
            addParameter(p,'continue_fit',struct()); % continue a previous run - structure
            parse(p,varargin{:});
            par_lb = p.Results.par_lb;
            par_ub = p.Results.par_ub;
            par_shared = p.Results.par_shared;
            data_weights = p.Results.data_weights;
            obj.error_function = p.Results.error_function;
            fun_args = p.Results.fun_args;
            continue_fit = p.Results.continue_fit;

            par_mat = zeros(obj.n_iter, length(par_init)); % record sampled parameter sets
            ss_mat = zeros(obj.n_iter,1); % records sum-of-squares for parameter sets
            par = par_init;
            if ~any(strcmp(fieldnames(continue_fit),'sigmas'))
                sigma = obj.sigma_init; % initial sigma - later adjusted to match the target acceptance ratio
            else
                sigma = continue_fit.sigmas(end);
            end
            if ~any(strcmp(fieldnames(continue_fit),'rng_state'))
                rng_init = rng('shuffle','combRecursive');
                rng(rng_init);
            else
                rng_init = continue_fit.rng_init;
                rng(continue_fit.rng_state);
            end
            par_mat(1,:) = par;
            yf = fun(par, xdata, fun_args{:});
            ss = obj.error_function(yf, ydata, data_weights);
            ss_mat(1,:) = ss;
            gamma_0 = 10;
            gamma = gamma_0;
            for iter=2:obj.n_iter
                if iter==2
                    par_new = obj.generate_new_params(sigma, par, par_lb, par_ub, par_shared);
                    par_new(setdiff(1:numel(par),ind_par)) = par(setdiff(1:numel(par),ind_par));
                else
                    par_new(ind_par) = par(ind_par) -gamma*dF;
                end
                par_new(par_shared==ind_par) = par_new(ind_par);
                if par_new(ind_par)<par_lb(ind_par) || par_new(ind_par)>par_ub(ind_par)
                    gamma = gamma/2;
                    ss_mat(iter,:) = nan;
                    par_mat(iter,:) = nan;
                    continue;
                end
                try
                    yf = fun(par_new, xdata, fun_args{:});
                catch ex
                    save('data\temp\par_problem.mat','par_new');
                    disp('parameter set problem');
                    gamma = gamma/2;
                    ss_mat(iter,:) = nan;
                    par_mat(iter,:) = nan;
                    continue;
                end
                ss_new = obj.error_function(yf,ydata,data_weights);
                dF = (ss_new -ss)/sign(par_new(ind_par) -par(ind_par));
                if ss_new < ss
                    ss = ss_new;
                    par = par_new;
                    gamma = gamma_0;
                elseif iter>2
                    disp([num2str(ind_par),' ',num2str(iter),' gamma reduced']);
                    gamma = gamma/2;
                end
                par_mat(iter,:) = par;
                ss_mat(iter,:) = ss;
                if iter>2 && (abs(dF) < 1e-10 || gamma<1e-2)
                    break;
                end
            end
            par_mat((iter+1):end,:) = [];
            ss_mat((iter+1):end,:) = [];
            par_mat(isnan(ss_mat),:) = [];
            ss_mat(isnan(ss_mat),:) = [];
            if ~any(strcmp(fieldnames(continue_fit),'par_mat'))
                results_fit = struct();
                results_fit.par_mat = par_mat;
                results_fit.ss_mat = ss_mat;
                results_fit.rng_init = rng_init;
                results_fit.par_ind = ind_par*ones(numel(ss_mat),1);
            else
                results_fit = continue_fit;
                results_fit.par_mat = cat(1,results_fit.par_mat,par_mat);
                results_fit.ss_mat = cat(1,results_fit.ss_mat,ss_mat);
                results_fit.par_ind = cat(1,results_fit.par_ind,ind_par*ones(numel(ss_mat),1));
            end
            results_fit.rng_state = rng();
            [ss_bestfit, ind_min] = min(results_fit.ss_mat);
            par_bestfit = results_fit.par_mat(ind_min, :);
            results_fit.ss_bestfit = ss_bestfit;
            results_fit.par_bestfit = par_bestfit;
        end

        function lar = log_acceptance_ratio(obj, ss_new, ss_old, par_new, par_old, beta, sigma, par_lb, par_ub, par_shared)
            lar = obj.log_likelihood_ratio(ss_new, ss_old, beta) + ...
                log(obj.correct_negative_rejection(par_new, par_old, sigma, par_lb, par_ub, par_shared)) + ...
                obj.log_proposal_ratio(par_new, par_old, sigma);
        end
        
        function llr = log_likelihood_ratio(obj, ss_new, ss_old, beta) %#ok<INUSL>
            %% ratio of stationary distributions 
            llr = -beta*.5*(ss_new-ss_old);
        end
        
        function cnr = correct_negative_rejection(obj, par_new, par_old, sigma, par_lb, par_ub, par_shared) %#ok<INUSL>
            %% for nonsymmetric MH - correct for the nonsymmetric MVN jumps
            % incorrect acceptance ratio is introduced when out-of-bound parameter sets are rejected
            % normalizing functions (normal CDFs) for the jump between the parameters do not
            % cancel out when there are bounds
            % Multivariate CDF around the par set is calculated for the
            % given par bounds - closer to the bounds it drops below 1
            % https://darrenjw.wordpress.com/2012/06/04/metropolis-hastings-mcmc-when-the-proposal-and-target-have-differing-support/
            par_fixed_shared = (par_lb==par_ub) | (par_shared>0); % do not take into account fixed or shared parameters
            if any(par_lb(~par_fixed_shared)>-inf) || any(par_ub(~par_fixed_shared)<inf) % if there are bounds
                % cnr = mvncdf(par_lb(~par_fixed_shared), par_ub(~par_fixed_shared), par_old(~par_fixed_shared), diag((sigma*par_old(~par_fixed_shared)).^2))/mvncdf(par_lb(~par_fixed_shared), par_ub(~par_fixed_shared), par_new(~par_fixed_shared), diag((sigma*par_new(~par_fixed_shared)).^2));
                % since Sigma is diagonal, the individual normal distrs are
                % independent, hence the cdfs can be calculated separately
                % and multiplied, instead of running mvncdf (which is
                % limited to 25 parameters and a bit slower)
                % normcdf(([0, 0.45]-par_old(5))./((sigma*par_old(5))),0,1) is
                % transformed from normcdf([0, 0.45],par_old(5),(sigma*par_old(5)))
                p_old_lb = normcdf((par_lb(~par_fixed_shared)-par_old(~par_fixed_shared))./((sigma*par_old(~par_fixed_shared))));
                p_old_ub = normcdf((par_ub(~par_fixed_shared)-par_old(~par_fixed_shared))./((sigma*par_old(~par_fixed_shared))));
                mvncdf_old = prod(p_old_ub-p_old_lb);
                p_new_lb = normcdf((par_lb(~par_fixed_shared)-par_new(~par_fixed_shared))./((sigma*par_new(~par_fixed_shared))));
                p_new_ub = normcdf((par_ub(~par_fixed_shared)-par_new(~par_fixed_shared))./((sigma*par_new(~par_fixed_shared))));
                mvncdf_new = prod(p_new_ub-p_new_lb);
                cnr = mvncdf_old/mvncdf_new;
            else % otherwise normalizing functions are equal
                cnr = 1;
            end
        end
        
        function lpr = log_proposal_ratio(obj, par_new, par_old, sigma) %#ok<INUSL>
            %% correct for State-Dependent Proposal Scalings
            % as in 'Optimal Proposal Distributions and Adaptive MCMC' by Jeffrey S. Rosenthal
            lpr = sum(log(abs(par_old))) -sum(log(abs(par_new))) -.5/sigma^2*sum((par_old-par_new).^2.*(1./par_new.^2-1./par_old.^2));
        end
        
        function beta = update_beta(obj, ss, n) %#ok<INUSL>
            %% to double-check implementation!!!
            %% error variance (of data) - Gibbs hyperparameter - update using Gamma conjugate prior
            % equivallent to sampling of scaled inverse chi-squared distribution
            % https://dovlab.wordpress.com/2013/02/07/inferring-the-error-variance-in-metropolis-hastings-mcmc/
            % https://rsmith.math.ncsu.edu/MA540_F20/LECTURES/Lecture6.pdf
            a = 0; b = 0;
            beta=gamrnd(a+.5*n,b+.5*ss);
        end
        
        function beta = calc_beta(obj, xdata, ydata) %#ok<INUSL>
            n_bins = 10;
            [~,~,bin] = histcounts(xdata,n_bins);
            beta_bins = accumarray(bin(:),ydata,[],@std);
            beta = 1/mean(beta_bins)^2;
        end

        function sigma = update_sigma(obj, sigma, params, t, n_t)
            %% proposal scaling - adjust sigma to yield the target acceptance ratio
            % similar as in http://probability.ca/jeff/ftpdir/galinart.pdf
            % or http://www.utstat.toronto.edu/wordpress/WSFiles/technicalreports/0610.pdf
            n_u = size(unique(params((t-n_t+1):t,:),'rows'),1);
            acceptance_ratio = n_u/n_t;
            sigma = sigma + sigma*(acceptance_ratio/obj.target_acceptance_ratio-1)*(t/50)^(-0.5); % it will decrease adapting till iteration
            sigma = max(sigma,obj.sigma_min);
        end
        
        function par_new = generate_new_params(obj, sigma, par, par_lb, par_ub, par_shared) %#ok<INUSL>
            par_fixed_shared = (par_lb==par_ub) | (par_shared>0);% do not sample fixed or shared parameters
            cov_mat = diag((sigma*par(~par_fixed_shared)).^2); % State-Dependent Proposal Scalings
            par_new = par;
            par_new(~par_fixed_shared) = mvnrnd(par(~par_fixed_shared), cov_mat);
            par_new(par_shared>0) = par_new(par_shared(par_shared>0));
            while any(par_new(~par_fixed_shared) < par_lb(~par_fixed_shared))...
                    || any(par_new(~par_fixed_shared) > par_ub(~par_fixed_shared))
                % resample if out of bounds
                par_new(~par_fixed_shared) = mvnrnd(par(~par_fixed_shared), cov_mat);
                par_new(par_shared>0) = par_new(par_shared(par_shared>0));
            end
        end
    end
    
    methods % error function
        function error = sum_of_squares(obj, ydata_est, ydata_target, data_weights)
            % ydatas are LxN matrices (columns are data points)
            % distance between columns is norm2, sum of squares is between
            % both data sets
            if nargin<4; data_weights = ones(1,size(ydata_est,2)); end
            dist = zeros(size(data_weights));
            for i=1:size(ydata_est,2)
                % exclude nan values in ydata_target (due to combination of
                % datasets of different size)
                dist(:,i) = norm(ydata_est(~isnan(ydata_target(:,i)),i)-ydata_target(~isnan(ydata_target(:,i)),i));
            end
            error = sum(data_weights.*(dist.^2));
        end
    end
    
    methods % plot
        function [fig, ax1, pp1, pp1s, ax2, pp2, pp2s, ax3, pp3] = progress_plot(obj, ss_fits, xdata, ydata, fun, par_init, fun_args, save_progress_png)
            if exist('vecnorm','file')>0; fn_vecnorm = @vecnorm; else; fn_vecnorm = @external_tools.VecNorm; end
            fig = figure('position',[230,380,920,430],'color','w');
            % Error function
            ax1 = subplot(4,2,1:2:7,'Parent',fig); hold(ax1,'on'); box(ax1,'on');
            pp1 = plot(ax1, 0, ss_fits, 's-');
            pp1s = plot(ax1, nan, nan, ':', 'linewidth', 2);
            xlabel(ax1,'Iteration number'); ylabel(ax1,'Sum-of-squares');
            % Data and fit
            ax2 = subplot(4,2,4:2:8,'Parent',fig); hold(ax2,'on'); box(ax2,'on');
            xfine = xdata;%linspace(min(xdata)-0.1*(max(xdata)-min(xdata)),max(xdata)+0.1*(max(xdata)-min(xdata)),1000);
            pp2 = plot(ax2, nan, nan, '-', 'linewidth', 2);
            pp2s = [];
            % Residuals
            ax3 = subplot(4,2,2,'Parent',fig); hold(ax3,'on'); box(ax3,'on');
            pp3 = plot(ax3, nan, nan, 'o', 'linewidth', 1);
            ylabel(ax3,'Residuals');
            set(ax3, 'children',flipud(get(ax3,'children')));
            set(ax3,'xticklabel',[]);
            set(ax3,'fontsize',get(ax2,'fontsize'));
            
            if size(ydata, 1)==1 && size(xdata, 1)==1 % if 1D-1D case
                plot(ax2, xdata, ydata, 'o');
                pp2.XData = xfine; pp2.YData = fun(par_init, xfine, fun_args{:});
                for i=1:obj.nr_random_walks
                    pp2s(i).XData = repmat([xfine,nan],1,obj.n_iter); pp2s(i).YData = nan*ones(1,obj.n_iter*(length(xfine)+1));
                end
                xlabel(ax2,'X-Data'); ylabel(ax2,'Y-Data');
                ax2_ch = get(ax2,'children');
                set(ax2,'children',ax2_ch([numel(ax2_ch)-1,numel(ax2_ch),1:(numel(ax2_ch)-2)]));
                plot(ax3,[xfine(1),xfine(end)],[0,0],'--');
                pp3.XData = xdata; pp3.YData = ydata-fun(par_init, xdata, fun_args{:});
                xlim(ax3,[xfine(1),xfine(end)]);
                linkaxes([ax2,ax3],'x');
            else
                if size(ydata, 2)==1
                    pp2.XData = 1:size(ydata,1); pp2.YData = fun(par_init, xdata, fun_args{:});
                    xlabel(ax2,'Dimension'); ylabel(ax2,'Y-Data');
                end
                plot(ax3,[.5,size(ydata,2)+.5],[0,0],'--');
                xlim(ax3,[.5,size(ydata,2)+.5]);
                pp3.XData = 1:size(ydata,2); pp3.YData = fn_vecnorm(ydata-fun(par_init, xdata, fun_args{:}),2,1);
            end
            
            drawnow;
            if save_progress_png; export_fig('data\temp\mcmc_0','-nocrop'); end
        end
        
        function [pp1, pp1s, pp2, pp2s, pp3] = progress_plot_update(obj, ss_fits, results_fit, pp1, pp1s, pp2, pp2s, pp3, i, i_burn_in, xdata, ydata, fun, par_fit, fun_args, save_progress_png)
            if exist('vecnorm','file')>0; fn_vecnorm = @vecnorm; else; fn_vecnorm = @external_tools.VecNorm; end
            % Error function
            pp1.XData = 0:(length(ss_fits)-1);...
            pp1.YData = ss_fits;
            pp1.Parent.XLim = [max(0,length(ss_fits)-obj.k_burn_in-1), length(ss_fits)];
            if length(ss_fits)-obj.k_burn_in-1 >= i_burn_in
                % convergence annotation
                pp1s.XData = nan; pp1s.YData = nan;
            end
            % Data and fit
            if false && ((size(xdata,1)==1 && size(ydata,1)==1) || (size(ydata, 2)==1)) % 1D-1D case or N=1 case
                %if false
                xfine = xdata;
                for j=1:obj.nr_random_walks
                    if iscell(results_fit)
                        ydata_fun = results_fit{j}.ydata_fun;
                    else
                        ydata_fun = results_fit.ydata_fun;
                    end
                    yf = unique(squeeze(ydata_fun),'rows');
                    if size(ydata, 1)==1
                        pp2s(j).XData = repmat([xfine,nan],1,size(yf,1));
                    end
                    yf(:,end+1) = nan; %#ok<AGROW>                    
                    pp2s(j).YData = reshape(yf',1,size(yf,1)*size(yf,2));
                end
                pp2.YData = fun(par_fit, xfine, fun_args{:});
                %end
                % Residuals
                if size(ydata, 2)==1
                    pp3.YData = fn_vecnorm(ydata-fun(par_fit, xdata, fun_args{:}),2,1);
                else
                    pp3.YData = ydata-fun(par_fit, xdata, fun_args{:});
                end
            else
                pp3.YData = fn_vecnorm(ydata-fun(par_fit, xdata, fun_args{:}),2,1);
            end
            drawnow;
            if save_progress_png; export_fig(strcat('data\temp\mcmc_',num2str(i)),'-nocrop'); end
        end
    end
        
    methods(Static)
        function mcmc_run_check(params_mat,ss_mat,sigs,burn)
            if nargin<4; burn = 1; end
            metropolis_hastings.check_ss_trend(ss_mat);
            if nargin>2; metropolis_hastings.check_sig_trend(sigs); end
            metropolis_hastings.check_acceptance_ratio(params_mat);
            metropolis_hastings.compare_params(params_mat, {'a','b','c','a1','a2','k5'}, ss_mat, burn);
        end
        
        function check_sig_trend(sigs,range_min,range_max)
            if nargin==1
                range_min = 1;
                range_max = length(sigs);
            end
            figure;hold on;
            colors = {'r','g','b','m'};
            for i=1:size(sigs,1)
                plot(range_min:range_max,sigs(i,range_min:range_max),'-','color',colors{i});
            end
            xlabel('Iteration');
            ylabel('Sigma');
        end
        
        function check_ss_trend(ss_mat,range_min,range_max)
            if nargin==1
                range_min = 1;
                range_max = length(ss_mat);
            end
            figure;hold on;
            colors = {'r','g','b','m','c','k','y','r'};
            for i=1:size(ss_mat,1)
                plot(range_min:range_max,ss_mat(i,range_min:range_max),'-','color',colors{i});
            end
            xlabel('Iteration');
            ylabel('Sum of squares');
        end
        
        function check_acceptance_ratio(params_mat)
            x = 1:(size(params_mat,2)/50);
            colors = {'r','g','b','m'};
            figure;hold on;
            for k=1:4
                asd = squeeze(params_mat(k,:,:));
                y = zeros(size(x));
                if length(size(params_mat))==3
                    for i=x-1; y(i+1) = length(unique(asd((50*i+1):(50+50*i),:),'rows'))/50; end
                else
                    for i=x-1; y(i+1) = length(unique(asd((50*i+1):(50+50*i))))/50; end
                end
                plot(x,y,'-','color',colors{k});
            end
        end
        
        function [params_min, ss_min] = compare_params(results_fit, par_lbl, burn_in_per, plot_traj) %(par_mat, par_lbl, ss_mat, beta_mat, burn)
            if nargin<3
                burn_in_per = 0.3;
            end
            if nargin<4
                plot_traj = false;
            end
            if iscell(results_fit)
                burn_in = floor(burn_in_per*length(results_fit{1}.par_mat));
            else
                burn_in = floor(burn_in_per*length(results_fit.par_mat));
            end
            params_s = [];
            if iscell(results_fit) 
                for cpu_id=1:size(results_fit,1)
                    [~, paridx] = unique(squeeze(results_fit{cpu_id}.par_mat(:,:)),'rows');
                    idx = sort(paridx(paridx>burn_in));
                    params_s{cpu_id} = squeeze(results_fit{cpu_id}.par_mat(idx,:));
                    ss_s{cpu_id} = squeeze(results_fit{cpu_id}.ss_mat(idx,:))';
                    beta_s{cpu_id} = squeeze(results_fit{cpu_id}.beta_mat(idx,:));
                    if (burn_in-floor(burn_in))~=0
                        params_s{cpu_id} = params_s{cpu_id}(ss_s{cpu_id}<=min(min(results_fit{cpu_id}.ss_mat))*burn_in,:);
                        ss_s{cpu_id} = ss_s{cpu_id}(ss_s{cpu_id}<=min(min(results_fit{cpu_id}.ss_mat))*burn_in);
                    end
                end
                params_m = vertcat(params_s{:});
                ss_m = horzcat(ss_s{:})';
                beta_m = vertcat(beta_s{:});
            else
                [~, paridx] = unique(results_fit.par_mat,'rows');
                idx = sort(paridx(paridx>burn_in));
                params_m = results_fit.par_mat(idx,:);
                ss_m = results_fit.ss_mat(idx,:);
                beta_m = results_fit.beta_mat(idx,:);
                if (burn_in-floor(burn_in))~=0
                    params_m = params_m(ss_m<=min(ss_m)*burn_in,:);
                    beta_m = beta_m(ss_m<=min(ss_m)*burn_in,:);
                    ss_m = ss_m(ss_m<=min(ss_m)*burn_in);
                end
            end
            beta = beta_m(end,:);
            [ss_min,minix] = min(ss_m);
            params_min = params_m(minix,:);
            lp = size(params_m,2);
            fig = figure();
            ax = gobjects(lp,lp);
            warning('off','MATLAB:linkaxes:RequireDataAxes');
            for i=1:lp
                for j=i:lp
                    set(0, 'CurrentFigure', fig); % no way to use subaxis with 'parent' option
                    ax(i,j) = subaxis(lp,lp,i+(j-1)*lp, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.05);
                    hold(ax(i,j),'on'); box(ax(i,j),'on');
                    xtickformat(ax(i,j),'%.2g'); ytickformat(ax(i,j),'%.2g');
                    if i==j
                        histogram(ax(i,j),params_m(:,i),50);
                        [a1,~] = histcounts(params_m(:,i));
                        plot(ax(i,j),[params_min(i) params_min(i)], [0 max(a1)],'--r','linewidth',2);
                        mm = nansum(params_m(:,i).*exp(-(ss_m)))/nansum(exp(-(ss_m)));
                        % print with 2 decimal places
                        title(ax(i,j),sprintf('wavg %s, min %s',num2str(mm,floor(log10(mm))+3),num2str(params_min(i),floor(log10(params_min(i)))+3)));
                    else
                        m = 5;
                        if iscell(params_s) 
                            for k=1:length(params_s)
                                if plot_traj
                                    plot(ax(i,j),params_s{k}(:,i),params_s{k}(:,j));
                                else
                                    scatter(ax(i,j),params_s{k}(1:m:end,i),params_s{k}(1:m:end,j),5,exp(-beta*.5*ss_s{k}(1:m:end)),'filled',...
                                        'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
                                end
                            end
                        else
                            if plot_traj
                                plot(ax(i,j),params_m(:,i),params_m(:,j));
                            else
                                scatter(ax(i,j),params_m(1:m:end,i),params_m(1:m:end,j),5,exp(-beta*.5*ss_m(1:m:end)),'filled',...
                                    'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
                            end
                        end
                        scatter(ax(i,j),params_min(i),params_min(j),100,'xk','LineWidth',1);
                        ylim(ax(i,j),[min(params_m(:,j))-0.1*abs(max(params_m(:,j))-min(params_m(:,j))),max(params_m(:,j))+0.1*abs(max(params_m(:,j))-min(params_m(:,j)))]);
                    end
                    xlim(ax(i,j),[min(params_m(:,i))-0.1*abs(max(params_m(:,i))-min(params_m(:,i))),max(params_m(:,i))+0.1*abs(max(params_m(:,i))-min(params_m(:,i)))]);
                end
                set(0, 'CurrentFigure', fig);
                sax = subaxis(lp,lp,(i-1)*lp+1);
                ylabel(sax,par_lbl{i});
                sax = subaxis(lp,lp,(lp-1)*lp+i);
                xlabel(sax,par_lbl{i});
            end
            for i=1:lp
                linkaxes(ax(i,:),'x');
                if i>=2
                    linkaxes(ax(1:(i-1),i),'y');
                end
            end
        end
        
        function [par_min, results_fit] = load_results(filename)
            res = load(filename);
            if any(strcmp(fieldnames(res),'par_new'))
                par_min = res.par_new;
                results_fit = res;
                return;
            end
            if any(strcmp(fieldnames(res.results),'converged'))
                results_fit = res.results.converged;
            else
                results_fit = res.results.burn_in{end};
            end
            if iscell(results_fit)
                par_min = [];
                ss_min = inf;
                for i=1:numel(results_fit)
                    [ss, ind] = min(results_fit{i}.ss_mat);
                    if ss<ss_min
                        ss_min = ss;
                        par_min = results_fit{i}.par_mat(ind, :);
                    end
                end
            else
                [~, ind] = min(results_fit.ss_mat);
                par_min = results_fit.par_mat(ind, :);
            end
        end
    end
end