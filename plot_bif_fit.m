clear;

plot_data();

function plot_data(varargin)
    p = inputParser;
    addOptional(p,'plot_par_table',true);
    addOptional(p,'export_data',false);
    parse(p,varargin{:});

    conditions = {'MCF7_WT', 'MCF7_PTPRG_KO', 'MCF7_PTPN2_KO', 'MCF7_PTPN2_rescue', 'MCF7_p22_KO', 'MCF7_WB',...
        'MCF7_PTPRG_verylow', 'MCF7_PTPRG_medium', 'MCF7_PTPRG_high'};

    % select model used for fitting
    model = models.egfr_ptprg_model;
    filename_fitting = 'model_v10_fit_v12.5.mat';
    %filename_fitting = 'temp\tpde623877_1685_4e94_84e4_79f532ddf201.mat';
    
    data_fit = data_fitting(filename_fitting, conditions, model);
    
    plot_bif_diagrams(data_fit, model, conditions, p.Results.plot_par_table);

    if p.Results.export_data
        % export continuation data to csv files
        export_data(par_master, model, conditions);
    end
end

function plot_bif_diagrams(data_fit, model, conditions, plot_par_table)
    % plot all conditions as 
    % 1) dose response profiles overlaid with the data, and residuals
    % 2) bifurcation diagrams for G1 parameter, in absence of stimulus
    % 3) 3-d bifurcation diagrams (G1 & stimulus) with DR profile for the estimated parameter values
    par_lgth = numel(data_fit.par_master)/numel(conditions); % calculate no. parameters per condition
    for j = 1:numel(conditions)
        pars = data_fit.par_master((1:par_lgth) +(j-1)*par_lgth);
        par_min_master = data_fit.par_min_master((1:par_lgth) +(j-1)*par_lgth);
        par_max_master = data_fit.par_max_master((1:par_lgth) +(j-1)*par_lgth);
        plot_bif_diagrams_condition(pars, data_fit.xx_master(:,j), data_fit.yy_master(:,j), model, conditions{j}, [0,1], plot_par_table, par_min_master, par_max_master);
        %break;
    end
end

function fig = plot_bif_diagrams_condition(par_fit, xx, yy, model, condition, range, plot_par_table, par_min_master, par_max_master)
    if nargin<6; range = [0,1]; end
    if nargin<7; plot_par_table = false; end
    % set up figure with one top residuals plot and three bottom bifurcation diagrams
    uts = get(0,'units'); set(0,'units','pixels');
    ss = get(0,'screensize');
    set(0,'units',uts);
    offset = 100;
    fig = figure('position',[offset,offset,ss(3)-2*offset,ss(4)-2*offset],'color','w');
    n_rows = 4; n_cols = 3;
    span_ax = @(col) (n_cols+col):n_cols:(n_cols+col+(n_rows-2)*n_cols);
    ax1 = subplot(n_rows,n_cols,1,'Parent',fig); hold(ax1,'on'); box(ax1,'on');
    title(ax1,strrep(condition,'_','\_'));
    ax2 = subplot(n_rows,n_cols,span_ax(1),'Parent',fig); hold(ax2,'on'); box(ax2,'on');
    ax3 = subplot(n_rows,n_cols,span_ax(2),'Parent',fig); hold(ax3,'on'); box(ax3,'on');
    ax4 = subplot(n_rows,n_cols,span_ax(3),'Parent',fig); hold(ax4,'on'); 
    
    plot_dr(par_fit, model, xx, yy, range, ax1, ax2);
    
    xmax = plot_bif(par_fit, model, ax2, ax3);

    plot_3d_bif(par_fit, model, ax4, 3*xmax);

    ann = [];
    if plot_par_table
        txt = {'Parameters:'};
        fn = fieldnames(model.par)';
        k = 2;
        for i=1:length(fn)
            if par_min_master(i) < par_max_master(i)
                txt{k} = [strrep(fn{i},'_','\_'),' = ',num2str(par_fit(i),3)];
                k = k+1;
            end
        end
        ann = annotation(fig, 'textbox', [0.5, 0.62, 0.11, 0.37], 'String', txt, 'FontSize',10, 'Color',.2*[1,1,1], 'FitBoxToText','on');
        set_pos_ann(ann, ax3);
    end

    set(fig, 'SizeChangedFcn', {@fig_size_change, [ax1,ax2,ax3,ann]});
end

function plot_dr(par_fit, model, xx, yy, range, ax1, ax2)
    % run and plot dose-response continuation (bifurcation over the liganded fraction EGF_EGFRt)
    x = range(1):.0001:range(2);
    y = dose_response_est(par_fit, x, 'model', model);
    yy_est = dose_response_est(par_fit, xx, 'model', model);
    % plot residuals in the first (top left) axes
    plot(ax1, xx, yy-yy_est, 'o');
    plot(ax1, [min([xx',range]),max([xx',range])], [0,0], '--');
    % plot DR profile (fit) in the second (bottom left) axes
    plot(ax2, xx, yy, 'o');
    p = plot(ax2, x, y, '-', 'linewidth', 2);
    if ~isequal(range,[0,1])
        x = 0:.001:1;
        y = dose_response_est(par_fit, x, 'model', model);
        plot(ax2, x, y, '--', 'linewidth', 2, 'color',p.Color);
    end
    linkaxes([ax1,ax2],'x');
    axis(ax2, 'equal');
    set(ax2, 'FontSize', 15);
    xlim(ax2,[0,1]); ylim(ax2,[min(0,min(yy)),max(1,max(yy))]);
    pos1 = external_tools.plotboxpos(ax1);
    pos2 = external_tools.plotboxpos(ax2);
    set(ax1,'position',[pos2(1),pos1(2),pos2(3),pos1(4)]); % align position of axes with residuals
    set(ax1,'fontsize',get(ax2,'fontsize'));
    xticklabels(ax1,{});
    ylabel(ax1, 'Residuals');
    ylabel(ax2, 'Response');
    xlabel(ax2, 'Dose');
    drawnow;
end

function xmax = plot_bif(par_fit, model, ax2, ax3)
    % run and plot bifurcation diagram over G1 parameter
    model.set_params(par_fit);
    model.par.EGF_EGFRtt = 0;
    ct = dynamics.continuation(model);
    ct.calc_profile('g1', [0, 10], [0, 1]);
    g1 = ct.prof_bif.vars_cont(1,:);
    Epp = ct.readout(ct.prof_bif.vars_cont);
    p1 = plot(ax3, g1, Epp, '-', 'linewidth', 2);
    ind = find(g1<=model.par.g1,1,'last');
    g1_fr = (model.par.g1-g1(ind))/(g1(ind+1)-g1(ind));
    Epp_ind = Epp(ind)*(1-g1_fr) +Epp(ind+1)*g1_fr;
    p2 = plot(ax3, model.par.g1, Epp_ind, 'o');
    p2.MarkerFaceColor = p2.Color;
    p1.Color = 0.7*[1,1,1];
    linkaxes([ax2,ax3],'y');
    xmax = max([model.par.g1, max(g1(Epp>0.03))*2]);
    xlim(ax3,[0,xmax]);
    pos2 = external_tools.plotboxpos(ax2);
    pos3 = external_tools.plotboxpos(ax3);
    set(ax3,'position',[pos3(1),pos2(2),pos3(3),pos2(4)]);
    set(ax3,'fontsize',get(ax2,'fontsize'));
    ylabel(ax3, 'Response');
    xlabel(ax3, '\Gamma_1');
    drawnow;
end

function plot_3d_bif(par_fit, model, ax, g1_max)
    if nargin<3
        fig = figure('position',[200,100,670,720],'color','w');
        ax = axes('parent',fig);
    end
    if nargin<4
        % this value is used in the paper figures, otherwise it's estimated from the g1-bif.-diag. for speed-up
        g1_max = 30;
    end
    hold(ax,'on');
    fn = fieldnames(model.par);
    for i=1:length(fn)
        model.par.(fn{i}) = par_fit(i);
    end
    EEts = linspace(0,0.5,20);
    ct = dynamics.continuation(model);
    for i=1:numel(EEts)
        EEt = EEts(i);
        model.par.EGF_EGFRtt = EEt;
        ct.calc_profile('g1', [0, g1_max], [0, 1]);
        g1 = ct.prof_bif.vars_cont(1,:);
        Epp = ct.readout(ct.prof_bif.vars_cont);
        if i==1
            xmax = max([model.par.g1, max(g1(Epp>0.03))*2]);
            n_samples = 1000; %round(g1_max/xmax*pts_visible); % set n_samples to stretch to g1_max with frequency so as to match the number of visible pts within xmax
            g1_interp = nan*zeros(numel(EEts),n_samples);
            Ep_interp = nan*zeros(numel(EEts),n_samples);
            plot3(ax, g1, EEt*ones(size(g1)), Epp, '-', 'linewidth', 2, 'Color',0.7*[1,1,1]);
        end
        [g1_interp(i,:), Ep_interp(i,:)] = get_curve_mean_dist_points(g1,Epp,n_samples);
    end

    step = 15;
    samples = 1:step:size(g1_interp,2);
    if samples(end)~=size(g1_interp,2); samples = [samples,size(g1_interp,2)]; end
    mesh(ax,g1_interp(:,samples),repmat(EEts',1,numel(samples)),Ep_interp(:,samples),'EdgeColor',0.9*[1,1,1],'FaceAlpha','0.0');
    surf(ax,g1_interp,repmat(EEts',1,n_samples),Ep_interp,'FaceAlpha',0.4,'EdgeColor','none','FaceLighting','phong');
    if exist('clim','file')>0; clim(ax,[-1 1]); else; caxis(ax,[-1 1]); end %#ok<CAXIS> 
    EEts = 0:.001:0.5;
    Ep_dr = dose_response_est(par_fit, EEts, 'model', model);
    plot3(ax,repmat(model.par.g1,1,numel(EEts)),EEts,Ep_dr,'-','Color',[216,82,24]./255,'linewidth',2);
    xlim(ax,[0,xmax]); ylim(ax,[0,EEts(end)]); zlim(ax,[0,1]);
    xlabel(ax,'\Gamma_1'); ylabel(ax,'\alpha_L'); zlabel(ax,'\alpha_P');
    set(ax, 'FontSize', 20);
    view(ax,-35,30);
    rotate3d(ax,'on');
    drawnow;
end

function y_est = dose_response_est(par, x_data, varargin)
    p = inputParser;
    addParameter(p,'model',[]);
    parse(p,varargin{:});
    model = p.Results.model;
    if isempty(model)
        model = models.egfr_ptprg_model;
    end
    model.set_params(par);
    [xx, ~, xx_ind] = unique(x_data);
    ct = dynamics.continuation(model);
    ct.calc_profile('EGF_EGFRtt', [0, 1], [0, 1]);
    x = ct.prof_bif.vars_cont(1,:);
    y = ct.readout(ct.prof_bif.vars_cont);
    if ~isempty(ct.prof_bif.LP_inxs)
        ind_1 = min(ct.prof_bif.LP_inxs);
        if length(ct.prof_bif.LP_inxs) == 2
            ind_2 = max(ct.prof_bif.LP_inxs);
        else
            ind_2 = find(isnan(x),1,'last')+1;
        end
        ind_21 = find(x(ind_2:end)>x(ind_1),1,'first') + ind_2 -1;
        y_21 = interp1(x(ind_2:end), y(ind_2:end), x(ind_1)+1e-10, 'PCHIP', nan);
        yy = interp1([x(1:ind_1),x(ind_1)+1e-10,x(ind_21:end)], [y(1:ind_1),y_21,y(ind_21:end)], xx, 'PCHIP', nan);
    else
        yy = interp1(x, y, xx, 'PCHIP', nan);
    end
    y_est = yy(xx_ind);
end

function [xx,yy] = get_curve_mean_dist_points(x,y,num_points)
    w_xy = [1,10];
    dists = sqrt(w_xy(1).*(x(2:end)-x(1:(end-1))).^2 +w_xy(2).*(y(2:end)-y(1:(end-1))).^2);
    distscum = cumsum(dists);
    distscum = distscum./distscum(end); % cumulative distance on profile
    distscum = [0, distscum]; % coordinates really..
    eqdist = linspace(0,1,num_points);
    xx = interp1(distscum,x,eqdist);
    yy = interp1(distscum,y,eqdist);
end

function export_data(par_master, model, conditions)
    par_lgth = numel(par_master)/numel(conditions); % calculate no. parameters per condition
    filename_dr = 'data_dr.xlsx';
    filename_cont = 'data_cont.xlsx';
    for j = 1:numel(conditions)
        pars = par_master((1:par_lgth) +(j-1)*par_lgth);

        fn = fieldnames(model.par);
        for i=1:length(fn)
            model.par.(fn{i}) = pars(i);
        end

        ct = dynamics.continuation(model);
        prof_bif = ct.calc_profile('EGF_EGFRtt', [0, 1], [0, 1]);
        if prof_bif.ret_status ~= 0
            disp('ret_status ~= 0');
            return;
        end
        x = prof_bif.vars_cont(1,:);
        if ~isempty(prof_bif.LP_inxs)
            ind_1 = min(prof_bif.LP_inxs); % first LP
            if length(prof_bif.LP_inxs) == 2
                ind_2 = max(prof_bif.LP_inxs); % last LP
            else
                ind_2 = find(isnan(x),1,'last')+1;
            end
            ind_21 = find(x(ind_2:end)>=x(ind_1),1,'first') + ind_2 -1; % find first point on the upper branch, right from the first LP
            vars_dr = zeros(size(prof_bif.vars_cont,1), ind_1+2+size(prof_bif.vars_cont,2)-ind_21+1);
            for k = 1:size(prof_bif.vars_cont, 1)
                vk = interp1(x(ind_2:end), prof_bif.vars_cont(k,ind_2:end), x(ind_1), 'PCHIP', nan);
                vars_dr(k,:) = [prof_bif.vars_cont(k,1:ind_1), nan, vk, prof_bif.vars_cont(k,ind_21:end)];
            end
        else
            vars_dr = prof_bif.vars_cont;
        end
        prof_bif.vars_cont((end+1):(end+4),:) = [1-prof_bif.vars_cont(1,:)-prof_bif.vars_cont(2,:); 1-prof_bif.vars_cont(3,:); prof_bif.vars_cont(1,:)-prof_bif.vars_cont(4,:); prof_bif.vars_cont(2,:)+prof_bif.vars_cont(4,:)];
        T = array2table(prof_bif.vars_cont','VariableNames',{'alpha_L','EGFRpt','PTPRGat','EGF_EGFRpt','EGFRnpt','PTPRGit','EGF_EGFRnpt','alpha_P'});
        writetable(T,filename_cont,'Sheet',conditions{j},'WriteVariableNames',true);

        vars_dr((end+1):(end+4),:) = [1-vars_dr(1,:)-vars_dr(2,:); 1-vars_dr(3,:); vars_dr(1,:)-vars_dr(4,:); vars_dr(2,:)+vars_dr(4,:)];
        T = array2table(vars_dr','VariableNames',{'alpha_L','EGFRpt','PTPRGat','EGF_EGFRpt','EGFRnpt','PTPRGit','EGF_EGFRnpt','alpha_P'});
        writetable(T,filename_dr,'Sheet',conditions{j},'WriteVariableNames',true);
    end
end

function fig_size_change(source, event, axn) %#ok<INUSL,INUSD> 
    ax1 = axn(1); ax2 = axn(2); ax3 = axn(3); ann = axn(4);
    pos1 = external_tools.plotboxpos(ax1);
    pos2 = external_tools.plotboxpos(ax2);
    pos3 = external_tools.plotboxpos(ax3);
    set(ax1,'position',[pos2(1),pos1(2),pos2(3),pos1(4)]); % align position of axes with residuals
    set(ax3,'position',[pos3(1),pos2(2),pos3(3),pos2(4)]);
    if ~isempty(ann)
        set_pos_ann(ann, ax3);
    end
end

function set_pos_ann(ann, ax3)
    an_pos_target = sum(ax3.Position([2,4]))+0.01;
    drawnow;
    count = 0;
    eta = 0.25;
    while abs(an_pos_target-ann.Position(2))>0.01 && count<5
        ann.FontSize = ann.FontSize*(1-eta+eta*(0.99-an_pos_target)/ann.Position(4));
        drawnow;
        ann.Position(1) = ax3.Position(1)+ax3.Position(3)/2-ann.Position(3)/2;
        count = count + 1;
    end
end
