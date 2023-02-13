classdef egfr_ptprg_model < handle
    properties
        % from WT condition
        par = struct(... 
            'g1',0.3563,'g2',0.0249,'g3',0.0107,'g4',0.012,...
            'e1',4.22e-4,'e2',0.0196,'e3',0.0145,'e4',0.3505,...
            'a1',0.3182,'a2',0.4339,'a3',0.1487,'a4',9.89e-4,...
            'b1',1.8535,'k21',0.087,'EGF_EGFRtt',0);
        par_kin = struct('EGFRt',0.8,'k1',0.01);
        labels = {};
        init_conds = [];
    end
    
    methods
        function obj = egfr_ptprg_model()
            obj.set_initial_cond();
            obj.labels = {'EGFRpt', 'PTPRGat', 'EGF\_EGFRpt'};
        end
        
        function obj = set_initial_cond(obj)
            obj.init_conds = [0.01, 0.4, 0.0];
        end

        function set_params(obj, par)
            fn = fieldnames(obj.par);
            for i=1:length(fn)
                obj.par.(fn{i}) = par(i);
            end
        end
        
        function f = f_continuation(obj, state, bif_par)
            pars = obj.par;
            if nargin < 3
                EGF_EGFRtt = state(1,:);
            elseif strcmp(bif_par,'EGF_EGFRtt')
                EGF_EGFRtt = state(1,:);
            else
                EGF_EGFRtt = obj.par.EGF_EGFRtt;
                pars.(bif_par) = state(1,:);
            end

            f = obj.eqs_model(state(2:end,:), pars, EGF_EGFRtt);
        end

        function fp = fp_continuation(obj, state, bif_par)
            pars = obj.par;
            if nargin < 3
                bif_par = '';
                EGF_EGFRtt = obj.par.EGF_EGFRtt;
            elseif strcmp(bif_par,'EGF_EGFRtt')
                EGF_EGFRtt = state(1,:);
            else
                EGF_EGFRtt = obj.par.EGF_EGFRtt;
                pars.(bif_par) = state(1,:);
            end
            EGFRpt = state(2,:);
            EGFRnpt = 1 -EGFRpt -EGF_EGFRtt;
            PTPRGat = state(3,:);
            PTPRGit = 1 -PTPRGat;
            EGF_EGFRpt = state(4,:);
            EGF_EGFRnpt = EGF_EGFRtt -EGF_EGFRpt;
            % for all terms define derivatives over the input parameter (for now) and all of the variables
            dEGF_EGFRtt = [1, 0, 0, 0];
            dEGFRpt = [0, 1, 0, 0];
            dEGFRnpt = [-1, -1, 0, 0];
            dPTPRGat = [0, 0, 1, 0];
            dPTPRGit = [0, 0, -1, 0];
            dEGF_EGFRpt = [0, 0, 0, 1];
            dEGF_EGFRnpt = [1, 0, 0, -1];
            % use derivatives of model equations (using product rule) to perform the calculations
            fp = zeros(numel(state)-1,numel(state));
            fp(1,:) = dEGFRnpt.*(pars.e1*EGFRnpt +pars.e2*EGF_EGFRnpt +pars.a1*EGFRpt +pars.a2*EGF_EGFRpt) +EGFRnpt.*(pars.e1*dEGFRnpt +pars.e2*dEGF_EGFRnpt +pars.a1*dEGFRpt +pars.a2*dEGF_EGFRpt) -pars.g1*(dPTPRGat.*EGFRpt +PTPRGat*dEGFRpt) -pars.g2*dEGFRpt;
            fp(2,:) = dPTPRGit -pars.k21*dPTPRGat -pars.b1.*(dEGFRpt +dEGF_EGFRpt).*PTPRGat -pars.b1*(EGFRpt +EGF_EGFRpt).*dPTPRGat;
            fp(3,:) = dEGF_EGFRnpt*(pars.e3*EGFRnpt +pars.e4*EGF_EGFRnpt +pars.a3*EGFRpt +pars.a4*EGF_EGFRpt) +EGF_EGFRnpt*(pars.e3*dEGFRnpt +pars.e4*dEGF_EGFRnpt +pars.a3*dEGFRpt +pars.a4*dEGF_EGFRpt) -pars.g3*(dPTPRGat*EGF_EGFRpt +PTPRGat*dEGF_EGFRpt) -pars.g4*dEGF_EGFRpt;

            if strcmp(bif_par,'') % for jacobian matrix calculation exclude the first column
                fp(:,1) = [];
            elseif ~strcmp(bif_par,'EGF_EGFRtt') % update first column with derivates over a different parameter
                % order: g1, g2, g3, g4, e1, e2, e3, e4, a1, a2, a3, a4, b1, k21
                ind_par = strcmp(fieldnames(obj.par),bif_par);
                fp_par = [-PTPRGat.*EGFRpt, -EGFRpt, 0, 0, EGFRnpt*EGFRnpt, EGFRnpt*EGF_EGFRnpt, 0, 0, EGFRnpt*EGFRpt, EGFRnpt*EGF_EGFRpt, 0, 0, 0, 0; ...
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -PTPRGat*(EGFRpt+EGF_EGFRpt), -PTPRGat; ...
                    0, 0, -PTPRGat*EGF_EGFRpt, -EGF_EGFRpt, 0, 0, EGF_EGFRnpt*EGFRnpt, EGF_EGFRnpt*EGF_EGFRnpt, 0, 0, EGF_EGFRnpt*EGFRpt, EGF_EGFRnpt*EGF_EGFRpt, 0, 0];
                fp(:,1) = fp_par(:,ind_par);
            end
        end

        function [J, stable, eigvals] = jacobian(obj, state, bif_par)
            fp = obj.fp_continuation(state, bif_par);
            J = fp(:,2:end);
            if any(isinf(J(:)))
                stable = 0;
            elseif issparse(J)
                eigvals = eig(J);
                stable = all(real(eigs(J))<=0);
            else
                eigvals = eig(J);
                stable = all(real(eigvals)<=0);
            end
        end

        function [dstate_dt] = df_model(obj, time, state, experiment)
            try
                EGF_EGFRtt = experiment.get_input(time);
            catch ex
                EGF_EGFRtt = obj.par.EGF_EGFRtt;
            end
            if isempty(EGF_EGFRtt); EGF_EGFRtt = obj.par.EGF_EGFRtt; end
            % use model equations scaled by kinetic parameters for an ODE execution
            dstate_dt = [obj.par_kin.EGFRt; obj.par_kin.k1; obj.par_kin.EGFRt].*obj.eqs_model(state, obj.par, EGF_EGFRtt);
        end

        function eqs = eqs_model(obj, state, pars, EGF_EGFRtt)
            % main model equations (without kinetic terms)
            eqs = zeros(size(state));
            EGFRpt = state(1);                        
            PTPRGat = state(2);
            EGF_EGFRpt = state(3);
            PTPRGit = 1 -PTPRGat;
            EGFRnpt = 1 -EGFRpt -EGF_EGFRtt;
            EGF_EGFRnpt = EGF_EGFRtt -EGF_EGFRpt;
            
            eqs(1) = EGFRnpt*(pars.e1*EGFRnpt +pars.e2*EGF_EGFRnpt +pars.a1*EGFRpt +pars.a2*EGF_EGFRpt) -pars.g1*PTPRGat*EGFRpt -pars.g2*EGFRpt; %EGFRpt
            eqs(2) = PTPRGit -pars.k21*PTPRGat -pars.b1*(EGFRpt +EGF_EGFRpt)*PTPRGat; %PTPRGat
            eqs(3) = EGF_EGFRnpt*(pars.e3*EGFRnpt +pars.e4*EGF_EGFRnpt +pars.a3*EGFRpt +pars.a4*EGF_EGFRpt) -pars.g3*PTPRGat*EGF_EGFRpt -pars.g4*EGF_EGFRpt; %EGF_EGFRpt
        end
    end
end