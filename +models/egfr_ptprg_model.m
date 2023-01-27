classdef egfr_ptprg_model < handle
    properties
        % from WT condition
        par = struct(... 
            'g1',0.3563,'g2',0.0249,'g3',0.0107,'g4',0.012,...
            'e1',4.22e-4,'e2',0.0196,'e3',0.0145,'e4',0.3505,...
            'a1',0.3182,'a2',0.4339,'a3',0.1487,'a4',9.89e-4,...
            'b1',1.8535,'b2',1.8535,...
            'k21',0.087,'k3',0.0,...
            'EGF_EGFRtt',0);
        par_kin = struct('EGFRt',0.8,'k1',0.01,'k4',2.25,'k5',1.0);
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
        
        function f = f_continuation(obj, x, bif_par)
            pars = obj.par;
            if nargin < 3
                EGF_EGFRtt = x(1,:);
            elseif strcmp(bif_par,'EGF_EGFRtt')
                EGF_EGFRtt = x(1,:);
            else
                EGF_EGFRtt = obj.par.EGF_EGFRtt;
                pars.(bif_par) = x(1,:);
            end
            EGFRpt = x(2,:);
            EGFRnpt = 1 -EGFRpt -EGF_EGFRtt;
            PTPRGat = x(3,:);
            PTPRGit = 1 -PTPRGat;
            EGF_EGFRpt = x(4,:);
            EGF_EGFRnpt = EGF_EGFRtt -EGF_EGFRpt;

            f = zeros(numel(x)-1,1);
            f(1,:) = EGFRnpt*(pars.e1*EGFRnpt +pars.e2*EGF_EGFRnpt +pars.a1*EGFRpt +pars.a2*EGF_EGFRpt) -pars.g1*PTPRGat*EGFRpt -pars.g2*EGFRpt; %EGFRpt
            f(2,:) = PTPRGit -pars.k21*PTPRGat -(pars.b1*EGFRpt +pars.b2*EGF_EGFRpt)*PTPRGat; %PTPRGat
            f(3,:) = EGF_EGFRnpt*(pars.e3*EGFRnpt +pars.e4*EGF_EGFRnpt +pars.a3*EGFRpt +pars.a4*EGF_EGFRpt +pars.k3) -pars.g3*PTPRGat*EGF_EGFRpt -pars.g4*EGF_EGFRpt; %EGF_EGFRpt
        end

        function fp = fp_continuation(obj, x, bif_par)
            pars = obj.par;
            if nargin < 3
                EGF_EGFRtt = x(1,:);
            elseif strcmp(bif_par,'EGF_EGFRtt')
                EGF_EGFRtt = x(1,:);
            else
                EGF_EGFRtt = obj.par.EGF_EGFRtt;
                pars.(bif_par) = x(1,:);
            end
            EGFRpt = x(2,:);
            EGFRnpt = 1 -EGFRpt -EGF_EGFRtt;
            PTPRGat = x(3,:);
            PTPRGit = 1 -PTPRGat;
            EGF_EGFRpt = x(4,:);
            EGF_EGFRnpt = EGF_EGFRtt -EGF_EGFRpt;
            dEGF_EGFRtt = [1, 0, 0, 0];
            dEGFRpt = [0, 1, 0, 0];
            dEGFRnpt = [-1, -1, 0, 0];
            dPTPRGat = [0, 0, 1, 0];
            dPTPRGit = [0, 0, -1, 0];
            dEGF_EGFRpt = [0, 0, 0, 1];
            dEGF_EGFRnpt = [1, 0, 0, -1];
            fp = zeros(numel(x)-1,numel(x));
            fp(1,:) = dEGFRnpt.*(pars.e1*EGFRnpt +pars.e2*EGF_EGFRnpt +pars.a1*EGFRpt +pars.a2*EGF_EGFRpt) +EGFRnpt.*(pars.e1*dEGFRnpt +pars.e2*dEGF_EGFRnpt +pars.a1*dEGFRpt +pars.a2*dEGF_EGFRpt) -pars.g1*(dPTPRGat.*EGFRpt +PTPRGat*dEGFRpt) -pars.g2*dEGFRpt;
            fp(2,:) = dPTPRGit -pars.k21*dPTPRGat -(pars.b1.*dEGFRpt +pars.b2*dEGF_EGFRpt).*PTPRGat -(pars.b1*EGFRpt +pars.b2*EGF_EGFRpt).*dPTPRGat;
            fp(3,:) = dEGF_EGFRnpt*(pars.e3*EGFRnpt +pars.e4*EGF_EGFRnpt +pars.a3*EGFRpt +pars.a4*EGF_EGFRpt +pars.k3) +EGF_EGFRnpt*(pars.e3*dEGFRnpt +pars.e4*dEGF_EGFRnpt +pars.a3*dEGFRpt +pars.a4*dEGF_EGFRpt) -pars.g3*(dPTPRGat*EGF_EGFRpt +PTPRGat*dEGF_EGFRpt) -pars.g4*dEGF_EGFRpt;

            if strcmp(bif_par,'g1')
                fp(:,1) = [-PTPRGat.*EGFRpt; 0; 0];
            elseif ~strcmp(bif_par,'EGF_EGFRtt')
                % g1,g2,g3,g4,e1,e2,e3,e4,a1,a2,a3,a4,b1,b2,k21,k3
                ind_par = strcmp(fieldnames(obj.par),bif_par);
                fp_par = [-PTPRGat.*EGFRpt, -EGFRpt, 0, 0, EGFRnpt*EGFRnpt, EGFRnpt*EGF_EGFRnpt, 0, 0, EGFRnpt*EGFRpt, EGFRnpt*EGF_EGFRpt, 0, 0, 0, 0, 0, 0; ...
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -PTPRGat*EGFRpt, -PTPRGat*EGF_EGFRpt, -PTPRGat, 0; ...
                    0, 0, -PTPRGat*EGF_EGFRpt, -EGF_EGFRpt, 0, 0, EGF_EGFRnpt*EGFRnpt, EGF_EGFRnpt*EGF_EGFRnpt, 0, 0, EGF_EGFRnpt*EGFRpt, EGF_EGFRnpt*EGF_EGFRpt, 0, 0, 0, EGF_EGFRnpt];
                fp(:,1) = fp_par(:,ind_par);
            end
        end

        function [dydt] = df_model(obj, t, y, experiment)
            try
                EGF_EGFRtt = experiment.get_input(t);
            catch ex
                EGF_EGFRtt = obj.par.EGF_EGFRtt;
            end
            if isempty(EGF_EGFRtt); EGF_EGFRtt = obj.par.EGF_EGFRtt; end
            
            dydt = zeros(size(y));
            EGFRpt = y(1);                        
            PTPRGat = y(2);
            EGF_EGFRpt = y(3);
            PTPRGit = 1 -PTPRGat;
            EGFRnpt = 1 -EGFRpt -EGF_EGFRtt;
            EGF_EGFRnpt = EGF_EGFRtt -EGF_EGFRpt;
            
            dydt(1) = obj.par_kin.EGFRt*(EGFRnpt*(obj.par.e1*EGFRnpt +obj.par.e2*EGF_EGFRnpt +obj.par.a1*EGFRpt +obj.par.a2*EGF_EGFRpt) -obj.par.g1*PTPRGat*EGFRpt -obj.par.g2*EGFRpt); %EGFRpt
            dydt(2) = obj.par_kin.k1*(PTPRGit -obj.par.k21*PTPRGat -(obj.par.b1*EGFRpt +obj.par.b2*EGF_EGFRpt)*PTPRGat); %PTPRGat
            dydt(3) = obj.par_kin.EGFRt*(EGF_EGFRnpt*(obj.par.e3*EGFRnpt +obj.par.e4*EGF_EGFRnpt +obj.par.a3*EGFRpt +obj.par.a4*EGF_EGFRpt) +obj.par.k3*EGF_EGFRnpt -obj.par.g3*PTPRGat*EGF_EGFRpt -obj.par.g4*EGF_EGFRpt); %EGF_EGFRpt
        end
    end
end