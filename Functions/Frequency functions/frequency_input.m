
% This function selects the corresponding concentration and calculates the EIS
% based on charge and discharge rate.

% Inputs:
%    - f is frequency, Hz.
%    - p is a MATLAB struct containing the initial parameters.
%    - volt is a MATLAB struct containing the output variables,
%      such as the output voltage.
%    - vari is a MATLAB struct containing the internal variables,
%      such as the concentrations and the potentials.
%    - Crates is the charge and discharge rate, dimensionless.
%    - soc_thr is the soc selected for calculating SEIS/DEIS, dimensionless.
%    - mode is the selection switch for charging or discharging mode.

% Outputs:
%    - state is a MATLAB struct containing the output variables,
%      such as the output impedance and internal state variables.

function state = frequency_input(f,p,volt,vari,Crates,soc_thr,mode)

tic()

p = fcn_system_vectors(p);
m = fcn_system_matrices(p);

if Crates == 0   % SEIS under the conditions of rate = 0 and a given SOC
    for n_thr = 1:length(soc_thr)
        % calculate static variables
        cs_soc0(1:p.nn*p.nrn,n_thr)       = p.cs_max_neg...
            * (p.NP.s0_neg + soc_thr(n_thr) * (p.NP.s100_neg - p.NP.s0_neg));  % Equation (54-2)
        cs_soc0(p.nn*p.nrn+1:p.nrt,n_thr) = p.cs_max_pos...
            * (p.NP.s0_pos + soc_thr(n_thr) * (p.NP.s100_pos - p.NP.s0_pos));
        % calculate static EIS
        [vari_soc{n_thr},Z_soc{n_thr}] = DFN_frequency(f,p,m,cs_soc0(:,n_thr),p.ce0 * ones(p.nx,1));
        % calculate static variables
        vari_soc{n_thr}.phis = [p.U_neg(p.NP.s0_neg + soc_thr(n_thr) * (p.NP.s100_neg - p.NP.s0_neg)) * ones(p.nn,1);...
                                NaN(p.ns,1);...
                                p.U_pos(p.NP.s0_pos + soc_thr(n_thr) * (p.NP.s100_pos - p.NP.s0_pos)) * ones(p.np,1)];
        vari_soc{n_thr}.phie = zeros(p.nx,1);
        vari_soc{n_thr}.j_dim = [zeros(p.nn,1); NaN(p.ns,1); zeros(p.np,1)];

        fprintf('At soc = %2.2f, finished the static EIS simulation \n',soc_thr(n_thr));
    end
else   % DEIS under the conditions of charging-discharging at rate ≠ 0 to the set SOC
    n_thr = 1;
    for n_t = 1:length(volt.soc)
        if mode == 1 && (volt.soc(n_t) < soc_thr(n_thr) && volt.soc(n_t) >= 0)...
        || mode == 2 && (volt.soc(n_t) > soc_thr(n_thr) && volt.soc(n_t) <= 1)
            % calculate dynamic EIS
            [vari_soc{n_thr},Z_soc{n_thr}] = DFN_frequency(f,p,m,vari.cs(:,n_t),vari.ce(:,n_t));
            % calculate dynamic variables
            vari_soc{n_thr}.phis = vari.phis(:,n_t);
            vari_soc{n_thr}.phie = vari.phie(:,n_t);
            vari_soc{n_thr}.j_dim = vari.j(:,n_t) * p.F / volt.i_app(n_t);
            
            fprintf('At soc = %2.2f, finished the dynamic EIS simulation \n',soc_thr(n_thr));
            
            n_thr = n_thr +1;
            if n_thr > length(soc_thr)
                break
            end
        end
    end
end

if Crates == 0
    dispstat(sprintf('At rate = %2.2f, finished the static EIS simulation in %2.2f s \n',Crates,toc()),'keepthis','timestamp');
else
    dispstat(sprintf('At rate = %2.2f, finished the dynamic EIS simulation in %2.2f s \n',Crates,toc()),'keepthis','timestamp');
end

%% Store states for the later plotting.
state.p = p;
state.f = f;
state.Crates  = Crates;
state.soc_thr = soc_thr;
state.vari = vari_soc;
state.Z   = Z_soc;

end
