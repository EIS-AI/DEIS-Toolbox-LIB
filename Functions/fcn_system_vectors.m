
% This function converts the specified parameters into vectors.

% Inputs:
%    - p is a MATLAB struct containing the initial parameters.

% Outputs:
%    - p is a MATLAB struct containing the initial parameters
%      and their corresponding vector form.

function p = fcn_system_vectors(p)

%% number of nodes (spatial discretization)
p.nn  = p.gridsize(1);    % across negative electrode
p.ns  = p.gridsize(2);    % across seperator
p.np  = p.gridsize(3);    % across positive electrode
p.nrn = p.gridsize(4);    % in the radial direction in the negative electrode
p.nrp = p.gridsize(5);    % in the radial direction in the positive electrode

p.nnp = p.nn + p.np;               % for φ_s or j
p.nns = p.nn + p.ns;
p.nx  = p.nn + p.ns + p.np;        % for c_e or φ_e
p.nrt = p.nrn*p.nn + p.nrp*p.np;   % for c_s
p.nt  = p.nrt + 2*p.nx + 2*p.nnp;  % for c_s + c_e + φ_s + φ_e + j

% interval size
p.dx_n = p.L_neg / p.nn;        % along the x-dimension in the negative electrode
p.dx_s = p.L_sep / p.ns;        % along the x-dimension in the seperator
p.dx_p = p.L_pos / p.np;        % along the x-dimension in the positive electrode
p.dr_n = p.rs_neg / (p.nrn-1);  % along the radial direction in the negative electrode
p.dr_p = p.rs_pos / (p.nrp-1);  % along the radial direction in the positive electrode

p.dx = [p.dx_n*ones(p.nn,1);...
        p.dx_s*ones(p.ns,1);...
        p.dx_p*ones(p.np,1)];

% the range of the electrodes throughout the cell
p.elec_range = [linspace(1,p.nn,p.nn)...
                linspace(p.nns+1,p.nx,p.np)];

%% parameter
p.k          = [p.k_neg*ones(p.nn,1);...
                p.k_pos*ones(p.np,1)];       % [mol / (m^2·s) / (mol/m^3)^(1+p.alpha_a)]
p.cs_bar_max = [p.cs_max_neg*ones(p.nn,1);...
                p.cs_max_pos*ones(p.np,1)];  % [mol/m^3]
p.alpha_a    = [p.alpha_a_neg*ones(p.nn,1);...
                p.alpha_a_pos*ones(p.np,1)]; % [-]
p.alpha_c    = [p.alpha_c_neg*ones(p.nn,1);...
                p.alpha_c_pos*ones(p.np,1)]; % [-]
p.Cdl        = [p.Cdl_neg*ones(p.nn,1);...
                p.Cdl_pos*ones(p.np,1)];     % [F/m^2]

p.eps_e = [p.epse_neg*ones(p.nn,1);...
           p.epse_sep*ones(p.ns,1);...
           p.epse_pos*ones(p.np,1)];         % [-]
p.brug  = [p.brug_neg*ones(p.nn,1);...
           p.brug_sep*ones(p.ns,1);...
           p.brug_pos*ones(p.np,1)];         % [-]

% sei resistance [Ωm^2]
p.Rsei_neg = p.rou_sei_neg * (p.rs_neg * p.delta_sei_neg) / (p.rs_neg + p.delta_sei_neg);
p.Rsei_pos = 0;
p.Rsei = [p.Rsei_neg*ones(p.nn,1);...
          p.Rsei_pos*ones(p.np,1)];

% sei capacitance [F/m^2]
p.Csei_neg = p.epse_sei_neg * (p.rs_neg + p.delta_sei_neg) / (p.rs_neg * p.delta_sei_neg);
p.Csei_pos = 0;
p.Csei = [p.Csei_neg*ones(p.nn,1);...
          p.Csei_pos*ones(p.np,1)];

% ε_s: active material volume fraction at the electrodes [-]
p.epss_neg = 1 - p.epse_neg - p.epsf_neg;
p.epss_pos = 1 - p.epse_pos - p.epsf_pos;
p.eps_s = [p.epss_neg*ones(p.nn,1);...
           p.epss_pos*ones(p.np,1)];

% specific surface area of the porous electrode [m^2/m^3]
p.a_s_neg = 3*p.epss_neg/p.rs_neg;
p.a_s_pos = 3*p.epss_pos/p.rs_pos;
p.a_s = [p.a_s_neg*ones(p.nn,1);...
         p.a_s_pos*ones(p.np,1)];

% σ_eff: effective solid conductivity [S/m]
p.sigma_eff_neg = p.sigma_neg*p.epss_neg * ones(p.nn,1);
p.sigma_eff_pos = p.sigma_pos*p.epss_pos * ones(p.np,1);
p.sigma_eff = [p.sigma_eff_neg;...
               p.sigma_eff_pos];

% effective electrolyte diffusivity [m^2/s]
p.De_eff = @(c,T) p.De(c,T).*p.eps_e.^p.brug;

% κ_eff: effective electrolyte conductivity [S/m]
p.kappa_eff = @(c,T) p.kappa(c,T).*p.eps_e.^p.brug;

% κ_D,eff: effective diffusional conductivity [A/m]
p.nu = @(c,T) -2*p.R*T./p.F.*p.kappa_eff(c,T)*(1-p.t_plus).*(1+p.dlnfdce(c,T));

% first derivative of open-circuit potential at the electrode [V·m^3/mol]
npoints = 1e5;
x = linspace(0,1,npoints);
U_neg_full = p.U_neg(x);
U_pos_full = p.U_pos(x);

dU_neg = (U_neg_full(2:npoints)-U_neg_full(1:npoints-1))./(x(2:npoints)-x(1:npoints-1));
dU_pos = (U_pos_full(2:npoints)-U_pos_full(1:npoints-1))./(x(2:npoints)-x(1:npoints-1));
p.dU_neg = @(w) qinterp1(x(2:npoints)',dU_neg',w);
p.dU_pos = @(z) qinterp1(x(2:npoints)',dU_pos',z);

end
