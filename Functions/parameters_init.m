
% This function defines the parameters used in the simulation.

% Inputs:
%    - no parameters need to be passed in.

% Outputs:
%    - p is a MATLAB struct containing the initial parameters.

% Parameter:
% [1] Bai et al., Decouple charge transfer reactions in the Li-ion battery,
% Journal of Energy Chemistry, 92 759-798 (2024)
% [2] Torchio et al., LIONSIMBA: a Matlab framework based on a finite volume
% model suitable for li-ion battery design, simulation, and control,
% Journal of The Electrochemical Society, 163 (7) A1192-A1205 (2016)

function p = parameters_init

% α_a and α_c: transfer coefficient of surface reaction [-]
p.alpha_a_neg = 0.5;               % anodic charge transfer coefficent in the negative electrode [-]
p.alpha_c_neg = 1 - p.alpha_a_neg; % cathodic charge transfer coefficent in the negative electrode [-]
p.alpha_a_pos = 0.5;               % anodic charge transfer coefficent in the positive electrode [-]
p.alpha_c_pos = 1 - p.alpha_a_pos; % cathodic charge transfer coefficent in the positive electrode [-]
% Note that the cathodic transfer coefficient alpha_c is automatically computed from alpha_a+alpha_c=1.

% σ: solid-phase conductivity [S/m]
p.sigma_neg = 100;
p.sigma_pos = 100;

% ε_e: porosity (electrolyte volume fraction) [-]
p.epse_neg = 0.485;
p.epse_pos = 0.385;
p.epse_sep = 0.724;

% ε_f: fill volume fraction [-]
p.epsf_neg = 0.0326;
p.epsf_pos = 0.0250;

% Bruggeman constant [-]
p.brug_neg = 1.5;
p.brug_pos = 1.5;
p.brug_sep = 1.5;

% solid-phase Li diffusion coefficient [m^2/s]
p.Ds_neg = 1.2e-14;
p.Ds_pos = 1.0e-14;

% particle radius [m]
p.rs_neg = 2.0e-6;
p.rs_pos = 2.0e-6;

% length of region of cell [m]
p.L_neg = 88e-6;
p.L_pos = 80e-6;
p.L_sep = 25e-6;

% rate constant for the electrochemical reaction [mol / (m^2·s) / (mol/m^3)^(1+p.alpha_a)]
p.k_neg = 5.031e-11;
p.k_pos = 2.334e-11;

% double-layer capacitance [F/m^2]
p.Cdl_neg = 0.1;
p.Cdl_pos = 0.1;

p.delta_sei_neg = p.rs_neg/50;   % δ_sei: sei thickness of the negative electrode [m]
p.rou_sei_neg = 1.4025e5;        % ρ_sei: sei resistivity of the negative electrode [Ω·m]
p.epse_sei_neg = 3.9216e-10;     % ε_sei: sei permittivity of the negative electrode [F/m]

%% Electrolyte parameter
p.ce0 = 1000;           % initial bulk electrolyte concentration [mol/m^3]

% The liquid/salt/polymer system consisting of a 10:27:63 v/v mixture of
% propylene carbonate/ethylene carbonate/dimethyl carbonate (PC/EC/DMC)
% is chosen as the electrolyte.

% Li+ transference number [-]
p.t_plus = 0.38;

% Conductivity of Li+ in the solution phase [S/m]
p.kappa = @(c,T) 1e-4*c.*((-10.5+0.668*1e-3*c+0.494*1e-6*c.^2)...
                         + (0.074-1.78*1e-5*c-8.86*1e-10*c.^2).*T...
                        + (-6.96*1e-5+2.8*1e-8*c).*T.^2).^2;

% Diffusion coefficient of the salt in the electrolyte phase [m^2/s]
p.De = @(c,T) 1e-4*10.^((-4.43-54./(T-229-5e-3*c)-0.22e-3*c));

% An electrolyte activity coefficient term of dlnf±/dlnce [-]
p.dlnfdce = @(c,T) (0.601-0.24*(c/1000).^0.5+0.982.*(1-0.0052*(T-294))*(c/1000).^1.5)*(1-p.t_plus)^-1-1;

%% Thermodynamic data
% θ: stoichiometries [-]
p.s0_neg   = 0.01429;   % at   0% SoC in the negative electrode of a fresh cell [-]
p.s100_neg = 0.85510;   % at 100% SoC in the negative electrode of a fresh cell [-]
p.s0_pos   = 0.99174;   % at   0% SoC in the positive electrode of a fresh cell [-]
p.s100_pos = 0.49950;   % at 100% SoC in the positive electrode of a fresh cell [-]

% the open-circuit potential function [V]
p.U_neg = @(w) 0.7222 + 0.1387*w + 0.029*w.^0.5 - 0.0172./w...
             + 0.0019./w.^1.5 + 0.2808*exp(0.9-15*w) - 0.7984*exp(0.4465*w - 0.4108);
p.U_pos = @(z) (- 4.656+88.669*z.^2 - 401.119*z.^4 ...
                + 342.909*z.^6 - 462.471*z.^8 + 433.434*z.^10)./...
               (- 1+18.933*z.^2 - 79.532*z.^4 ...
                + 37.311*z.^6 - 73.083*z.^8 + 95.96*z.^10);

% maximum solid-phase concentration [mol/m^3]
p.cs_max_neg = 30555;
p.cs_max_pos = 51554;

%% Constant 
p.F = 96487;            % Faraday's constant [C/mol]
p.R = 8.314;            % ideal gas constant [J/mol/K]
p.T = 298.15;           % reference temperature [K]

end
