
% This function defines the initial values ​​of the internal state variables,
% ie, c_s, c_e, φ_s, φ_e, and j, when t = 0 and i = 0.

% Inputs:
%    - p is a MATLAB struct containing the initial parameters
%      and their corresponding vector form.
%    - soc_init is the initial state of charge (soc) of the full cell, [-].

% Outputs:
%    - cs is initial solid electrolyte, [mol/m^3].
%    - ce is initial electrolyte electrolyte, [mol/m^3].
%    - phis is initial solid potential, [V].
%    - phie is initial electrolyte potential, [V].
%    - j is initial pore wall flux, [mol/m^2/s].

function [cs,ce,phis,phie,j] = variable_init(p,soc_init)

% initial concentration of lithium in the porous electrode [mol/m^3]
cs0_neg = p.cs_max_neg * (p.s0_neg + soc_init * (p.s100_neg - p.s0_neg));  % Equation (15-1)
cs0_pos = p.cs_max_pos * (p.s0_pos + soc_init * (p.s100_pos - p.s0_pos));  % Equation (15-2)

% solid phase Li concentration in the porous electrode: Equation (46-1)
cs(1:p.nn*p.nrn,1)       = cs0_neg;
cs(p.nn*p.nrn+1:p.nrt,1) = cs0_pos;

% electrolyte Li+ concentration in the porous electrode: Equation (46-2)
ce(:,1) = p.ce0 * ones(p.nx,1);

% solid phase potential in the porous electrode: Equations (46-3)
w = cs(p.nrn:p.nrn:p.nrn*p.nn)/p.cs_max_neg;
z = cs(p.nrn*p.nn+p.nrp:p.nrp:p.nrt)/p.cs_max_pos;
phis = [p.U_neg(w); p.U_pos(z)];

% electrolyte phase potential in the porous electrode: Equation (46-4)
phie = zeros(p.nx,1);

% reaction flux between the electrolyte and the porous electrode particle: Equation (46-5)
j = zeros(p.nnp,1);

end
