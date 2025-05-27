
% This function computes the static EIS and dynamic EIS of the DFN-like
% impedance models based on the input concentration variables.

% Inputs:
%    - f is frequency, Hz.
%    - p is a MATLAB struct containing the initial parameters
%      and their corresponding vector form.
%    - m is a MATLAB struct containing the system matrix,
%      such as the discretization matrix of the governing equations.
%    - cs is solid electrolyte, [mol/m^3].
%    - ce is electrolyte electrolyte, [mol/m^3].

% Outputs:
%    - vari is a MATLAB struct containing the internal variables,
%      such as the concentrations, the potentials, and charge transfer resistances.
%    - Z is a MATLAB struct containing the output variables,
%      such as the impedance of different components and diffusion impedance.

function [vari,Z] = DFN_frequency(f,p,m,cs,ce)

%% time domain variable
cs_bar = m.Acs_bar*cs;   % Equation (39-1)
ce_bar = m.Ace_bar*ce;   % Equation (39-2)

% stoichiometry of Li [-]
w_t = cs_bar(1:p.nn)       / p.cs_max_neg;
z_t = cs_bar(p.nn+1:p.nnp) / p.cs_max_pos;

% first derivative of open-circuit potential at the electrode [V·m^3/mol]
syms w z
dUdcs_bar = diag([double(subs(diff(p.U_neg,w), w, w_t) / p.cs_max_neg);...
                  double(subs(diff(p.U_pos,z), z, z_t) / p.cs_max_pos)]);     % Equation (48-17)

% exchange current density for the insertion process [A/m^2]
i0 = p.F*p.k .* ce_bar.^p.alpha_a .* (p.cs_bar_max-cs_bar).^p.alpha_a .* cs_bar.^p.alpha_c; % Equation (38-1)

% charge transfer resistance, [Ω·m^2].
DFN.Rct = sparse(diag(p.R*p.T/p.F ./ (p.alpha_a+p.alpha_c) ./ i0));  % Equation (75-4)

%% Store the calculated internal state variables for the later plotting.
vari.x = [               p.dx_n/2*(1:2:(p.nn*2-1))...
                 p.L_neg+p.dx_s/2*(1:2:(p.ns*2-1))...
         p.L_neg+p.L_sep+p.dx_p/2*(1:2:(p.np*2-1))]'...
      / (p.L_neg+p.L_sep+p.L_pos);

vari.stoich    = [w_t; NaN(p.ns,1); z_t];
vari.ce_dim    = ce / p.ce0;
vari.Rct = full([diag(DFN.Rct(1:p.nn,1:p.nn));...
                 NaN(p.ns,1);...
                 diag(DFN.Rct(p.nn+1:p.nnp,p.nn+1:p.nnp))]);

%% The system matrices that do not change in the inner loop.
m.Acs_bar_at   = m.Acs_bar;    % Equation (71-14)
m.Aphie_bar_at = m.Aphie_bar;  % ~
m.Acs_at = m.Acs;              % Equation (71-10)
m.Bcs_at = m.Bcs;              % ~
m.Bce_at = m.Bce;              % Equation (71-11)
m.Aphis_at = m.Aphis;          % Equation (71-12)
m.Bphis_at = m.Bphis;          % ~
m.Cphis_at = m.Cphis;          % ~
m.Bphie_at = m.Bphie;          % Equation (71-13)

%% The system matrices that do change in the inner loop.
change.Ace_at = sparse(fcn_Aw(p.De_eff(ce,p.T),p.eps_e.*p.dx,p.dx));       % Equations (30-4) and (71-11)

change.Aphie_at = sparse(fcn_Aw(p.kappa_eff(ce,p.T),ones(p.nx,1),p.dx));   % Equations (35-4) and (71-13)
change.Aphie_at(1,1:2) = [1 0];   % reference potential: φ_e1(x=0) = 0.   % Equation (34)

change.Dphie_at = sparse(fcn_Aw(p.nu(ce,p.T)./ce,ones(p.nx,1),p.dx));      % Equations (35-8)
change.Dphie_at(1,1:2) = [0 0];   % reference potential: φ_e1(x=0) = 0.   % Equation (34)

%% computer static electrochemical impedance spectroscopy (EIS)
for k = 1:length(f)

    s = i * 2 * pi * f(k);

%% Model DFN
% single particle model on the particle scale
    DFN.z_d    = sparse(1 / p.F * dUdcs_bar * m.Acs_bar_at / (s*eye(p.nrt)-m.Acs_at) * m.Bcs_at);  % Equation (75-5)
    DFN.z_F    = sparse(DFN.Rct + DFN.z_d);                                                        % Equation (75-1)
    DFN.z_R    = sparse(inv(s*diag(p.Cdl)  + inv(DFN.z_F)));                                       % Equation (75-2)
    DFN.z_int  = sparse(inv(s*diag(p.Csei) + inv(DFN.z_R + diag(p.Rsei))));                        % Equation (75-3)

% on the electrode scale
    DFN.dCedJ   = sparse(eye(p.nx) / (s*eye(p.nx) - change.Ace_at) * m.Bce_at);                % Equation (77-1): M_1
    DFN.dPhiedJ = sparse(-inv(change.Aphie_at) * (m.Bphie_at + change.Dphie_at * DFN.dCedJ));  % Equation (77-2): M_2
    DFN.dJdPhis = sparse(inv(m.Aphie_bar_at * DFN.dPhiedJ + p.F * DFN.z_int));                 % Equation (77-3): M_3

% internal state variables
    DFN.Phis = sparse(-inv(m.Aphis_at + m.Bphis_at * DFN.dJdPhis) * m.Cphis_at);  % Equation (76-4)
    DFN.J    = DFN.dJdPhis * DFN.Phis;                                            % Equation (76-3)
    DFN.Phie = DFN.dPhiedJ * DFN.J;                                               % Equation (76-2)
    DFN.Ce   = DFN.dCedJ * DFN.J;                                                 % Equation (76-1)

% the center of the control volume elements
    DFN = variables_full_frequency(p,ce,DFN,'DFN');

%% Model B
% single particle model on the particle scale
    B.z_int = DFN.z_int;

% on the electrode scale
    B.dCedJ   = sparse(zeros(p.nx,p.nnp));                                                 % Equation (77-1): M_1
    B.dPhiedJ = sparse(-inv(change.Aphie_at) * (m.Bphie_at + change.Dphie_at * B.dCedJ));  % Equation (77-2): M_2
    B.dJdPhis = sparse(inv(m.Aphie_bar_at * B.dPhiedJ + p.F * B.z_int));                   % Equation (77-3): M_3

% internal state variables
    B.Phis = sparse(-inv(m.Aphis_at + m.Bphis_at * B.dJdPhis) * m.Cphis_at);      % Equation (76-4)
    B.J    = B.dJdPhis * B.Phis;                                                  % Equation (76-3)
    B.Phie = B.dPhiedJ * B.J;                                                     % Equation (76-2)

% the center of the control volume elements
    B = variables_full_frequency(p,ce,B,'sim');

%% Model E
% single particle model on the particle scale
    E.z_d   = sparse(zeros(p.nnp));                                      % Equation (75-5)
    E.z_F   = sparse(DFN.Rct + E.z_d);                                   % Equation (75-1)
    E.z_R   = sparse(inv(s*diag(p.Cdl)  + inv(E.z_F)));                  % Equation (75-2)
    E.z_int = sparse(inv(s*diag(p.Csei) + inv(E.z_R + diag(p.Rsei))));   % Equation (75-3)

% on the electrode scale
    E.dCedJ   = sparse(zeros(p.nx,p.nnp));                                                 % Equation (77-1): M_1
    E.dPhiedJ = sparse(-inv(change.Aphie_at) * (m.Bphie_at + change.Dphie_at * E.dCedJ));  % Equation (77-2): M_2
    E.dJdPhis = sparse(inv(m.Aphie_bar_at * E.dPhiedJ + p.F * E.z_int));                   % Equation (77-3): M_3

% internal state variables
    E.Phis = sparse(-inv(m.Aphis_at + m.Bphis_at * E.dJdPhis) * m.Cphis_at);   % Equation (76-4)
    E.J    = E.dJdPhis * E.Phis;                                               % Equation (76-3)
    E.Phie = E.dPhiedJ * E.J;                                                  % Equation (76-2)

% the center of the control volume elements
    E = variables_full_frequency(p,ce,E,'sim');

%% impedance: Equations (78-1)-(78-4)
% Model DFN
    DFN.Z_neg(k,1)  = - (DFN.Phie_full(p.nn+1)   - DFN.Phis_full(1));
    DFN.Z_pos(k,1)  = - (DFN.Phis_full(p.nx+1)   - DFN.Phie_full(p.nns+1));
    DFN.Z_sep(k,1)  = - (DFN.Phie_full(p.nns+1)  - DFN.Phie_full(p.nn+1));
    DFN.Z_cell(k,1) = - (DFN.Phis_full(p.nx+1)   - DFN.Phis_full(1));

% Model B
    B.Z_neg(k,1)  = - (B.Phie_full(p.nn+1)   - B.Phis_full(1));
    B.Z_pos(k,1)  = - (B.Phis_full(p.nx+1)   - B.Phie_full(p.nns+1));
    B.Z_sep(k,1)  = - (B.Phie_full(p.nns+1)  - B.Phie_full(p.nn+1));
    B.Z_cell(k,1) = - (B.Phis_full(p.nx+1)   - B.Phis_full(1));

% Model E
    E.Z_neg(k,1)  = - (E.Phie_full(p.nn+1)   - E.Phis_full(1));
    E.Z_pos(k,1)  = - (E.Phis_full(p.nx+1)   - E.Phie_full(p.nns+1));
    E.Z_sep(k,1)  = - (E.Phie_full(p.nns+1)  - E.Phie_full(p.nn+1));
    E.Z_cell(k,1) = - (E.Phis_full(p.nx+1)   - E.Phis_full(1));

end

%% Store the calculated EIS for the later plotting.
% the impedance of Model DFN and Model E
Z.DFN_neg  = DFN.Z_neg;
Z.DFN_pos  = DFN.Z_pos;
Z.DFN_sep  = DFN.Z_sep;
Z.DFN_cell = DFN.Z_cell;

Z.MoE_neg  = E.Z_neg;
Z.MoE_pos  = E.Z_pos;
Z.MoE_cell = E.Z_cell;

% the impedance of sub process
if f(end) > 10^0
    for m = 1:length(f)
        if f(m) > 10^0
            M = m - 1;
            break
        end
    end
else
    M = length(f);
end

% solid diffusion impedance: Equations (80-1)-(80-2)
Z.Ds_neg  = B.Z_neg(1:M)  - E.Z_neg(1:M);
Z.Ds_pos  = B.Z_pos(1:M)  - E.Z_pos(1:M);
Z.Ds_cell = B.Z_cell(1:M) - E.Z_cell(1:M);

% electrolyte diffusion impedance: Equations (80-3)-(80-5)
Z.De_neg  = DFN.Z_neg(1:M)  - B.Z_neg(1:M);
Z.De_pos  = DFN.Z_pos(1:M)  - B.Z_pos(1:M);
Z.De_sep  = DFN.Z_sep(1:M)  - B.Z_sep(1:M);
Z.De_cell = DFN.Z_cell(1:M) - B.Z_cell(1:M);

% the sum of the charge transfer impedance at the solid/electrolyte interface Z_CT
% and the charge transfer impedance through the SEI film Z_SEI.
Z.CS_neg  = E.Z_neg(M+1:end)  - real(E.Z_neg(end));
Z.CS_pos  = E.Z_pos(M+1:end)  - real(E.Z_pos(end));
Z.CS_cell = E.Z_cell(M+1:end) - real(E.Z_cell(end));

end
