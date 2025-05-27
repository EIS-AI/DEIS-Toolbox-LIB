
% This function computes F_x.

% Inputs:
%    - p is a MATLAB struct containing the initial parameters
%      and their corresponding vector form.
%    - m is a MATLAB struct containing the system matrix,
%      such as the discretization matrix of the governing equations.
%    - x is the dependent variable being solved for.
%    - cs_prevt is solid electrolyte of the previous moment, [mol/m^3].
%    - ce_prevt is electrolyte electrolyte of the previous moment, [mol/m^3].
%    - phis_prevt is solid potential of the previous moment, [V].
%    - phie_prevt is electrolyte potential of the previous moment, [V].
%    - phisf_prevt is potential in the film close to the solid phase
%      of the previous moment, [V].
%    - i_app is the applied (dis) charge current density, [A/m^2].

% Outputs:
%    - F_x is function matrix.
%    - cs is solid electrolyte, [mol/m^3].
%    - ce is electrolyte electrolyte, [mol/m^3].
%    - phis is solid potential, [V].
%    - phie is electrolyte potential, [V].
%    - j is pore wall flux, [mol/m^2/s].
%    - phisf is potential in the film close to the solid phase, [V].
%    - jF is molar flux of Li-ions de-intercalating out of the solid phase, [mol/m^2/s].
%    - eta is overpotential, η, [V].
%    - Uocp is open-circuit potential, [V].
%    - m is a MATLAB struct containing the system matrix,
%      such as the discretization matrix of the governing equations.

function [F_x,cs,ce,phis,phie,j,phisf,jF,eta,Uocp,m]...
        = fcn_F(p,m,x,cs_prevt,ce_prevt,phis_prevt,phie_prevt,phisf_prevt,i_app)

% state variables: Equation (47-3)
cs  = x(1:p.nrt,1);
ce   = x((p.nrt+1):...
         (p.nrt+p.nx),1);
phis = x((p.nrt+p.nx+1):...
         (p.nrt+p.nx+p.nnp),1);
phie = x((p.nrt+p.nx+p.nnp+1)...
        :(p.nrt+p.nx+p.nnp+p.nx),1);
j    = x((p.nrt+p.nx+p.nnp+p.nx+1:p.nt),1);

m = fcn_system_matrices_change(p,m,ce,p.T,phis_prevt,phie_prevt,phisf_prevt);

% the vector of solid-phase surface concentrations
cs_bar = m.Acs_bar*cs;                % Equation (39-1): c-_s
% the parts of c_e(t_k) and φ_e(t_k) given in the electrodes
ce_bar = m.Ace_bar*ce;                % Equation (39-2): c-_e
phie_bar = m.Aphie_bar*phie;          % Equation (39-3): φ-_e

% analytical expression for φ_sf, j_dl and j_F: Equations (43-1)-(43-3)
phisf = m.Aphisf_hat*phis + m.Bphisf_hat*phie_bar + m.Cphisf_hat*j + m.Dphisf_hat;
jdl   = m.Ajdl_hat  *phis + m.Bjdl_hat  *phie_bar + m.Cjdl_hat  *j + m.Djdl_hat;
jF    = m.AjF_hat   *phis + m.BjF_hat   *phie_bar + m.CjF_hat   *j + m.DjF_hat;

% exchange current density for the insertion process: Equation (38-1)
i0 = p.F*p.k .* ce_bar.^p.alpha_a .* (p.cs_bar_max-cs_bar).^p.alpha_a .* cs_bar.^p.alpha_c;

w = cs_bar(1:p.nn)      /p.cs_max_neg;
z = cs_bar(p.nn+1:p.nnp)/p.cs_max_pos;

% open-circuit potential of the porous electrode: Equations (40-4)
Uocp = [p.U_neg(w); p.U_pos(z)];

% overpotential in the porous electrode: Equations (38-4)
eta = phis - phisf - Uocp;

% Equations (38-2)-(38-3)
expI  = exp( p.alpha_a*p.F.*eta/(p.R*p.T));
expII = exp(-p.alpha_c*p.F.*eta/(p.R*p.T));

% Equations (45-1)-(45-5)
F_I_cs     = m.Acs_hat*cs + m.Bcs_hat*jF + cs_prevt;
F_II_ce    = m.Ace_hat*ce + m.Bce_hat*j  + ce_prevt;
F_III_phis = m.Aphis*phis + m.Bphis*j + m.Cphis*i_app;
F_IV_phie  = m.Aphie*phie + m.Bphie*j + m.Dphie*log(ce);
F_V_j      = p.F*jF./i0 - (expI-expII);

% the function matrix: Equations (47-2)
F_x = [F_I_cs;F_II_ce;F_III_phis;F_IV_phie;F_V_j];

end
