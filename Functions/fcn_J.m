
% This function computes the Jacobian matrix J_x of F_x.

% Inputs:
%    - p is a MATLAB struct containing the initial parameters
%      and their corresponding vector form.
%    - m is a MATLAB struct containing the system matrix,
%      such as the discretization matrix of the governing equations.
%    - cs is solid electrolyte, [mol/m^3].
%    - ce is electrolyte electrolyte, [mol/m^3].
%    - jF is molar flux of Li-ions de-intercalating out of the solid phase, [mol/m^2/s].
%    - eta is overpotential, η, [V].

% Outputs:
%    - J_x is Jacobian matrix of F_x.

function J_x = fcn_J(p,m,cs,ce,jF,eta)

% the vector of solid-phase surface concentrations
cs_bar = m.Acs_bar*cs;   % Equation (39-1)
% the parts of c_e(t_k) given in the electrodes
ce_bar = m.Ace_bar*ce;   % Equation (39-2)

% exchange current density for the insertion process: Equation (38-1)
i0 = p.F*p.k .* ce_bar.^p.alpha_a .* (p.cs_bar_max-cs_bar).^p.alpha_a .* cs_bar.^p.alpha_c;

di0dcs_bar = diag((-p.alpha_a./(p.cs_bar_max-cs_bar)+p.alpha_c./cs_bar).*i0);   % Equation (48-15)
di0dce_bar = diag(p.alpha_a./ce_bar.*i0);                                       % Equation (48-16)

% first derivative of open-circuit potential at the porous electrode [V·m^3/mol]
w = cs_bar(1:p.nn)      /p.cs_max_neg;
z = cs_bar(p.nn+1:p.nnp)/p.cs_max_pos;
dUdcs_bar = diag([p.dU_neg(w)/p.cs_max_neg; p.dU_pos(z)/p.cs_max_pos]);         % Equation (48-17)

% Butler-Volmer electrochemical kinetic expression
expI  = exp( p.alpha_a*p.F.*eta/(p.R*p.T));      % Equation (38-2)
expII = exp(-p.alpha_c*p.F.*eta/(p.R*p.T));      % Equation (38-3)

dexpIdeta  =  diag(p.alpha_a*p.F/(p.R*p.T).*expI);  % based on Equation (38-2)
dexpIIdeta = -diag(p.alpha_c*p.F/(p.R*p.T).*expII); % based on Equation (38-3)

dexpIdcs_bar  = dexpIdeta * (-dUdcs_bar);        % Equation (48-7)
dexpIIdcs_bar = dexpIIdeta * (-dUdcs_bar);       % Equation (48-8)

dexpIdphis  = dexpIdeta * (eye(p.nnp) - m.Aphisf_hat);       % Equation (48-9)
dexpIIdphis = dexpIIdeta * (eye(p.nnp) - m.Aphisf_hat);      % Equation (48-10)

dexpIdphie_bar  = dexpIdeta * (-m.Bphisf_hat);     % Equation (48-11)
dexpIIdphie_bar = dexpIIdeta * (-m.Bphisf_hat);    % Equation (48-12)

dexpIdj  = dexpIdeta * (-m.Cphisf_hat);          % Equation (48-13)
dexpIIdj = dexpIIdeta * (-m.Cphisf_hat);         % Equation (48-14)

% Equations (45-1) and (48-1)
dFIdcs   = m.Acs_hat;
dFIdce   = sparse(zeros(p.nrt,p.nx));
dFIdphis = m.Bcs_hat * m.AjF_hat;
dFIdphie = m.Bcs_hat * m.BjF_hat * m.Aphie_bar;
dFIdj    = m.Bcs_hat * m.CjF_hat;

% Equations (45-2) and (48-1)
dFIIdcs   = zeros(p.nx,p.nrt);
dFIIdce   = m.Ace_hat;
dFIIdphis = zeros(p.nx,p.nnp);
dFIIdphie = zeros(p.nx,p.nx);
dFIIdj    = m.Bce_hat;

% Equations (45-3) and (48-1)
dFIIIdcs   = zeros(p.nnp,p.nrt);
dFIIIdce   = zeros(p.nnp,p.nx);
dFIIIdphis = m.Aphis;
dFIIIdphie = zeros(p.nnp,p.nx);
dFIIIdj    = m.Bphis;

% Equations (45-4) and (48-1)
dFIVdcs   = zeros(p.nx,p.nrt);
dFIVdce   = m.Dphie./ce;
dFIVdphis = zeros(p.nx,p.nnp);
dFIVdphie = m.Aphie;
dFIVdj    = m.Bphie;

% Equations (45-5) and (48-2)-(48-6)
dFVdcs   = (-p.F*diag(jF./i0.^2)*di0dcs_bar - (dexpIdcs_bar - dexpIIdcs_bar))...
                                             * m.Acs_bar;                        % Equation (48-2)
dFVdce   = (-p.F*diag(jF./i0.^2)*di0dce_bar) * m.Ace_bar;                        % Equation (48-3)
dFVdphis = p.F*diag(1./i0)*m.AjF_hat - (dexpIdphis - dexpIIdphis);               % Equation (48-4)
dFVdphie = (p.F*diag(1./i0)*m.BjF_hat - (dexpIdphie_bar - dexpIIdphie_bar))...
                                             * m.Aphie_bar;                      % Equation (48-5)
dFVdj    = p.F*diag(1./i0)*m.CjF_hat - (dexpIdj - dexpIIdj);                     % Equation (48-6)

% the Jacobian matrix: Equation (47-4)
J_x = [dFIdcs   dFIdce   dFIdphis   dFIdphie   dFIdj;...
       dFIIdcs  dFIIdce  dFIIdphis  dFIIdphie  dFIIdj;...
       dFIIIdcs dFIIIdce dFIIIdphis dFIIIdphie dFIIIdj;...
       dFIVdcs  dFIVdce  dFIVdphis  dFIVdphie  dFIVdj;...
       dFVdcs   dFVdce   dFVdphis   dFVdphie   dFVdj];

end
