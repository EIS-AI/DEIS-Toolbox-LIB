
% This function defines the system matrices that do change in the inner loop.

% Inputs:
%    - p is a MATLAB struct containing the initial parameters
%      and their corresponding vector form.
%    - m is a MATLAB struct containing the system matrix,
%      such as the discretization matrix of the governing equations.
%    - ce is electrolyte electrolyte, [mol/m^3].
%    - T is ambient temperature, [K].
%    - phis_prevt is solid potential of the previous moment, [V].
%    - phie_prevt is electrolyte potential of the previous moment, [V].
%    - phisf_prevt is potential in the film close to the solid phase
%      of the previous moment, [V].

% Outputs:
%    - m is a MATLAB struct containing the system matrix,
%      such as the discretization matrix of the governing equations.

function m = fcn_system_matrices_change(p,m,ce,T,phis_prevt,phie_prevt,phisf_prevt)

Ace = sparse(fcn_Aw(p.De_eff(ce,T),p.eps_e.*p.dx,p.dx));         % Equation (30-4)

m.Aphie = sparse(fcn_Aw(p.kappa_eff(ce,T),ones(p.nx,1),p.dx));   % Equation (35-4)
m.Aphie(1,1:2) = [1 0];   % reference potential: φ_e1(x=0) = 0. % Equation (34)

m.Dphie = sparse(fcn_Aw(p.nu(ce,T),ones(p.nx,1),p.dx));          % Equation (35-8)
m.Dphie(1,1:2) = [0 0];   % reference potential: φ_e1(x=0) = 0. % Equation (34)

% temporal discretization
m.Acs_hat = p.dt*m.Acs - speye(p.nrt);    % Equation (42-1): A^_cs
m.Bcs_hat = p.dt*m.Bcs;                   % Equation (42-2): B^_cs
m.Ace_hat = p.dt*  Ace - speye(p.nx);     % Equation (42-3): A^_ce
m.Bce_hat = p.dt*m.Bce;                   % Equation (42-4): B^_ce

% the parts of φ_e(t_k-1) given in the electrodes: Equation (39-3)
phie_bar_prevt = m.Aphie_bar*phie_prevt;

% Equations (44-1)-(44-4)
m.Aphisf_hat = sparse(- 1 / p.dt * diag(p.Rsei) * diag(p.Csei));
m.Bphisf_hat = sparse(  1 / p.dt * diag(p.Rsei) * diag(p.Csei) + eye(p.nnp));
m.Cphisf_hat = sparse(diag(p.F * p.Rsei));
m.Dphisf_hat = sparse(  1 / p.dt * diag(p.Rsei) * diag(p.Csei) * (phis_prevt-phie_bar_prevt));

% Equations (44-5)-(44-8)
m.Ajdl_hat = sparse(  1 / p.dt * diag(p.Cdl/p.F) * (eye(p.nnp) - m.Aphisf_hat));
m.Bjdl_hat = sparse(- 1 / p.dt * diag(p.Cdl/p.F) * m.Bphisf_hat);
m.Cjdl_hat = sparse(- 1 / p.dt * diag(p.Cdl/p.F) * m.Cphisf_hat);
m.Djdl_hat = sparse(- 1 / p.dt * diag(p.Cdl/p.F) * (phis_prevt-phisf_prevt)...
                    - 1 / p.dt * diag(p.Cdl/p.F) * m.Dphisf_hat);

% Equations (44-9)-(44-12)
m.AjF_hat = sparse(- 1 / p.dt * diag(p.Csei/p.F) - m.Ajdl_hat);
m.BjF_hat = sparse(  1 / p.dt * diag(p.Csei/p.F) - m.Bjdl_hat);
m.CjF_hat = sparse(eye(p.nnp) - m.Cjdl_hat);
m.DjF_hat = sparse(  1 / p.dt * diag(p.Csei/p.F) * (phis_prevt-phie_bar_prevt) - m.Djdl_hat);

end
