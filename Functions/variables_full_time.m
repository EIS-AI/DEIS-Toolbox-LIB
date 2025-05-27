
% This function computes the variables on the interface of control volumes
% by using continuity boundary conditions based on the variables
% at the center of the corresponding control volume.

% Inputs:
%    - p is a MATLAB struct containing the initial parameters
%      and their corresponding vector form.
%    - i_app is the applied (dis) charge current density, [A/m^2].
%    - ce is electrolyte electrolyte, [mol/m^3].
%    - phie is electrolyte potential, [V].
%    - phis is solid potential, [V].
%    - n_t is the number of simulation steps.

% Outputs:
%    - ce_full is electrolyte electrolyte, [mol/m^3].
%    - phie_full is electrolyte potential, [V].
%    - phis_full is solid potential, [V].
%    - ie_full is electrolyte current density, [A/m^2].

function [ce_full,phie_full,phis_full,ie_full] = variables_full_time(p,i_app,ce,phie,phis,n_t)

%% ce_full
De_eff = p.De_eff(ce,p.T);
D_delta = De_eff(2:end,:) ./ De_eff(1:end-1,:) .* p.dx(1:end-1) ./ p.dx(2:end);
ce_middle = (ce(1:end-1,:) + D_delta.*ce(2:end,:)) ./ (1 + D_delta);  % Equation (51-1)
ce_full = [ce(1,:); ce_middle; ce(end,:)];

%% φe_full
nu = p.nu(ce,p.T);
kappa_eff = p.kappa_eff(ce,p.T);

nu_delta = nu(2:end,:)   .* (log(ce(2:end,:)) - log(ce_middle))     ./ p.dx(2:end)...
         - nu(1:end-1,:) .* (log(ce_middle)   - log(ce(1:end-1,:))) ./ p.dx(1:end-1);
kappa_delta = kappa_eff(2:end,:) ./ kappa_eff(1:end-1,:) .* p.dx(1:end-1) ./ p.dx(2:end);

phie_middle = (phie(1:end-1,:) + kappa_delta .* phie(2:end,:)...
             + p.dx(1:end-1) ./ kappa_eff(1:end-1,:) .* nu_delta) ./ (1 + kappa_delta);  % Equation (51-2)
phie_full = [phie(1,:); phie_middle; phie(end,:)];

% verify boundary condition
% i_app_ver = - kappa_eff(p.nn,:) .* (phie_middle(p.nn,:)    - phie(p.nn,:))    ./ p.dx(p.nn,:)*2 ...
%             - nu(p.nn,:)        .* (log(ce_middle(p.nn,:)) - log(ce(p.nn,:))) ./ p.dx(p.nn,:)*2;

%% φs_full
sigma_delta_neg = p.sigma_eff_neg(2:end) ./ p.sigma_eff_neg(1:end-1)...
               .* p.dx(1:p.nn-1)         ./ p.dx(2:p.nn);
sigma_delta_pos = p.sigma_eff_pos(2:end) ./ p.sigma_eff_pos(1:end-1)...
               .* p.dx(p.nns+1:end-1)    ./ p.dx(p.nns+2:end);

phis_middle_neg = (phis(1:p.nn-1,:)     + sigma_delta_neg.*phis(2:p.nn,:))     ./ (1 + sigma_delta_neg);
phis_middle_pos = (phis(p.nn+1:end-1,:) + sigma_delta_pos.*phis(p.nn+2:end,:)) ./ (1 + sigma_delta_pos);

phis_full = [phis(1,:)+i_app*p.dx(1)/2/p.sigma_eff(1);...
             phis_middle_neg; phis(p.nn,:);...
             NaN(p.ns-1,n_t-1);
             phis(p.nn+1,:); phis_middle_pos;...
             phis(end,:)+i_app*p.dx(end)/2/p.sigma_eff(end)];     % Equation (51-3) and (51-4)

%% ie_full
kappa_eff = p.kappa_eff(ce,p.T);
nu = p.nu(ce,p.T);

% harmonic mean (HM)
kappa_eff_middle = kappa_eff(2:end,:) .* kappa_eff(1:end-1,:) .* (p.dx(2:end)+p.dx(1:end-1))...
               ./ (kappa_eff(2:end,:) .* p.dx(1:end-1) + kappa_eff(1:end-1,:) .* p.dx(2:end));  % Equation (22)
nu_middle = nu(2:end,:) .* nu(1:end-1,:) .* (p.dx(2:end)+p.dx(1:end-1))...
        ./ (nu(2:end,:) .* p.dx(1:end-1) + nu(1:end-1,:) .* p.dx(2:end));                       % Equation (22)

% central difference
dlncedx_middle = (log(ce(2:end,:))-log(ce(1:end-1,:))) ./ ((p.dx(2:end)+p.dx(1:end-1))/2);
dphiedx_middle =   (phie(2:end,:)-   phie(1:end-1,:))  ./ ((p.dx(2:end)+p.dx(1:end-1))/2);

% electrolyte current density
ie_O = - kappa_eff_middle .* dphiedx_middle;
ie_D = -        nu_middle .* dlncedx_middle;
ie_full = [zeros(1,n_t-1); ie_O+ie_D; zeros(1,n_t-1)];

end
