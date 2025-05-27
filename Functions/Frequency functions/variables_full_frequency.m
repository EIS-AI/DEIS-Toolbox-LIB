
% This function computes the variables on the interface of control volumes
% by using continuity boundary conditions based on the variables
% at the center of the corresponding control volume.

% Inputs:
%    - p is a MATLAB struct containing the initial parameters
%      and their corresponding vector form.
%    - ce is electrolyte electrolyte, [mol/m^3].
%    - out is a MATLAB struct containing the internal variables,
%      such as the concentrations and the potentials.
%    - Model is the selection switch for the complete DFN (CDFN) model or simplified model.

% Outputs:
%    - out is a MATLAB struct containing the internal variables,
%      such as the concentrations and the potentials.

function out = variables_full_frequency(p,ce,out,Model)

kappa_eff = p.kappa_eff(ce,p.T);
kappa_delta = kappa_eff(2:end,:) ./ kappa_eff(1:end-1,:) .* p.dx(1:end-1,:) ./ p.dx(2:end,:);

% about Ie_middle
% harmonic mean (HM)
kappa_eff_middle = kappa_eff(2:end,:) .* kappa_eff(1:end-1,:) .* (p.dx(2:end)+p.dx(1:end-1))...
               ./ (kappa_eff(2:end,:) .* p.dx(1:end-1) + kappa_eff(1:end-1,:) .* p.dx(2:end));  % Equation (22)
% central difference
dPhiedx_middle = (out.Phie(2:end,:)-out.Phie(1:end-1,:)) ./ ((p.dx(2:end)+p.dx(1:end-1))/2);
Ie_O = - kappa_eff_middle .* dPhiedx_middle;

if Model == 'DFN'

    % Ce_full
    De_eff = p.De_eff(ce,p.T);
    D_delta = De_eff(2:end,:) ./ De_eff(1:end-1,:) .* p.dx(1:end-1,:) ./ p.dx(2:end,:);
    Ce_middle = (out.Ce(1:end-1,:) + D_delta.*out.Ce(2:end,:)) ./ (1 + D_delta);  % Equation (79-1)
    out.Ce_full = [out.Ce(1,:); Ce_middle; out.Ce(end,:)];

    % Φe_middle
    nu = p.nu(ce,p.T);
    nu_delta = nu(2:end,:)  ./ce(2:end,:)   .* (out.Ce(2:end,:) - Ce_middle)     ./ p.dx(2:end,:)...
             - nu(1:end-1,:)./ce(1:end-1,:) .* (Ce_middle   - out.Ce(1:end-1,:)) ./ p.dx(1:end-1,:);
    Phie_middle = (out.Phie(1:end-1,:) + kappa_delta .* out.Phie(2:end,:)...
                 + p.dx(1:end-1,:) ./ kappa_eff(1:end-1,:) .* nu_delta) ./ (1 + kappa_delta);  % Equation (79-2)

    % Ie_middle
    nu_ce = p.nu(ce,p.T) ./ ce;
    nu_ce_middle = nu_ce(2:end,:) .* nu_ce(1:end-1,:) .* (p.dx(2:end)+p.dx(1:end-1))...
               ./ (nu_ce(2:end,:) .* p.dx(1:end-1) + nu_ce(1:end-1,:) .* p.dx(2:end));              % Equation (22)
    dlnCedx_middle = (  out.Ce(2:end,:)-  out.Ce(1:end-1,:)) ./ ((p.dx(2:end)+p.dx(1:end-1))/2);
    Ie_D = - nu_ce_middle .* dlnCedx_middle;
    out.Ie_middle = Ie_O + Ie_D;

elseif Model == 'sim'
    
    % Φe_middle
    Phie_middle = (out.Phie(1:end-1,:) + kappa_delta .* out.Phie(2:end,:)) ./ (1 + kappa_delta);  % Equation (79-2)

    % Ie_middle
    out.Ie_middle = Ie_O;

end

% Φe_full
out.Phie_full = full([out.Phie(1,:); Phie_middle; out.Phie(end,:)]);

% Φs_full
sigma_delta_neg = p.sigma_eff_neg(2:end) ./ p.sigma_eff_neg(1:end-1)...
               .* p.dx(1:p.nn-1)         ./ p.dx(2:p.nn);
sigma_delta_pos = p.sigma_eff_pos(2:end) ./ p.sigma_eff_pos(1:end-1)...
               .* p.dx(p.nns+1:end-1)    ./ p.dx(p.nns+2:end);
Phis_middle_neg = (out.Phis(1:p.nn-1,:)     + sigma_delta_neg.*out.Phis(2:p.nn,:))     ./ (1 + sigma_delta_neg);
Phis_middle_pos = (out.Phis(p.nn+1:end-1,:) + sigma_delta_pos.*out.Phis(p.nn+2:end,:)) ./ (1 + sigma_delta_pos);
out.Phis_full = full([out.Phis(1,1)+1*p.dx(1)/2/p.sigma_eff(1);...
                      Phis_middle_neg; out.Phis(p.nn,1);...
                      NaN(p.ns-1,1);
                      out.Phis(p.nn+1,1); Phis_middle_pos;...
                      out.Phis(end,1)+1*p.dx(end)/2/p.sigma_eff(end)]);     % Equation (79-3) and (79-4)

end
