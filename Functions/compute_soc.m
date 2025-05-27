
% This function computes the state of charge (soc) of the full cell
% based on the solid phase concentration of the negative electrode.

% It is worth noting that if the capacity of the negative electrode
% is lower than that of the positive electrode, the total Li average
% concentration of the negative electrode is selected to calculate the SOC;
% if the capacity of the negative electrode is higher than that of the
% positive electrode, the total average concentration of the positive
% electrode is selected to calculate the SOC.

% Inputs:
%    - p is a MATLAB struct containing the initial parameters.
%    - cs is solid electrolyte, [mol/m^3].

% Outputs:
%    - soc is the state of charge (soc) of the full cell, [-].

function soc = compute_soc(p,cs,NP)

if NP.Capa_neg <= NP.Capa_pos

    % negative electrode
    ri_n = 1:p.nrn;
    A_mean_n = (1/p.nn)*kron(ones(1,p.nn),eye(p.nrn));   % cs(x,r=1:nr,t)
    cs_avg_r_n = A_mean_n*cs(1:p.nn*p.nrn);              % p.nrn×1, similar to Equation (86-2)
    cs_f_n = @(r) interp1(((ri_n-1)*p.dr_n),cs_avg_r_n,r,'linear','extrap');

    n_int = 100;
    r_neg = linspace(0,p.rs_neg,n_int);
    dr_n = p.rs_neg/(n_int-1);

    fx_kp1_n = diag(r_neg(2:n_int).^2)*cs_f_n(r_neg(2:n_int))';     % rs^2 * cs
    fx_k_n = diag(r_neg(1:n_int-1).^2)*cs_f_n(r_neg(1:n_int-1))';   % rs^2 * cs

    % total average Li concentration in the solid particle, [mol/m^3].
    cs_avg_n = 3/p.rs_neg^3 * sum((fx_kp1_n+fx_k_n)/2*dr_n);      % similar to Equation (86-1)

    soc = (cs_avg_n/p.cs_max_neg-p.s0_neg)/(p.s100_neg-p.s0_neg); % Equation (86-4)

elseif NP.Capa_neg > NP.Capa_pos

    % positive electrode
    ri_p = 1:p.nrp;
    A_mean_p = (1/p.np)*kron(ones(1,p.np),eye(p.nrp));                   % cs(x,r=1:nr,t)
    cs_avg_r_p = A_mean_p*cs(p.nn*p.nrn+1:p.nrn*p.nn+p.nrp*p.np);        % p.nrp×1, similar to Equation (86-3)
    cs_f_p = @(r) interp1(((ri_p-1)*p.dr_p),cs_avg_r_p,r,'linear','extrap'); 

    % n_int = 100;
    r_pos = linspace(0,p.rs_pos,n_int);
    dr_p = p.rs_pos/(n_int-1);

    fx_kp1_p = diag(r_pos(2:end).^2)*cs_f_p(r_pos(2:end))';     % rs^2 * cs
    fx_k_p = diag(r_pos(1:end-1).^2)*cs_f_p(r_pos(1:end-1))';   % rs^2 * cs

    % total average Li concentration in the solid particle, [mol/m^3].
    cs_avg_p = 3/p.rs_pos^3 * sum((fx_kp1_p+fx_k_p)/2*dr_p);          % similar to Equation (86-1)

    soc = (cs_avg_p/p.cs_max_pos-p.s100_pos)/(p.s0_pos-p.s100_pos);   % Equation (86-4)

end

end
