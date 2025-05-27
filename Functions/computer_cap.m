
% This function computes the ampere-hour total capacity of the cell
% and the 1C charge/discharge current density.

% Inputs:
%    - p is a MATLAB struct containing the initial parameters.
%    - mode is the selection switch for charging or discharging mode.

% Outputs:
%    - NP is a MATLAB struct containing the output variables,
%      such as the 1C charge/discharge current density.

function NP = computer_cap(p,mode)

% electrode plate area [m^2]
NP.A_surf = 1;

% ε_s: active material volume fraction at the electrodes [-]
epss_neg = 1 - p.epse_neg - p.epsf_neg;
epss_pos = 1 - p.epse_pos - p.epsf_pos;

% the battery ampere-hour capacity [A·h]
NP.Capa_neg =  (p.s100_neg-p.s0_neg)*epss_neg*p.L_neg*NP.A_surf*p.F*p.cs_max_neg/3600;  % Equation (85-1)
NP.Capa_pos = -(p.s100_pos-p.s0_pos)*epss_pos*p.L_pos*NP.A_surf*p.F*p.cs_max_pos/3600;  % Equation (85-2)

% the overall cell capacity [A·h]
NP.Capa = min(NP.Capa_neg,NP.Capa_pos);   % Equation (85-3)

% 1C (dis) charge current density [A/m^2]
NP.i_1C = NP.Capa / 1 / NP.A_surf;  % Equation (85-4)

% During the charge and discharge process, the electrode capacity is conserved,
% that is, the capacity consumed by one electrode is
% equal to the capacity accumulated by the other electrode.
% stoichiometries [-]
if NP.Capa_neg <= NP.Capa_pos
    NP.s100_neg = p.s100_neg;
    NP.s0_neg   = p.s0_neg;
    if mode == 1
        NP.s100_pos = p.s100_pos;
        NP.s0_pos   = p.s100_pos + NP.Capa_neg / (epss_pos*p.L_pos*NP.A_surf*p.F*p.cs_max_pos/3600);
    elseif mode == 2
        NP.s0_pos   = p.s0_pos;
        NP.s100_pos = p.s0_pos   - NP.Capa_neg / (epss_pos*p.L_pos*NP.A_surf*p.F*p.cs_max_pos/3600);
    end
elseif NP.Capa_neg > NP.Capa_pos
    NP.s100_pos = p.s100_pos;
    NP.s0_pos   = p.s0_pos;
    if mode == 1
        NP.s100_neg = p.s100_neg;
        NP.s0_neg   = p.s100_neg - NP.Capa_pos / (epss_neg*p.L_neg*NP.A_surf*p.F*p.cs_max_neg/3600);
    elseif mode == 2
        NP.s0_neg   = p.s0_neg;
        NP.s100_neg = p.s0_neg   + NP.Capa_pos / (epss_neg*p.L_neg*NP.A_surf*p.F*p.cs_max_neg/3600);
    end
end

end
