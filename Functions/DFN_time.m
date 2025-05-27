
% This function simulates the time domain charge/discharge curve of the DFN model.

% Inputs:
%    - p is a MATLAB struct containing the initial parameters.
%    - soc_init is the initial state of charge (soc) of the full cell, [-].
%    - Crate is the charge/discharge rate, [-].
%    - mode is the selection switch for charging or discharging mode.

% Outputs:
%    - volt is a MATLAB struct containing the output variables,
%      such as the output voltage.
%    - vari is a MATLAB struct containing the internal variables,
%      such as the concentrations and the potentials.

function [volt,vari] = DFN_time(p,soc_init,Crate,mode)

%% Define variables for simulation.
p.Vmin = 2.70;      % minimum terminal voltages [V]
p.Vmax = 4.20;      % maximum terminal voltages [V]
% p.OCV_min = p.U_pos(p.s0_pos)  -p.U_neg(p.s0_neg);    2.7010 in the present study
% p.OCV_max = p.U_pos(p.s100_pos)-p.U_neg(p.s100_neg);  4.1618 in the present study

p = fcn_system_vectors(p);
m = fcn_system_matrices(p);

[cs,ce,phis,phie,j] = variable_init(p,soc_init);
i_app(1) = 0;
v_cell(1) = phis(p.nnp) - phis(1);  % Equation (50-4)
soc(1) = compute_soc(p,cs,p.NP);

%% Newton's method related parameters
% Measure simulation time for benchmark purposes.
dispstat('Starting the time domain simulation...','keepthis','timestamp')

% prev indicates the condition of the state at k-1.
x_prev = [cs; ce; phis; phie; j];  % p.nt×1: Equation (47-3)

% the value of the previous moment
cs_prevt = cs(:,1);            % p.nrt×1
ce_prevt = ce(:,1);            % p.nx ×1
phis_prevt = phis(:,1);        % p.nnp×1
phie_prevt = phie(:,1);        % p.nx ×1
phisf_prevt = zeros(p.nnp,1);  % p.nnp×1

% Parameters that remain unchanged in Newton's method.
Nm.tol = 1e-3;           % tolerance for convergence
Nm.iter_max = 1e2;       % Maximum iterations for the Newton's algorithm for solving the algebraic equations.
Nm.end_simulation = 0;   % Flag indicating that the end of simulation.

% Parameters that vary in Newton's method.
t_vec = 0;               % 𝑡_𝑘 = 𝑘·Δ𝑡
t = 1;                   % 𝑘 ∈ {1, …, 𝑁}, in which 𝑁 is the number of simulation steps.
solution_time = 0;       % the time taken by Newton's method

%% Simulation
while not(Nm.end_simulation)
% Inner loop
    tic()

    % charging/discharging time at time t
    t_vec(t+1) = t_vec(t) + p.dt;

    % charging/discharging current at time t: Equation (17)
    if mode == 1
        i_app(t+1) =   p.NP.i_1C * Crate * (1 + 1/50 * sin(2*pi*0*t_vec(t+1)));
    elseif mode == 2
        i_app(t+1) = - p.NP.i_1C * Crate * (1 + 1/50 * sin(2*pi*0*t_vec(t+1)));
    end

% Solve set of algebraic equations (AEs) using Newton's method at time t.
    for k = 1:Nm.iter_max
        % obtain function matrix: Equation (47-2)
        [F_x,cs(:,t+1),ce(:,t+1),phis(:,t+1),phie(:,t+1),j(:,t+1),...
         phisf(:,t+1),jF(:,t+1),eta(:,t+1),Uocp(:,t+1),m]...
       = fcn_F(p,m,x_prev,cs_prevt,ce_prevt,phis_prevt,phie_prevt,phisf_prevt,i_app(t+1));

        % norm(A,2) = sqrt(sum(abs(A.^2)))
        conv_check = norm(F_x,2);

        % If criterium for convergence is met, then break loop.
        if conv_check < Nm.tol  % Equation (49)
            break
        end

        % obtain Jacobian matrix: Equation (47-4)
        J_x = fcn_J(p,m,cs(:,t+1),ce(:,t+1),jF(:,t+1),eta(:,t+1));

        % solving using Newton's method: Equation (47-1)
        x_prev = real(x_prev - (J_x\F_x));

    end

    v_cell(t+1) = phis(p.nnp,t+1)+i_app(t+1)*p.dx(end)/2/p.sigma_eff(p.nnp)...
                - phis(1,t+1)    -i_app(t+1)*p.dx(1)  /2/p.sigma_eff(1);  % Equation (50-4)
    soc(t+1) = compute_soc(p,cs(:,t+1),p.NP);
    stoich(:,t+1) = [cs(p.nrn:p.nrn:p.nrn*p.nn,t+1)/p.cs_max_neg;...
                     cs(p.nrn*p.nn+p.nrp:p.nrp:p.nrt,t+1)/p.cs_max_pos];

    solution_time = solution_time + toc();

    if (mode == 1 && v_cell(t+1)<p.Vmin) || (mode == 2 && v_cell(t+1)>p.Vmax)
        warning('Voltage = %d exceeded specified bounds. Stopping the simulation.',v_cell(t+1))
        break
    end
    if any(stoich(:,t+1)>1) || any(stoich(:,t+1)<0)
        warning('cs exceeeding either cs_max or is lower than 0. Stopping the simulation.')
        break
    end
    if any(ce(:,t+1)<0)
        warning('ce is lower than 0. Stopping the simulation.')
        break
    end

    % the value of the previous moment
    cs_prevt = cs(:,t+1);
    ce_prevt = ce(:,t+1);
    phis_prevt = phis(:,t+1);
    phie_prevt = phie(:,t+1);
    phisf_prevt = phisf(:,t+1);
    
    t = t+1;
end

dispstat(sprintf('At rate = %2.2f, finished the time domain simulation in %2.2f s \n',Crate,solution_time),'keepthis','timestamp');

%% Store states
vari.x = [p.dx_n/2*(1:2:(p.nn*2-1))...
          p.L_neg+p.dx_s/2*(1:2:(p.ns*2-1))...
          p.L_neg+p.L_sep+p.dx_p/2*(1:2:(p.np*2-1))]'...
       / (p.L_neg+p.L_sep+p.L_pos);

n_t = t+1;

vari.cs = cs(:,2:n_t);
vari.ce = ce(:,2:n_t);
vari.phis = [phis(1:p.nn,2:n_t); NaN(p.ns,n_t-1); phis(p.nn+1:p.nnp,2:n_t)];
vari.phie = phie(:,2:n_t);
vari.j    = [j(1:p.nn,2:n_t);    NaN(p.ns,n_t-1); j(p.nn+1:p.nnp,2:n_t)];

vari.U    = [Uocp(1:p.nn,2:n_t); NaN(p.ns,n_t-1); Uocp(p.nn+1:p.nnp,2:n_t)];
vari.eta  = [eta(1:p.nn,2:n_t);  NaN(p.ns,n_t-1); eta(p.nn+1:p.nnp,2:n_t)];
[vari.ce_full,vari.phie_full,vari.phis_full,vari.ie_full]...
    = variables_full_time(p,i_app(2:n_t),ce(:,2:n_t),phie(:,2:n_t),phis(:,2:n_t),n_t);

volt.Crate = Crate;
volt.p = p;
volt.t = t_vec(2:n_t);
volt.i_app = i_app(2:n_t);
volt.soc = soc(2:n_t);

volt.v_neg = vari.phie_full(p.nn+1,:)  - vari.phis_full(1,:);         % Equation (50-1)
volt.v_pos = vari.phis_full(p.nnp+1,:) - vari.phie_full(p.nns+1,:);   % Equation (50-2)
volt.v_sep = vari.phie_full(p.nns+1,:) - vari.phie_full(p.nn+1,:);    % Equation (50-3)
volt.v_cell = v_cell(2:n_t);

% % stoichiometric coefficient for porous electrode [-]
% w = p.NP.s0_neg + soc(2:n_t) * (p.NP.s100_neg - p.NP.s0_neg);
% z = p.NP.s0_pos + soc(2:n_t) * (p.NP.s100_pos - p.NP.s0_pos);
% % the equilibrium potential of the full cell and its components [V]
% volt.Eq_neg  = - p.U_neg(w);
% volt.Eq_pos  = p.U_pos(z);
% volt.Eq_sep  = zeros(1,t);
% volt.Eq_cell = p.U_pos(z) - p.U_neg(w);

end
