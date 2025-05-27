
% This function defines the system matrices that do not change in the
% inner loop, which involves the 10 coefficient matrices
% in Equations (36-1)-(36-4) and (39-1)-(39-3).

% Inputs:
%    - p is a MATLAB struct containing the initial parameters
%      and their corresponding vector form.

% Outputs:
%    - m is a MATLAB struct containing the system matrix,
%      such as the discretization matrix of the governing equations.

function m = fcn_system_matrices(p)

%% Mass conservation in the solid.
Acs_jn = sparse(fcn_Acs_j(p.nrn));   % Equation (20-3)
Acs_jp = sparse(fcn_Acs_j(p.nrp));
Acs_n = sparse(kron(speye(p.nn),Acs_jn)*diag((1/p.dr_n^2)*p.Ds_neg*ones(1,p.nn)...
                                             *kron(speye(p.nn),ones(1,p.nrn))));
Acs_p = sparse(kron(speye(p.np),Acs_jp)*diag((1/p.dr_p^2)*p.Ds_pos*ones(1,p.np)...
                                             *kron(speye(p.np),ones(1,p.nrp))));
m.Acs = blkdiag(Acs_n,Acs_p);        % Equation (21-4)

Bcs_jn = -2*(p.dr_n+p.rs_neg)/(p.dr_n*p.rs_neg)*[zeros(p.nrn-1,1); 1];   % Equation (20-4)
Bcs_jp = -2*(p.dr_p+p.rs_pos)/(p.dr_p*p.rs_pos)*[zeros(p.nrp-1,1); 1];
Bcs_n = sparse(kron(eye(p.nn), Bcs_jn));
Bcs_p = sparse(kron(eye(p.np), Bcs_jp));
m.Bcs = blkdiag(Bcs_n, Bcs_p);       % Equation (21-5)

%% Mass conservation in the electrolyte.
Bce = (diag(p.a_s.*(1-p.t_plus)./(p.eps_e(p.elec_range))));                    % Equations (30-5) and (30-6)
m.Bce = sparse([Bce(1:p.nn,:); zeros(p.ns,p.nnp); Bce(p.nn+1:p.nnp,:)]);       % Equation (34-7)

%% Charge conservation in the solid.
Aphis_neg = fcn_Aw(p.sigma_eff(1:p.nn),ones(p.nn,1),p.dx(1:p.nn));             % Equation (32-4)
Aphis_pos = fcn_Aw(p.sigma_eff(p.nn+1:p.nnp),ones(p.np,1),p.dx(p.nns+1:p.nx)); % Equation (32-5)
m.Aphis = sparse(blkdiag(Aphis_neg, Aphis_pos));                               % Equation (32-6)

m.Bphis = sparse(diag(-p.a_s*p.F.*p.dx(p.elec_range)));                        % Equation (32-7)
m.Cphis = sparse([1; zeros(p.nnp-2,1); -1]);                                   % Equation (32-8)

%% Charge conservation in the electrolyte.
Bphie = (diag(p.a_s*p.F.*p.dx(p.elec_range)));                                 % Equations (35-5) and (35-6)
m.Bphie = sparse([Bphie(1:p.nn,:); zeros(p.ns,p.nnp); Bphie(p.nn+1:p.nnp,:)]); % Equation (35-7)
m.Bphie(1) = 0;   % reference potential: φ_e1(x=0) = 0.                       % Equation (34)

%% used to calculate the parts of c_ss, c_e and φ_e given in the electrodes.
Acs_temp_n = kron(eye(p.nn),[zeros(1,p.nrn-1) 1]);
Acs_temp_p = kron(eye(p.np),[zeros(1,p.nrp-1) 1]);
m.Acs_bar = sparse(blkdiag(Acs_temp_n,Acs_temp_p));                 % Equation (39-4): A-_cs

Ace_temp = blkdiag(eye(p.nn), zeros(p.ns), eye(p.np));
m.Ace_bar = sparse([Ace_temp(1:p.nn,:); Ace_temp(p.nns+1:p.nx,:)]); % Equation (39-5): A-_ce

m.Aphie_bar = m.Ace_bar;                                            % Equation (39-6): A-_φe

end
