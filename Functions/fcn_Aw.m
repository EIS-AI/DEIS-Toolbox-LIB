
% This function computes the matrix A_w in the full cell or battery components.

% Inputs:
%    - alpha1 is α_I, refers to D_e,eff, κ_eff, κ_D,eff or σ_eff.
%    - alpha0 is α_O×δ_x, see to Equation (28-1).
%    - dx is interval size along the x-dimension.

% Outputs:
%    - Aw is the matrix A_w for c_e, φ_s and φ_e.

function Aw = fcn_Aw(alpha1,alpha0,dx)
% harmonic mean: Equation (22)
alpha1_face= alpha1(2:end).*alpha1(1:end-1).*(dx(2:end)+dx(1:end-1))...
          ./(alpha1(2:end).*dx(1:end-1)+alpha1(1:end-1).*dx(2:end));

% linear interpolation technique: similar to Equation (27-3)
dwdx_face = 2./(dx(2:end)+dx(1:end-1));

% approximate the fluxes on the surfaces of the control volume element:
% similar to Equation (27-4), also see to Equations (30-8), (32-9) and (35-9).
gamma = alpha1_face.*dwdx_face;

Am = diag(gamma,-1)...
   + diag([-gamma(1); -gamma(2:end)-gamma(1:end-1); -gamma(end)])...
   + diag(gamma,1);  % see to Equation (28-1)
Aw = diag(1./(alpha0)) * Am;

end
