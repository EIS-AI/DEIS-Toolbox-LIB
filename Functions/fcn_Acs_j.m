
% This function computes the matrix A_cs at the position x_j in the electrode.

% Inputs:
%    - nr is the number of grid points along the radial coordinate, [-].

% Outputs:
%    - Acs_j is the matrix A_cs at the position x_j for c_s.

function Acs_j = fcn_Acs_j(nr)
temp1 = linspace(0,nr-2,nr-1)./linspace(1,nr-1,nr-1);
temp1(nr-1) = 2;
temp2 = [2 linspace(2,nr-1,nr-2)./linspace(1,nr-2,nr-2)];

Acs_j = diag(temp1,-1)...
      - 2*eye(nr)...
      + diag(temp2,1);  % Equation (20-3)

end
