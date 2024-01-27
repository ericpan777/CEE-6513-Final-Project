function E = calculateEnergy(u, x, L)
% calculateEnergy: Calculate the energy of a numerical solution to a PDE
% Inputs:
%   u - Numerical solution vector
%   x - Vector of spatial coordinates corresponding to the solution u
%   L - Half-length of the domain
%
% Output:
%   E - Calculated energy of the numerical solution

% Number of points in the solution
N = length(u);

% Calculate the spatial step size
dx = (2 * L) / (N - 1);

% Initialize the derivative of u (dudx) with zeros
dudx = zeros(size(u));

% Calculate the derivative using central differences
% For interior points
dudx(2:end-1) = (u(3:end) - u(1:end-2)) / (2 * dx);

% For the first and last points (accounting for periodic boundary conditions)
dudx(1) = (u(2) - u(end-1)) / (2 * dx);
dudx(end) = (u(1) - u(end-1)) / (2 * dx);

% For the first and last points (accounting for periodic boundary conditions)
f_x = exp(-x.^2) - sqrt(pi)/(2*L);

% Calculate the energy using the trapezoidal rule for numerical integration
% The energy functional includes a term for the gradient (0.5 * dudx^2)
% and a term for the potential energy (f_x .* u)
E_u = trapz(x, 0.5 * dudx.^2 - f_x .* u);

% Return the calculated energy
E = E_u(1);
end
