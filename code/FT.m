function [us, xs] = FT(L, n)
% Fourier Transform Method
% x-coordinates at which to evaluate the solution
x = linspace(-L, L, 1000);

% Initialize cell arrays to store solutions and x-coordinates
us = cell(1, n);
xs = cell(1, n);

% Define the function on the right-hand side of the PDE
f = @(x) exp(-x.^2) - sqrt(pi)/(2*L);

% Loop over different numbers of Fourier series terms
for i = 1:n
    N = 2^(i+1); % Number of terms in the Fourier series

    % Initialize the solution u(x)
    u = zeros(size(x));

    % Compute and sum the Fourier series terms
    for j = -N:N
        % Skip j = 0 to avoid division by zero
        if j == 0
            continue;
        end

        % Compute the jth Fourier coefficient using numerical integration
        integrand = @(x) f(x) .* exp(-1i * j * pi * x / L);
        a_n = integral(integrand, -L, L) / (2*L);

        % Add the jth term to the solution
        u = u + a_n * exp(1i * j * pi * x / L) / ((j * pi / L)^2);
    end

    % Take the real part of the solution (the imaginary part should be negligible)
    u = real(u);
    u = transpose(u);
    us{i} = u; % Store the solution
    xs{i} = x; % Store the x-coordinates
end
