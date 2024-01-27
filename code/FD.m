function [us, xs] = FD(L, n)
% 2nd-order Finite Difference Method
% Initialize arrays for mesh sizes, solutions, and x-coordinates
m = zeros(n, 1);
h = zeros(n, 1);
us = cell(1, n);
xs = cell(1, n);

% Loop over mesh refinements to compute numerical solutions
for i = 1:n
    % Calculate number of internal nodes and mesh size
    m(i) = 2^(i+1);
    h(i) = (2*L)/m(i);

    % Initialize matrices and vectors
    A = zeros(m(i));        % Stiffness matrix
    B = zeros(m(i), 1);     % Load vector
    x = zeros(m(i) + 1, 1); % x coordinates

    % Populate A matrix for internal nodes
    % Setting up the tridiagonal matrix with periodic boundary conditions
    A(1, 1) = -2;
    A(1, 2) = 1;
    A(1, m(i)) = 1;

    for j = 2:m(i)-1
        A(j, j-1) = 1;
        A(j, j) = -2;
        A(j, j+1) = 1;
    end

    A(m(i), m(i)-1) = 1;
    A(m(i), m(i)) = -2;
    A(m(i), 1) = 1;

    A = A/(h(i)^2); % Scale A matrix by mesh size squared

    % Populate B vector with the source term
    for j = 1:m(i)
        x(j) = -L + (j-1)*h(i);
        B(j) = -(exp(-x(j)^2) - sqrt(pi)/(2*L));
    end

    % Solve the linear system A*u = B
    u = pinv(A) * B;
    u(end+1) = u(1); % Ensure periodicity

    x(m(i) + 1) = L;

    % Store the solutions and the x-coordinates
    us{i} = u;
    xs{i} = x;
end