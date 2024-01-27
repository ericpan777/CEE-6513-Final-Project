function [us, xs] = FE_HO(L, n)
% 3rd-order Finite Elements Method
% Initialize arrays to store mesh sizes, solutions, and x-coordinates
m = zeros(n, 1);
us = cell(1, n);
xs = cell(1, n);   % x-coordinates of global nodes and internal nodes of each element
x_Gs = cell(1, n); % x-coordinates of global nodes only

% Loop over different mesh sizes
for i = 1:n
    m(i) = 2^(i+1);

    % Initialize matrices for stiffness matrix K, force vector F, and solution vector u
    K = zeros(3 * m(i) + 1);
    F = zeros(3 * m(i) + 1, 1);
    u = zeros(3 * m(i) + 1, 1);
    x = zeros(3 * m(i) + 1, 1);
    x_G = zeros(m(i) + 1, 1);

    % Calculate the step size
    h = (2*L)/m(i);

    % Fill the x-coordinate array
    for j = 1:(m(i) + 1)
        x_G(j) = (-L) + (j-1)*h;
    end
    x_Gs{i} = x_G;

    for j = 1:m(i)
        % Define the element nodes
        x1 = x_G(j);
        x2 = x_G(j) + (x_G(j+1) - x_G(j))*(1/3);
        x3 = x_G(j) + (x_G(j+1) - x_G(j))*(2/3);
        x4 = x_G(j+1);

        x(3*j - 2) = x1;
        x(3*j - 1) = x2;
        x(3*j) = x3;

        % Define the symbolic variable for the interpolation functions
        syms x_N;

        % Define the 2nd-order interpolation functions
        N1 = ((x_N-x2)*(x_N-x3)*(x_N-x4))/((x1-x2)*(x1-x3)*(x1-x4));
        N2 = ((x_N-x1)*(x_N-x3)*(x_N-x4))/((x2-x1)*(x2-x3)*(x2-x4));
        N3 = ((x_N-x1)*(x_N-x2)*(x_N-x4))/((x3-x1)*(x3-x2)*(x3-x4));
        N4 = ((x_N-x1)*(x_N-x2)*(x_N-x3))/((x4-x1)*(x4-x2)*(x4-x3));

        % Calculate derivatives of the shape functions
        dN1dx = diff(N1, x_N);
        dN2dx = diff(N2, x_N);
        dN3dx = diff(N3, x_N);
        dN4dx = diff(N4, x_N);

        dN1dx_dN1dx = dN1dx*dN1dx;
        dN1dx_dN2dx = dN1dx*dN2dx;
        dN1dx_dN3dx = dN1dx*dN3dx;
        dN1dx_dN4dx = dN1dx*dN4dx;

        dN2dx_dN1dx = dN2dx*dN1dx;
        dN2dx_dN2dx = dN2dx*dN2dx;
        dN2dx_dN3dx = dN2dx*dN3dx;
        dN2dx_dN4dx = dN2dx*dN4dx;

        dN3dx_dN1dx = dN3dx*dN1dx;
        dN3dx_dN2dx = dN3dx*dN2dx;
        dN3dx_dN3dx = dN3dx*dN3dx;
        dN3dx_dN4dx = dN3dx*dN4dx;

        dN4dx_dN1dx = dN4dx*dN1dx;
        dN4dx_dN2dx = dN4dx*dN2dx;
        dN4dx_dN3dx = dN4dx*dN3dx;
        dN4dx_dN4dx = dN4dx*dN4dx;

        % Define the products of derivatives for stiffness matrix computation
        % (converting to functions for numerical integration)
        dN1dx_dN1dx_func = matlabFunction(dN1dx_dN1dx);
        dN1dx_dN2dx_func = matlabFunction(dN1dx_dN2dx);
        dN1dx_dN3dx_func = matlabFunction(dN1dx_dN3dx);
        dN1dx_dN4dx_func = matlabFunction(dN1dx_dN4dx);

        dN2dx_dN1dx_func = matlabFunction(dN2dx_dN1dx);
        dN2dx_dN2dx_func = matlabFunction(dN2dx_dN2dx);
        dN2dx_dN3dx_func = matlabFunction(dN2dx_dN3dx);
        dN2dx_dN4dx_func = matlabFunction(dN2dx_dN4dx);

        dN3dx_dN1dx_func = matlabFunction(dN3dx_dN1dx);
        dN3dx_dN2dx_func = matlabFunction(dN3dx_dN2dx);
        dN3dx_dN3dx_func = matlabFunction(dN3dx_dN3dx);
        dN3dx_dN4dx_func = matlabFunction(dN3dx_dN4dx);

        dN4dx_dN1dx_func = matlabFunction(dN4dx_dN1dx);
        dN4dx_dN2dx_func = matlabFunction(dN4dx_dN2dx);
        dN4dx_dN3dx_func = matlabFunction(dN4dx_dN3dx);
        dN4dx_dN4dx_func = matlabFunction(dN4dx_dN4dx);

        % Update the stiffness matrix K with the integrals of the products of derivatives
        % (element-wise stiffness matrix contribution)
        K((3*j-2), (3*j-2)) = K((3*j-2), (3*j-2)) + integral(dN1dx_dN1dx_func, x1, x4);
        K((3*j-2), (3*j-1)) = K((3*j-2), (3*j-1)) + integral(dN1dx_dN2dx_func, x1, x4);
        K((3*j-2), (3*j  )) = K((3*j-2), (3*j  )) + integral(dN1dx_dN3dx_func, x1, x4);
        K((3*j-2), (3*j+1)) = K((3*j-2), (3*j+1)) + integral(dN1dx_dN4dx_func, x1, x4);

        K((3*j-1), (3*j-2)) = K((3*j-1), (3*j-2)) + integral(dN2dx_dN1dx_func, x1, x4);
        K((3*j-1), (3*j-1)) = K((3*j-1), (3*j-1)) + integral(dN2dx_dN2dx_func, x1, x4);
        K((3*j-1), (3*j  )) = K((3*j-1), (3*j  )) + integral(dN2dx_dN3dx_func, x1, x4);
        K((3*j-1), (3*j+1)) = K((3*j-1), (3*j+1)) + integral(dN2dx_dN4dx_func, x1, x4);

        K((3*j  ), (3*j-2)) = K((3*j  ), (3*j-2)) + integral(dN3dx_dN1dx_func, x1, x4);
        K((3*j  ), (3*j-1)) = K((3*j  ), (3*j-1)) + integral(dN3dx_dN2dx_func, x1, x4);
        K((3*j  ), (3*j  )) = K((3*j  ), (3*j  )) + integral(dN3dx_dN3dx_func, x1, x4);
        K((3*j  ), (3*j+1)) = K((3*j  ), (3*j+1)) + integral(dN3dx_dN4dx_func, x1, x4);

        K((3*j+1), (3*j-2)) = K((3*j+1), (3*j-2)) + integral(dN4dx_dN1dx_func, x1, x4);
        K((3*j+1), (3*j-1)) = K((3*j+1), (3*j-1)) + integral(dN4dx_dN2dx_func, x1, x4);
        K((3*j+1), (3*j  )) = K((3*j+1), (3*j  )) + integral(dN4dx_dN3dx_func, x1, x4);
        K((3*j+1), (3*j+1)) = K((3*j+1), (3*j+1)) + integral(dN4dx_dN4dx_func, x1, x4);

        % Define the source term for each shape function and convert to function handles
        f1 = N1 * (exp(-x_N^2) - sqrt(pi)/(2*L));
        f2 = N2 * (exp(-x_N^2) - sqrt(pi)/(2*L));
        f3 = N3 * (exp(-x_N^2) - sqrt(pi)/(2*L));
        f4 = N4 * (exp(-x_N^2) - sqrt(pi)/(2*L));

        f1_func = matlabFunction(f1);
        f2_func = matlabFunction(f2);
        f3_func = matlabFunction(f3);
        f4_func = matlabFunction(f4);

        % Update the force vector F with the integrals of the source terms
        F(3*j-2) = F(3*j-2) + integral(f1_func, x1, x4);
        F(3*j-1) = F(3*j-1) + integral(f2_func, x1, x4);
        F(3*j)   = F(3*j)   + integral(f3_func, x1, x4);
        F(3*j+1) = F(3*j+1) + integral(f4_func, x1, x4);
    end
     % Apply periodic boundary conditions by wrapping the matrix and force vector
    K(1, :) = K(1, :) + K(end, :);
    K(:, 1) = K(:, 1) + K(:, end);
    F(1) = F(1) + F(end);

    % Remove the last row and column to account for the periodic boundary condition
    K(end, :) = [];
    K(:, end) = [];
    F(end) = [];

    % Solve the linear system to find the solution vector u
    u = pinv(K) * F;
    u(3 * m(i) + 1) = u(1);
    x(3 * m(i) + 1) = L;

    % Store the solutions and the x-coordinates
    us{i} = u;
    xs{i} = x;
end
