clear;
clc;
close all;

% Define the problem parameters
L = 10; % Length of the domain
n = 8; % Number of mesh refinements or Fourier series terms
n_HO = n-2;

% Initialize arrays to store energies and mesh sizes/number of terms
energies_FE = zeros(n, 1);
energies_FD = zeros(n, 1);
energies_FT = zeros(n, 1);
energies_FE_HO = zeros(n, 1);
energies_FD_HO = zeros(n_HO, 1);

mesh_sizes_FE = zeros(n, 1);
mesh_sizes_FD = zeros(n, 1);
terms_FT = zeros(n, 1);
mesh_sizes_FE_HO = zeros(n, 1);
mesh_sizes_FD_HO = zeros(n_HO, 1);

mesh_num = zeros(n, 1);
for i = 1:n
    mesh_num(i) = 2^(i+1);
end

mesh_num_HO = zeros(n_HO, 1);
for i = 1:n_HO
    mesh_num_HO(i) = 2^(i+3);
end

% Function calls
tic;
[us_FE, xs_FE] = FE(L, n);
toc;
tic;
[us_FD, xs_FD] = FD(L, n);
toc;
tic;
[us_FT, xs_FT] = FT(L, n);
toc;
tic;
[us_FE_HO, xs_FE_HO] = FE_HO(L, n);
toc;
tic;
[us_FD_HO, xs_FD_HO] = FD_HO(L, n_HO);
toc;

% Plot for Finite Element Method
figure;
hold on;
title('Finite Element Method Solutions (2nd-order)');
xlabel('x');
ylabel('u');
ylim([-2, 2.5]);
for i = 1:n
    plot(xs_FE{i}, us_FE{i}, 'DisplayName', sprintf('mesh refinements = %d', 2^(i+1)));
end
legend('FontSize', 8, 'Location', 'south');
hold off;

% Plot for Finite Difference Method
figure;
hold on;
title('Finite Difference Method Solutions (2nd-order)');
xlabel('x');
ylabel('u');
ylim([-2, 2.5]);
for i = 1:n
    plot(xs_FD{i}, us_FD{i}, 'DisplayName', sprintf('mesh refinements = %d', 2^(i+1)));
end
legend('FontSize', 8, 'Location', 'south');
hold off;

% Plot for Fourier Series Method
figure;
hold on;
title('Fourier Series Method Solutions');
xlabel('x');
ylabel('u');
ylim([-2, 2.5]);
for i = 1:n
    plot(xs_FT{i}, us_FT{i}, 'DisplayName', sprintf('Fourier series terms = %d', 2^(i+1)));
end
legend('FontSize', 8, 'Location', 'south');
hold off;

% Plot for Finite Element Method HO
figure;
hold on;
title('Finite Element Method Solutions (3rd-order)');
xlabel('x');
ylabel('u');
ylim([-2, 2.5]);
for i = 1:n
    plot(xs_FE_HO{i}, us_FE_HO{i}, 'DisplayName', sprintf('mesh refinements = %d', 2^(i+1)));
end
legend('FontSize', 8, 'Location', 'south');
hold off;

% Plot for Finite Difference Method HO
figure;
hold on;
title('Finite Difference Method Solutions (4th-order)');
xlabel('x');
ylabel('u');
ylim([-2, 2.5]);
for i = 1:n_HO
    plot(xs_FE_HO{i}, us_FE_HO{i}, 'DisplayName', sprintf('mesh refinements = %d', 2^(i+3)));
end
legend('FontSize', 8, 'Location', 'south');
hold off;

% Calculate energies and mesh sizes/number of terms for each method
for i = 1:n
    energies_FE(i) = calculateEnergy(us_FE{i}, xs_FE{i}, L);
    energies_FD(i) = calculateEnergy(us_FD{i}, xs_FD{i}, L);
    energies_FT(i) = calculateEnergy(us_FT{i}, xs_FT{i}, L);
    energies_FE_HO(i) = calculateEnergy(us_FE_HO{i}, xs_FE_HO{i}, L);

    mesh_sizes_FE(i) = (2 * L) / length(xs_FE{i});
    mesh_sizes_FD(i) = (2 * L) / length(xs_FD{i});
    terms_FT(i) = 2^i;
    mesh_sizes_FE_HO(i) = (2 * L) / length(xs_FE_HO{i});
end

for i = 1:n_HO
    energies_FD_HO(i) = calculateEnergy(us_FD_HO{i}, xs_FD_HO{i}, L);

    mesh_sizes_FD_HO(i) = (2 * L) / length(xs_FD_HO{i});
end

% Calculate the energy errors for FEM and FDM by taking absolute differences
% For FT, since it's a series expansion, we take the absolute energy directly
energy_errors_FE = abs(energies_FE - [0; energies_FE(1:end-1)]);
energy_errors_FD = abs(energies_FD - [0; energies_FD(1:end-1)]);
energy_errors_FT = abs(energies_FT);
energy_errors_FE_HO = abs(energies_FE_HO - [0; energies_FE_HO(1:end-1)]);
energy_errors_FD_HO = abs(energies_FD_HO - [0; energies_FD_HO(1:end-1)]);

% Plot the energy errors in log-log scale
figure;
loglog(mesh_num, energy_errors_FE, 'o-', 'DisplayName', 'FEM 2rd');
hold on;
loglog(mesh_num, energy_errors_FD, 's-', 'DisplayName', 'FDM 2rd');
loglog(mesh_num, energy_errors_FT, 'x-', 'DisplayName', 'FT');
loglog(mesh_num, energy_errors_FE_HO, '*-', 'DisplayName', 'FEM 3rd');
loglog(mesh_num_HO, energy_errors_FD_HO, '^-', 'DisplayName', 'FDM 4th');
xlabel('Mesh Refinements / # Fourier terms');
ylabel('Energy Error');
title('Energy Error vs Mesh Refinements / # Fourier terms');
legend('Location', 'southwest');
grid on;

% Perform linear regression on log-log scale to determine convergence rates
% FE
log_mesh_sizes_FE = log(mesh_sizes_FE(2:end)); % Start from the second element to avoid log(0)
log_energy_errors_FE = log(energy_errors_FE(2:end));
p_FE = polyfit(log_mesh_sizes_FE, log_energy_errors_FE, 1);
slope_FE = p_FE(1);

% FD
log_mesh_sizes_FD = log(mesh_sizes_FD(2:end));
log_energy_errors_FD = log(energy_errors_FD(2:end));
p_FD = polyfit(log_mesh_sizes_FD, log_energy_errors_FD, 1);
slope_FD = p_FD(1);

% FT
log_terms_FT = log(terms_FT);
log_energy_errors_FT = log(energy_errors_FT);
p_FT = polyfit(log_terms_FT, log_energy_errors_FT, 1);
slope_FT = p_FT(1);

% FE_HO
log_mesh_sizes_FE_HO = log(mesh_sizes_FE_HO(2:end)); % Start from the second element to avoid log(0)
log_energy_errors_FE_HO = log(energy_errors_FE_HO(2:end));
p_FE_HO = polyfit(log_mesh_sizes_FE_HO, log_energy_errors_FE_HO, 1);
slope_FE_HO = p_FE_HO(1);

% FD_HO
log_mesh_sizes_FD_HO = log(mesh_sizes_FD_HO(2:end)); % Start from the second element to avoid log(0)
log_energy_errors_FD_HO = log(energy_errors_FD_HO(2:end));
p_FD_HO = polyfit(log_mesh_sizes_FD_HO, log_energy_errors_FD_HO, 1);
slope_FD_HO = p_FD_HO(1);

% Display the convergence rates
fprintf('Convergence rate for FEM 2rd: %f\n', slope_FE);
fprintf('Convergence rate for FDM 2rd: %f\n', slope_FD);
fprintf('Convergence rate for FT: %f\n', slope_FT);
fprintf('Convergence rate for FEM 3rd: %f\n', slope_FE_HO);
fprintf('Convergence rate for FDM 4th: %f\n', slope_FD_HO);
