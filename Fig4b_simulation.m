clearvars; clc; close all; 

%% Define parameters
g_TR = 0.4; gamma0 = 0.1;
U1 = [1, -1; 1, 1]/sqrt(2); U2 = [1, 0; 0, -1]; 
psi_in_list = [ 1/sqrt(2), 1, 1/sqrt(2), 0,  1/sqrt(2),   1/sqrt(2);
               1i/sqrt(2), 0, 1/sqrt(2), 1, -1/sqrt(2), -1i/sqrt(2)];
var_phi = -pi/2; % Change as needed

%% Calculate spectrum
num_pts_within_roundtrip = 250; num_roundtrips = 800;
spectrum_raw = zeros(num_pts_within_roundtrip * num_roundtrips, 6); 
Omegat_list = linspace(0, num_roundtrips, size(spectrum_raw, 1) + 1) * 2*pi; Omegat_list(end) = [];
delta_omega_list = linspace(-1, 1, size(spectrum_raw, 1) + 1) * 0.2 * 2*pi; delta_omega_list(end) = [];

T = @(x) transmittance(x, g_TR, var_phi);

for psi_in_index = 1 : 6
    psi_in = psi_in_list(:, psi_in_index);

    for loop_index = 1 : size(spectrum_raw, 1)
        Omegat = Omegat_list(loop_index);
        delta_omega = delta_omega_list(loop_index);
        P_00 =  exp(1i * delta_omega) * exp(-gamma0);
        P_05 = -exp(1i * delta_omega) * exp(-gamma0);

        psi_ss_00 = inv(eye(2) - U1 * T(Omegat) * U2 * U1 * T(Omegat - 2*pi) * U2 * P_00^2) * ...
                    (eye(2) + P_00 * U1 * T(Omegat) * U2) * psi_in;
        psi_ss_05 = inv(eye(2) - U1 * T(Omegat) * U2 * U1 * T(Omegat - 2*pi) * U2 * P_05^2) * ...
                    (eye(2) + P_05 * U1 * T(Omegat) * U2) * psi_in;
        
        spectrum_raw(loop_index, psi_in_index) = ...
            psi_ss_00' * psi_ss_00 + psi_ss_05' * psi_ss_05;
    end
end

% Reshape everything for plotting
spectrum_raw = reshape(spectrum_raw, num_pts_within_roundtrip, num_roundtrips, 6);

%% Solving Eqs. (19) and (25) of Methods section
odd_list = 1 : 2 : num_roundtrips;
if mod(num_roundtrips, 2) == 1, odd_list(end) = []; end
even_list = odd_list + 1;

spectrum_odd  = spectrum_raw(:, odd_list, :);
spectrum_even = spectrum_raw(:, even_list, :);
spectrum = zeros(size(spectrum_even));

for fast_time_index = 1 : num_pts_within_roundtrip
    for roundtrip_index = 1 : numel(odd_list)
        for psi_in_index = 1 : 6
            % Solving Eqs. (19) and (25) in the Methods section
            temp = 1 / (2 * (1 - exp(-4*gamma0))) * [1, -exp(-2*gamma0); -exp(-2*gamma0), 1] * ...
                   [spectrum_odd(fast_time_index, roundtrip_index, psi_in_index); ...
                    spectrum_even(fast_time_index, roundtrip_index, psi_in_index)];
            spectrum(fast_time_index, roundtrip_index, psi_in_index) = temp(2);
        end
    end
end

%% Calculate eigenstate
S = zeros(3, size(spectrum, 1)); % Rows: Sx, Sy, Sz

for Omegat_index = 1 : size(spectrum, 1)
    Sx = max(spectrum(Omegat_index, 1:round(end/2), 3)) ...
       - max(spectrum(Omegat_index, 1:round(end/2), 5));
    Sy = max(spectrum(Omegat_index, 1:round(end/2), 1)) ...
       - max(spectrum(Omegat_index, 1:round(end/2), 6));
    Sz = max(spectrum(Omegat_index, 1:round(end/2), 2)) ...
       - max(spectrum(Omegat_index, 1:round(end/2), 4));

    S(1, Omegat_index) = Sx / sqrt(Sx^2 + Sy^2 + Sz^2);
    S(2, Omegat_index) = Sy / sqrt(Sx^2 + Sy^2 + Sz^2);
    S(3, Omegat_index) = Sz / sqrt(Sx^2 + Sy^2 + Sz^2);
end

S(:, 1) = [-sign(var_phi); 0; 0]; % Remove numerical defect
S = transpose(S);

% In the simulation, we have kx = t going from 0 to 2pi. But in theory we
% have kx going from -pi to pi. So we need to re-arrange the sequence in
% which the simulation results are shown to match the theoretical results.
S = [S; S(1 : round(num_pts_within_roundtrip/2), :)];
S(1 : round(num_pts_within_roundtrip/2), :) = [];

%% Plot
[xx, yy, zz] = sphere(30); color_list = hsv(size(S, 1)); 
figure; lw = 3; lw2 = 3; sz = 80; ftsz = 16; lfsm = 1; ftnm = 'Arial'; fpos = [0.05 0.05 0.5 0.8];

surf(xx, yy, zz, 'FaceAlpha', 0.3, 'FaceColor', [1, 1, 1] * 0.95, 'EdgeColor', [1, 1, 1] * 0.7); hold on;

for loop_index = 1 : size(S, 1) - 1
    line([S(loop_index, 1); S(loop_index+1, 1)], ...
         [S(loop_index, 2); S(loop_index+1, 2)], ...
         [S(loop_index, 3); S(loop_index+1, 3)], 'LineWidth', lw2, 'Color', color_list(loop_index, :));
    hold on;
end
line([S(size(S, 1), 1); S(1, 1)], ...
     [S(size(S, 1), 2); S(1, 2)], ...
     [S(size(S, 1), 3); S(1, 3)], 'LineWidth', lw2, 'Color', color_list(loop_index, :)); hold on;

scatter3(S(:, 1), S(:, 2), S(:, 3), sz, color_list, "filled"); hold off;

lightangle(180, 20); lighting gouraud;
axis equal; axis([-1, 1, -1, 1, -1, 1]); view([150, 15]); grid off;
xticks([-1, 0, 1]); yticks([-1, 0, 1]); zticks([-1, 0, 1]);
xlabel('\sigma_x'); ylabel('\sigma_y'); zlabel('\sigma_z');
set(gca, 'fontname', ftnm, 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', lfsm, 'linewidth', lw, 'Layer', 'Top', 'Box', 'on');
set(gcf, 'unit', 'normalized', 'Position', fpos);


function y = transmittance(Omegat, g_TR, var_phi)
M = 3; delta_x = pi/2; delta_y = pi/2;
if mod(floor(Omegat/2/pi), 2) == 0
    y = [exp(1i * g_TR * cos(Omegat + delta_x)), 0;
         0, exp(1i * g_TR * cos(Omegat - delta_x))];
else
    y = [exp(1i * g_TR * cos(M * Omegat + delta_y + var_phi)), 0;
         0, exp(1i * g_TR * cos(M * Omegat - delta_y + var_phi))];
end
end