clearvars; clc; close all; 

M = 3; g_TR = 0.9; gamma0 = 0.2;
U1 = [1, -1; 1, 1] / sqrt(2); U2 = [1, 0; 0, -1];
varphi = pi * 0.0; % Change as needed

num_pts_within_roundtrip = 250; num_roundtrips = 801;
spectrum = zeros(num_pts_within_roundtrip * num_roundtrips, 1);
Omegat_list = linspace(0, num_roundtrips, numel(spectrum) + 1) * 2*pi; Omegat_list(end) = [];
delta_omega_list = linspace(-1, 1, numel(spectrum) + 1) * 0.2 * 2*pi; delta_omega_list(end) = [];

T = @(x) transmittance(x, M, g_TR, varphi);
psi_in = [0; 1];

for loop_index = 1 : numel(spectrum)
    Omegat = Omegat_list(loop_index);
    delta_omega = delta_omega_list(loop_index);
    P_00 =  exp(1i * delta_omega) * exp(-gamma0);
    P_05 = -exp(1i * delta_omega) * exp(-gamma0);

    psi_ss_00 = inv(eye(2) - U1 * T(Omegat) * U2 * U1 * T(Omegat - 2*pi) * U2 * P_00^2) * ...
        (eye(2) + P_00 * U1 * T(Omegat) * U2) * psi_in;
    psi_ss_05 = inv(eye(2) - U1 * T(Omegat) * U2 * U1 * T(Omegat - 2*pi) * U2 * P_05^2) * ...
        (eye(2) + P_05 * U1 * T(Omegat) * U2) * psi_in;

    spectrum(loop_index) = psi_ss_00' * psi_ss_00 + psi_ss_05' * psi_ss_05;
end
% The time origin should start from half roundtrip in the plot
spectrum(1 : num_pts_within_roundtrip/2) = [];
spectrum(end - num_pts_within_roundtrip/2 + 1 : end) = [];
spectrum = reshape(spectrum, num_pts_within_roundtrip, num_roundtrips - 1);
% Remove a numerical defect
if varphi == 0.5 * pi
    spectrum(126, :) = spectrum(126 - 1, :);
end
% Normalize
spectrum = spectrum / max(max(spectrum));

figure; lw = 3; ftsz = 30;
[Omegat_plot, delta_omega_plot] = meshgrid(linspace(-0.5, 0.5, num_pts_within_roundtrip), linspace(-1, 1, num_roundtrips - 1) * 0.2);
surf(Omegat_plot, delta_omega_plot, transpose(spectrum), 'EdgeColor', 'none');
colormap jet; view([0, 0, 1]); clim([min(min(spectrum)), max(max(spectrum))]); grid off;

xlim([-0.5, 0.5]); ylim([-0.2, 0.2]); xticks([-0.5, 0, 0.5]); yticks([-0.2, 0, 0.2]);
xlabel('t / T_R'); ylabel('\delta\omega / \Omega_R');
set(gca, 'fontname', 'Arial', 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', lw, 'Layer', 'Top', 'box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.3 0.6]);


function y = transmittance(Omegat, M, g_TR, var_phi)
delta_x = pi/2; delta_y = pi/2;
if mod(floor(Omegat/2/pi), 2) == 0
    y = [exp(1i * g_TR * cos(Omegat + delta_x)), 0;
         0, exp(1i * g_TR * cos(Omegat - delta_x))];
else
    y = [exp(1i * g_TR * cos(M * Omegat + delta_y + var_phi)), 0;
         0, exp(1i * g_TR * cos(M * Omegat - delta_y + var_phi))];
end
end