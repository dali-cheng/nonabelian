clearvars; clc; close all;
m0 = 0.5; % Change as needed

figure; lw = 2; sz = 20; ftsz = 22; lfsm = 1; ftnm = 'Arial';
fpos = [0.2 0.05 0.4 0.8];

[Omegat_plot, delta_omega_plot, spectrum1] = get_spectrum([1; 0], m0);
[~, ~, spectrum2] = get_spectrum([0; 1], m0);
spectrum = spectrum1 + spectrum2;

surf(Omegat_plot, delta_omega_plot, transpose(spectrum), 'EdgeColor', 'none');
colormap jet; view([0, 0, 1]); clim([min(min(spectrum)), max(max(spectrum))]); grid off;
xlim([-0.5, 0.5]); xticks([-0.5, 0, 0.5]);
ylim([-1, 1] * 0.015); yticks([-1, 0, 1] * 0.015);

xlabel('t / T_R'); ylabel('\delta\omega / \Omega_R');
set(gca, 'fontname', ftnm, 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', lfsm, 'linewidth', lw, 'Layer', 'Top', 'Box', 'on');
set(gcf, 'unit', 'normalized', 'Position', fpos);


function [Omegat_plot, delta_omega_plot, spectrum] = get_spectrum(psi_in, m0)
g_TR = 0.1; gamma0 = 0.005;
M = 3; varphi = pi * 0.0;
delta_omega_range = 0.015;

num_pts_within_roundtrip = 250; num_roundtrips = 901;
spectrum = zeros(num_pts_within_roundtrip * num_roundtrips, 1);
Omegat_list = linspace(0, num_roundtrips, numel(spectrum) + 1) * 2*pi; Omegat_list(end) = [];
delta_omega_list = linspace(-1, 1, numel(spectrum) + 1) * delta_omega_range * 2*pi; delta_omega_list(end) = [];
[Omegat_plot, delta_omega_plot] = meshgrid(linspace(-0.5, 0.5, num_pts_within_roundtrip), linspace(-1, 1, num_roundtrips - 1) * delta_omega_range);

T = @(x) transmittance(x, M, g_TR, varphi, m0);

for loop_index = 1 : numel(spectrum)
    Omegat = Omegat_list(loop_index);
    delta_omega = delta_omega_list(loop_index);
    P = exp(1i * delta_omega) * exp(-gamma0);

    psi_ss = inv(eye(2) - T(Omegat) * P) * psi_in;

    spectrum(loop_index) = psi_ss' * psi_ss;
end
% Time origin should start from half roundtrip in the plot
spectrum(1 : num_pts_within_roundtrip/2) = [];
spectrum(end - num_pts_within_roundtrip/2 + 1 : end) = [];
spectrum = reshape(spectrum, num_pts_within_roundtrip, num_roundtrips - 1);
% Normalize
spectrum = spectrum / max(max(spectrum));
end

function y = transmittance(Omegat, M, g_TR, var_phi, m0)
U1 = [1, -1; 1, 1]/sqrt(2); 
U2 = [-1i, -1i; 1, -1] / sqrt(2);

y1 = [exp(-1i * g_TR * (cos(M*Omegat + var_phi))), 0;
      0, exp(+1i * g_TR * (cos(M*Omegat + var_phi)))];
y2 = [exp(-1i * g_TR * (cos(Omegat))), 0;
      0, exp(+1i * g_TR * (cos(Omegat)))];
y3 = [exp(-1i * g_TR * m0), 0;
      0, exp(+1i * g_TR * m0)];
y4 = [exp(-1i * g_TR * (sin(M*Omegat + var_phi))), 0;
      0, exp(+1i * g_TR * (sin(M*Omegat + var_phi)))];
y5 = [exp(-1i * g_TR * (sin(Omegat))), 0;
      0, exp(+1i * g_TR * (sin(Omegat)))];

y = 0.2 * U1 * y1 * U1' + ...
    0.2 * U2 * y2 * U2' + ...
    0.2 * y3 + ...
    0.2 * y4 + ...
    0.2 * y5;
end