clearvars; clc; close all;

M = 3; g_TR = 0.9; gamma0 = 0.1;
U1 = [exp(1i*pi/4), 0; 0, exp(-1i*pi/4)]; 
U2 = [-1, -1; 1, -1] / sqrt(2);
var_phi = pi * 0.0;

phi = [0; 0.4; 0]; % Change as needed

num_pts_within_roundtrip = 250; num_roundtrips = 901;
spectrum = zeros(num_pts_within_roundtrip * num_roundtrips, 1);
Omegat_list = linspace(0, num_roundtrips, numel(spectrum) + 1) * 2*pi; Omegat_list(end) = [];

delta_omega_range = 0.15;
delta_omega_list = linspace(-1, 1, numel(spectrum) + 1) * delta_omega_range * 2*pi; delta_omega_list(end) = [];

T = @(x) transmittance(x, M, g_TR, var_phi, phi);
psi_in = [1; 1i]/sqrt(2);

for loop_index = 1 : numel(spectrum)
    Omegat = Omegat_list(loop_index);
    delta_omega = delta_omega_list(loop_index);
    P_00 = exp(1i * delta_omega) * exp(-gamma0);
    P_01 = P_00 * exp(1i * 2*pi/3);
    P_02 = P_00 * exp(1i * 4*pi/3);

    psi_ss_00 = inv(eye(2) - ...
                    U1 * T(Omegat) * U2 * ...
                    U1 * T(Omegat - 2*pi) * U2 * ...
                    U1 * T(Omegat - 4*pi) * U2 * P_00^3) * ...
                   (eye(2) + ...
                    U1 * T(Omegat) * U2 * P_00 + ...
                    U1 * T(Omegat) * U2 * U1 * T(Omegat - 2*pi) * U2 * P_00^2) * ...
                    psi_in;
    psi_ss_01 = inv(eye(2) - ...
                    U1 * T(Omegat) * U2 * ...
                    U1 * T(Omegat - 2*pi) * U2 * ...
                    U1 * T(Omegat - 4*pi) * U2 * P_01^3) * ...
                   (eye(2) + ...
                    U1 * T(Omegat) * U2 * P_01 + ...
                    U1 * T(Omegat) * U2 * U1 * T(Omegat - 2*pi) * U2 * P_01^2) * ...
                    psi_in;
    psi_ss_02 = inv(eye(2) - ...
                    U1 * T(Omegat) * U2 * ...
                    U1 * T(Omegat - 2*pi) * U2 * ...
                    U1 * T(Omegat - 4*pi) * U2 * P_02^3) * ...
                   (eye(2) + ...
                    U1 * T(Omegat) * U2 * P_02 + ...
                    U1 * T(Omegat) * U2 * U1 * T(Omegat - 2*pi) * U2 * P_02^2) * ...
                    psi_in;

    spectrum(loop_index) = psi_ss_00' * psi_ss_00 ...
                         + psi_ss_01' * psi_ss_01 ...
                         + psi_ss_02' * psi_ss_02;
end
% Time origin should start from half roundtrip in the plots
spectrum(1 : num_pts_within_roundtrip/2) = [];
spectrum(end - num_pts_within_roundtrip/2 + 1 : end) = [];
spectrum = reshape(spectrum, num_pts_within_roundtrip, num_roundtrips - 1);
spectrum(126, :) = (spectrum(125, :) + spectrum(127, :)) / 2; % Remove numerical defect
% Normalize
spectrum = spectrum / max(max(spectrum));

figure; lw = 2; sz = 20; ftsz = 22; lfsm = 1; ftnm = 'Arial';
fpos = [0.2 0.05 0.4 0.8];
[Omegat_plot, delta_omega_plot] = meshgrid(linspace(-0.5, 0.5, num_pts_within_roundtrip), linspace(-1, 1, num_roundtrips - 1) * delta_omega_range);
surf(Omegat_plot, delta_omega_plot, transpose(spectrum), 'EdgeColor', 'none');
colormap jet; view([0, 0, 1]); clim([min(min(spectrum)), max(max(spectrum))]); grid off;
xlim([-0.5, 0.5]); xticks([-0.5, 0, 0.5]);
ylim([-1, 1] * delta_omega_range); yticks([-1, 0, 1] * delta_omega_range);
xlabel('t / T_R'); ylabel('\delta\omega / \Omega_R');

set(gca, 'fontname', ftnm, 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', lfsm, 'linewidth', lw, 'Layer', 'Top', 'Box', 'on');
set(gcf, 'unit', 'normalized', 'Position', fpos);


function y = transmittance(Omegat, M, g_TR, var_phi, phi)
switch mod(floor(Omegat/2/pi), 3)
    case 0
        y = [exp(+1i * g_TR * (sin(M*Omegat + var_phi) + phi(1))), 0;
             0, exp(-1i * g_TR * (sin(M*Omegat + var_phi) + phi(1)))];
    case 1
        y = [exp(+1i * g_TR * phi(2)), 0;
             0, exp(-1i * g_TR * phi(2))];
    case 2
        y = [exp(+1i * g_TR * (sin(Omegat) + phi(3))), 0;
             0, exp(-1i * g_TR * (sin(Omegat) + phi(3)))];
    otherwise, error("Error!");
end
end