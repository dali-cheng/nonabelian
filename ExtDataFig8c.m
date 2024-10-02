clearvars; clc; close all; 
filename = '3'; % Change as needed
load(['Data_ExtDataFig8c', filename, '.mat']);
% In these data files, 'FSR' is the free spectral range of the cavity,
% 'phi' is [phi_x, phi_y, phi_z] where the scalar potential is A0 = phi dot
% sigma, 'sweep_sig' is the signal that sweeps the laser frequency, and 
% 'detector_sig' is the time-dependent transmission signal from the 
% photodetector.

%% Segment transmission spectrum to create band structure
chunk_size = detector_sig.desc.fs / FSR;
bandstr = zeros(floor(numel(detector_sig.y) / chunk_size), round(chunk_size));
for row_index = 1 : size(bandstr, 1)
    bandstr(row_index, :) = detector_sig.y(...
       (round(row_index * chunk_size) - round(chunk_size) + 1) : round(row_index * chunk_size));
end

%% Calculate frequency (energy) range
fast_time = linspace(-0.5, 0.5, size(bandstr, 2) + 1); fast_time(end) = [];
delta_omega_range = range(sweep_sig.y) * 1000 * 4/30 / (FSR/1e6);
delta_omega = linspace(-delta_omega_range/2, delta_omega_range/2, size(bandstr, 1));

%% Find time and energy origins
switch filename
    case '1'
        omega_offset = -0.525; t_offset = 0.520;
    case '2'
        omega_offset = -0.195; t_offset = 0.570;
    case '3'
        omega_offset = -0.740; t_offset = 0.635;
end

t_offset_int = round(t_offset * size(bandstr, 2));
bandstr = [bandstr, bandstr(:, 1 : t_offset_int)]; bandstr(:, 1 : t_offset_int) = [];

%% Genaralization of Eq. (19) in the Methods section
delta_omega_00 = (delta_omega >= omega_offset - 1/6) & (delta_omega < omega_offset + 1/6);
delta_omega_01 = (delta_omega >= omega_offset + 1/6) & (delta_omega < omega_offset + 1/2);
delta_omega_02 = (delta_omega >= omega_offset + 1/2) & (delta_omega < omega_offset + 5/6);
bandstr_00 = bandstr(delta_omega_00, :);
bandstr_01 = bandstr(delta_omega_01, :);
bandstr_02 = bandstr(delta_omega_02, :);
delta_omega_high_reso = delta_omega(delta_omega_00) - omega_offset;
% In case the dimensions do not match:
new_bandstr_dimension = min([size(bandstr_00, 1), size(bandstr_01, 1), size(bandstr_02, 1)]);
bandstr_00 = bandstr_00(1 : new_bandstr_dimension, :);
bandstr_01 = bandstr_01(1 : new_bandstr_dimension, :);
bandstr_02 = bandstr_02(1 : new_bandstr_dimension, :);
delta_omega_high_reso = delta_omega_high_reso(1 : new_bandstr_dimension);

bandstr = bandstr_00 + bandstr_01 + bandstr_02;

%% Plot band structure
[fast_time_plot, delta_omega_plot] = meshgrid(fast_time, delta_omega_high_reso);
figure; lw = 2; sz = 20; ftsz = 22; lfsm = 1; ftnm = 'Arial';
fpos = [0.2 0.05 0.4 0.8];
surf(fast_time_plot, delta_omega_plot, bandstr, 'EdgeColor', 'none'); grid off;
xlim([-0.5, 0.5]); ylim([-1, 1] * 0.15); 
xticks([-0.5, 0, 0.5]); yticks([-1, 0, 1] * 0.15); 
view([0, 0, 1]); colormap jet;
xlabel('t / T_R'); ylabel('\delta\omega / \Omega_R'); 

set(gca, 'fontname', ftnm, 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', lfsm, 'linewidth', lw, 'Layer', 'Top', 'Box', 'on');
set(gcf, 'unit', 'normalized', 'Position', fpos);