clearvars; clc; close all; 
filename = '3'; % Change as needed
load(['Data_Fig3d', filename, '.mat']);
% In these data files, 'FSR' is the free spectral range of the cavity,
% 'var_phi' is the modulation phase parameter, 'sweep_sig' is the signal
% that sweeps the laser frequency, and 'detector_sig' is the time-dependent
% transmission signal from the photodetector.

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
        omega_offset = -0.105; t_offset = 0.660;
    case '2'
        omega_offset = +0.227; t_offset = 0.810;
    case '3'
        omega_offset = -0.220; t_offset = 0.890;
end

t_offset_int = round(t_offset * size(bandstr, 2));
bandstr = [bandstr, bandstr(:, 1 : t_offset_int)]; bandstr(:, 1 : t_offset_int) = [];

%% Eq. (19) in the Methods section
delta_omega_00 = (delta_omega >= omega_offset - 0.25) & (delta_omega < omega_offset + 0.25);
delta_omega_05 = (delta_omega >= omega_offset + 0.25) & (delta_omega < omega_offset + 0.75);
bandstr_00 = bandstr(delta_omega_00, :);
bandstr_05 = bandstr(delta_omega_05, :);
delta_omega = delta_omega(delta_omega_00) - omega_offset;
if size(bandstr_00, 1) ~= size(bandstr_05, 1)
    bandstr_00 = bandstr_00(1 : min(size(bandstr_00, 1), size(bandstr_05, 1)), :);
    bandstr_05 = bandstr_05(1 : min(size(bandstr_00, 1), size(bandstr_05, 1)), :);
    delta_omega = delta_omega(1 : min(size(bandstr_00, 1), size(bandstr_05, 1)));
end
bandstr = bandstr_00 + bandstr_05;

%% Plot band structure
figure; lw = 3; ftsz = 30;
[fast_time_plot, delta_omega_plot] = meshgrid(fast_time, delta_omega);
surf(fast_time_plot, delta_omega_plot, bandstr, 'EdgeColor', 'none');
xlim([-0.5, 0.5]); ylim([-0.2, 0.2]); view([0, 0, 1]); colormap jet; grid off;
xticks([-0.5, 0, 0.5]); yticks([-0.2, 0, 0.2]);
xlabel('t / T_R'); ylabel('\delta\omega / \Omega_R');

set(gca, 'fontname', 'Arial', 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', lw, 'Layer', 'Top', 'box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.3 0.6]);