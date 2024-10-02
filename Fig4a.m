clearvars; clc; close all;
load('Data_Fig4_1.mat');

g_TR = 0.4; gamma = 0.1; 
sync_offset = 0.62; omega_offset = [0.66, 1.685, -3.355, -2.332, -1.33, -0.323]; 
% Sync osc data acquisition with modulation

%% Generate band structure
chunk_size = detector_sig.desc.fs / FSR;
detector_sig.y(1 : round(sync_offset * chunk_size)) = [];
bandstr_raw = zeros(floor(numel(detector_sig.y) / chunk_size), round(chunk_size));
for row_index = 1 : size(bandstr_raw, 1)
    bandstr_raw(row_index, :) = detector_sig.y(...
       (round(row_index * chunk_size) - round(chunk_size) + 1) : round(row_index * chunk_size));
end
bandstr_raw = bandstr_raw - 0.012; 
% Deduct the background, the background noise from pre-amp

%% State segmentation
fast_time = linspace(0, 1, size(bandstr_raw, 2) + 1); fast_time(end) = [];
% Time origin up to a constant shift
delta_omega_range = 10;
delta_omega = linspace(-delta_omega_range/2, delta_omega_range/2, size(bandstr_raw, 1));

bandstr_proj = {}; delta_omega_proj = {};
for state_index = 1 : 6
    delta_omega_00 = (delta_omega >= omega_offset(state_index) + 0.0) & (delta_omega < omega_offset(state_index) + 0.5);
    delta_omega_05 = (delta_omega >= omega_offset(state_index) + 0.5) & (delta_omega < omega_offset(state_index) + 1.0);
    bandstr_00 = bandstr_raw(delta_omega_00, :);
    bandstr_05 = bandstr_raw(delta_omega_05, :);
    delta_omega_proj{state_index} = delta_omega(delta_omega_00) - omega_offset(state_index) - 0.25;
    if size(bandstr_00, 1) ~= size(bandstr_05, 1)
        bandstr_00 = bandstr_00(1 : min(size(bandstr_00, 1), size(bandstr_05, 1)), :);
        bandstr_05 = bandstr_05(1 : min(size(bandstr_00, 1), size(bandstr_05, 1)), :);
        delta_omega_proj{state_index} = delta_omega_proj{state_index}(1 : min([size(bandstr_00, 1), size(bandstr_05, 1)]));
    end
    bandstr_proj{state_index} = bandstr_00 + bandstr_05; % Eq. (19) of Methods section
end

% Normalize. See section "Calibration of transmittance of the free-space 
% beam splitter" for explanations. These numbers are measured separately.
bandstr_proj{1} = bandstr_proj{1} / 1.000;
bandstr_proj{2} = bandstr_proj{2} / 1.211;
bandstr_proj{3} = bandstr_proj{3} / 0.911;
bandstr_proj{4} = bandstr_proj{4} / 0.621;
bandstr_proj{5} = bandstr_proj{5} / 0.941;
bandstr_proj{6} = bandstr_proj{6} / 0.836;

%% Solving Eqs. (19) and (25) of Methods section
num_roundtrips = min([size(bandstr_proj{1}, 1), size(bandstr_proj{2}, 1), ...
                      size(bandstr_proj{3}, 1), size(bandstr_proj{4}, 1), ...
                      size(bandstr_proj{5}, 1), size(bandstr_proj{6}, 1)]);
odd_list = 1 : 2 : num_roundtrips;
if mod(num_roundtrips, 2) == 1, odd_list(end) = []; end
even_list = odd_list + 1;

for state_index = 1 : 6
    bandstr_odd{state_index}  = bandstr_proj{state_index}(odd_list, :);
    bandstr_even{state_index} = bandstr_proj{state_index}(even_list, :);
    bandstr{state_index}   = zeros(size(bandstr_even{state_index}));

    for fast_time_index = 1 : size(bandstr_proj{1}, 2)
        for roundtrip_index = 1 : numel(odd_list)
            temp = 1 / (2 * (1 - exp(-4*gamma))) * [1, -exp(-2*gamma); -exp(-2*gamma), 1] * ...
                   [bandstr_odd{state_index}(roundtrip_index, fast_time_index); ...
                    bandstr_even{state_index}(roundtrip_index, fast_time_index)];
            bandstr{state_index}(roundtrip_index, fast_time_index) = temp(2);
        end
    end
end

for state_index = 1 : 6
    delta_omega_proj{state_index} = delta_omega_proj{state_index}(odd_list);
end

%% Plot band structure
lw = 3; sz = 1920; ftsz = 16; lfsm = 1; ftnm = 'Arial'; fpos = [0.05 0.05 0.6 0.8];
for state_index = 1 : 6
    [fast_time_plot, delta_omega_proj_plot] = meshgrid(fast_time, delta_omega_proj{state_index});

    figure; surf(fast_time_plot, delta_omega_proj_plot, bandstr{state_index}, 'EdgeColor', 'none'); hold on;
    
    xlabel('t / T_R'); xticklabels([]); ylabel('\delta\omega / \Omega_R');
    xlim([0, 1]); xticks([0, 0.5, 1]); ylim([-0.1, 0.1]); yticks([-0.1, 0, 0.1]);
    view([0, 0, 1]); colormap jet; clim([0, 0.03]); grid off;

    ax(state_index) = gca;
    set(gca, 'fontname', ftnm, 'fontsize', ftsz, 'fontweight', 'normal', ...
        'labelfontsizemultiplier', lfsm, 'linewidth', lw, 'Layer', 'Top', 'Box', 'on');
    set(gcf, 'unit', 'normalized', 'Position', fpos);
end

%% Find peak intensity
fast_time_index = 176;
% Look at a narrower window to avoid noise outside
energy_upper_bound = 0.48;
energy_lower_bound = 0.34;

[~, I1] = max(bandstr{1}(round(end * energy_lower_bound) : round(end * energy_upper_bound), fast_time_index));
[~, I2] = max(bandstr{2}(round(end * energy_lower_bound) : round(end * energy_upper_bound), fast_time_index));
[~, I3] = max(bandstr{3}(round(end * energy_lower_bound) : round(end * energy_upper_bound), fast_time_index));
[~, I4] = max(bandstr{4}(round(end * energy_lower_bound) : round(end * energy_upper_bound), fast_time_index));
[~, I5] = max(bandstr{5}(round(end * energy_lower_bound) : round(end * energy_upper_bound), fast_time_index));
[~, I6] = max(bandstr{6}(round(end * energy_lower_bound) : round(end * energy_upper_bound), fast_time_index));

I1 = I1 - 1 + round(size(bandstr{1}, 1) * energy_lower_bound);
I2 = I2 - 1 + round(size(bandstr{2}, 1) * energy_lower_bound);
I3 = I3 - 1 + round(size(bandstr{3}, 1) * energy_lower_bound);
I4 = I4 - 1 + round(size(bandstr{4}, 1) * energy_lower_bound);
I5 = I5 - 1 + round(size(bandstr{5}, 1) * energy_lower_bound);
I6 = I6 - 1 + round(size(bandstr{6}, 1) * energy_lower_bound);

hold(ax(1), 'on'); scatter3(ax(1), fast_time(fast_time_index), delta_omega_proj{1}(I1), 10, sz, 'k', 'filled', "pentagram");
hold(ax(2), 'on'); scatter3(ax(2), fast_time(fast_time_index), delta_omega_proj{2}(I2), 10, sz, 'k', 'filled', "pentagram");
hold(ax(3), 'on'); scatter3(ax(3), fast_time(fast_time_index), delta_omega_proj{3}(I3), 10, sz, 'k', 'filled', "pentagram");
hold(ax(4), 'on'); scatter3(ax(4), fast_time(fast_time_index), delta_omega_proj{4}(I4), 10, sz, 'k', 'filled', "pentagram");
hold(ax(5), 'on'); scatter3(ax(5), fast_time(fast_time_index), delta_omega_proj{5}(I5), 10, sz, 'k', 'filled', "pentagram");
hold(ax(6), 'on'); scatter3(ax(6), fast_time(fast_time_index), delta_omega_proj{6}(I6), 10, sz, 'k', 'filled', "pentagram");

lw = 8;
hold(ax(1), 'on'); line(ax(1), [1, 1] * fast_time(fast_time_index), [-0.1, 0.1], [1, 1], 'LineWidth', lw, 'Color', 'w', 'LineStyle', ':');
hold(ax(2), 'on'); line(ax(2), [1, 1] * fast_time(fast_time_index), [-0.1, 0.1], [1, 1], 'LineWidth', lw, 'Color', 'w', 'LineStyle', ':');
hold(ax(3), 'on'); line(ax(3), [1, 1] * fast_time(fast_time_index), [-0.1, 0.1], [1, 1], 'LineWidth', lw, 'Color', 'w', 'LineStyle', ':');
hold(ax(4), 'on'); line(ax(4), [1, 1] * fast_time(fast_time_index), [-0.1, 0.1], [1, 1], 'LineWidth', lw, 'Color', 'w', 'LineStyle', ':');
hold(ax(5), 'on'); line(ax(5), [1, 1] * fast_time(fast_time_index), [-0.1, 0.1], [1, 1], 'LineWidth', lw, 'Color', 'w', 'LineStyle', ':');
hold(ax(6), 'on'); line(ax(6), [1, 1] * fast_time(fast_time_index), [-0.1, 0.1], [1, 1], 'LineWidth', lw, 'Color', 'w', 'LineStyle', ':');