clearvars; clc; close all;
filename = '1'; load(['Data_Fig4_', filename, '.mat']); % Change as needed

g_TR = 0.4; gamma = 0.1; 
switch filename
    case '1'
        sync_offset = 0.62; t_offset = 175; omega_offset = [0.66, 1.685, -3.355, -2.332, -1.33, -0.323]; 
        % Sync osc data acquisition with modulation
    case '2'
        sync_offset = 1.01; t_offset = 8; omega_offset = [-4.792, -3.76, -2.79, -1.757, -0.754, 0.252]; 
end

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

%% Calculate eigenstate
% Look at a narrower window to avoid noise outside
energy_upper_bound = 0.48;
energy_lower_bound = 0.34;

S = zeros(3, numel(fast_time));
% Rows: Sx, Sy, Sz
for fast_time_index = 1 : numel(fast_time)
    % Sx = X - Y = State 3 - State 5
    % Sy = L - R = State 1 - State 6
    % Sz = H - V = State 2 - State 4
    Sx = max(bandstr{3}(round(end * energy_lower_bound) : round(end * energy_upper_bound), fast_time_index)) ...
        - max(bandstr{5}(round(end * energy_lower_bound) : round(end * energy_upper_bound), fast_time_index));
    Sy = max(bandstr{1}(round(end * energy_lower_bound) : round(end * energy_upper_bound), fast_time_index)) ...
        - max(bandstr{6}(round(end * energy_lower_bound) : round(end * energy_upper_bound), fast_time_index));
    Sz = max(bandstr{2}(round(end * energy_lower_bound) : round(end * energy_upper_bound), fast_time_index)) ...
        - max(bandstr{4}(round(end * energy_lower_bound) : round(end * energy_upper_bound), fast_time_index));

    S(1, fast_time_index) = Sx / sqrt(Sx^2 + Sy^2 + Sz^2);
    S(2, fast_time_index) = Sy / sqrt(Sx^2 + Sy^2 + Sz^2);
    S(3, fast_time_index) = Sz / sqrt(Sx^2 + Sy^2 + Sz^2);
end

%% Post-processing on the trajectory
S = transpose(S);

% Apply global rotation. This additional SO(3) rotation comes from the SM
% pigtail of the PSG. 
U_rot = rotm(so3([0.03, 0.75, 0.15] * pi, "eul", "YZX"));
for loop_index = 1 : size(S, 1)
    S(loop_index, :) = transpose(U_rot * transpose(S(loop_index, :)));
end

S = [S; S(1 : t_offset, :)];
S(1 : t_offset, :) = [];

save(['Data_ExtDataFig6_', filename, '.mat'], 'S');

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

scatter3(S(:, 1), S(:, 2), S(:, 3), sz, color_list, "filled"); hold on;
if filename == '1'
    sz = 1440; scatter3(S(1, 1), S(1, 2), S(1, 3), sz, "k", "filled", "pentagram"); hold off;
end

lightangle(180, 20); lighting gouraud;
axis equal; axis([-1, 1, -1, 1, -1, 1]); view([150, 15]); grid off;
xticks([-1, 0, 1]); yticks([-1, 0, 1]); zticks([-1, 0, 1]);
xlabel('\sigma_x'); ylabel('\sigma_y'); zlabel('\sigma_z');
set(gca, 'fontname', ftnm, 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', lfsm, 'linewidth', lw, 'Layer', 'Top', 'Box', 'on');
set(gcf, 'unit', 'normalized', 'Position', fpos);