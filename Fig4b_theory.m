clearvars; close all; clc;
varphi = -pi/2; % Change as needed
Ax = pi/2 * pauli(3); Ay = pi/2 * pauli(1);

%% Define trajectory
kx_traj = linspace(-pi, pi, 251); kx_traj(end) = [];
ky_traj = 3 * kx_traj + varphi;

%% Calculate eigenstate
color_list = hsv(numel(kx_traj));
S = zeros(numel(kx_traj), 3);
for k_index = 1 : numel(kx_traj)
    kx = kx_traj(k_index); ky = ky_traj(k_index);
    H = cosm(kx * eye(2) - Ax) + cosm(ky * eye(2) - Ay);
    [eig_vec, eig_val] = eig(H, 'vector');
    [eig_val, eig_sort_index] = sort(eig_val, 'ascend');
    eig_vec = eig_vec(:, eig_sort_index);
    V_minus = eig_vec(:, 1); % Look at the lower band

    S(k_index, 1) = V_minus' * pauli(1) * V_minus;
    S(k_index, 2) = V_minus' * pauli(2) * V_minus;
    S(k_index, 3) = V_minus' * pauli(3) * V_minus;
end

%% Plot
[xx, yy, zz] = sphere(30);
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
xticks([-1, 0, 1]); xlabel('\sigma_x');
yticks([-1, 0, 1]); ylabel('\sigma_y');
zticks([-1, 0, 1]); zlabel('\sigma_z');
set(gca, 'fontname', ftnm, 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', lfsm, 'linewidth', lw, 'Layer', 'Top', 'Box', 'on');
set(gcf, 'unit', 'normalized', 'Position', fpos);