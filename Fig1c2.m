clearvars; close all; clc;
grid_num = 50;

Ax = pi/2 * pauli(3); Ay = pi/2 * pauli(1);

t_list = 0 : round(grid_num * (2 + sqrt(2))) - 1;
E = zeros(numel(t_list), 2);

for t_index = 1 : numel(t_list)
    t = t_list(t_index);
    if t < grid_num
        kx = t / grid_num * pi; ky = 0; % Gamma to X
    else if t < 2 * grid_num
        kx = pi; ky = (t - grid_num) / grid_num * pi; % X to M
    else
        kx = (1 - (t - 2 * grid_num) / grid_num / sqrt(2)) * pi;
        ky = kx; % M to Gamma
    end
    end
    H = cosm(kx * eye(2) - Ax) + cosm(ky * eye(2) - Ay);
    E(t_index, :) = sort(eig(H));
end

figure; lw = 3; ftsz = 30;
plot(t_list, E(:, 1), 'LineWidth', lw, 'Color', 'k'); hold on;
plot(t_list, E(:, 2), 'LineWidth', lw, 'Color', 'k'); hold on;
line([50, 50], [-4, 4], 'LineWidth', lw, 'Color', 'k'); hold on;
line([100, 100], [-4, 4], 'LineWidth', lw, 'Color', 'k'); hold off;

xlim([0, 170]); ylim([-2, 2]); yticks([-2, 0, 2]); xticks([]); ylabel('E');

set(gca, 'fontname', 'Arial', 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', lw, 'Layer', 'Top', 'box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.35 0.6]);