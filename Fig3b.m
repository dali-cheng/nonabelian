clearvars; close all; clc;
varphi = pi * 1.0; % Change as needed

Ax = pi/2 * pauli(3); Ay = pi/2 * pauli(1);

kx_list = linspace(-1, 1, 101) * pi;
E = zeros(numel(kx_list), 2);

for k_index = 1 : numel(kx_list)
    kx = kx_list(k_index);
    ky = 3 * kx + varphi;
    H = cosm(kx * eye(2) - Ax) + cosm(ky * eye(2) - Ay);
    E(k_index, :) = sort(eig(H));
end

figure; lw = 3; ftsz = 30;
plot(kx_list/pi, E(:, 1), 'LineWidth', lw, 'Color', 'k'); hold on;
plot(kx_list/pi, E(:, 2), 'LineWidth', lw, 'Color', 'k'); hold off;

xlim([-1, 1]); ylim([-2, 2]); xticks([-1, 0, 1]); yticks([-2, 0, 2]); 
xlabel('kx / pi'); ylabel('E');

set(gca, 'fontname', 'Arial', 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', lw, 'Layer', 'Top', 'box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.3 0.6]);