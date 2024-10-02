clearvars; close all; clc;
varphi = 0.0;
Ax = pi/2 * pauli(3); Ay = pi/2 * pauli(1);
A0 = 0.9 * pauli(2); % Change as needed

kx_list = linspace(-1, 1, 301) * pi;
E = zeros(numel(kx_list), 2);

for k_index = 1 : numel(kx_list)
    kx = kx_list(k_index);
    ky = 3 * kx + varphi;
    H = cosm(kx * eye(2) - Ax) + cosm(ky * eye(2) - Ay) + A0;
    H = 2/3 * H; % For consistency with calculations in the two-roundtrip interleaving modulation scheme
    H = (H + H') / 2; % Forcing Hermiticity
    E(k_index, :) = sort(eig(H));
end

figure; lw = 2; sz = 20; ftsz = 22; lfsm = 1; ftnm = 'Arial';
fpos = [0.2 0.05 0.4 0.8];
plot(kx_list/pi, E(:, 1), 'LineWidth', lw, 'Color', 'k'); hold on;
plot(kx_list/pi, E(:, 2), 'LineWidth', lw, 'Color', 'k'); hold off;

xlim([-1, 1]); ylim([-1, 1] * 2); xticks([-1, 0, 1]); yticks([-1, 0, 1] * 2); 
xlabel('kx / pi'); ylabel('E');
set(gca, 'fontname', ftnm, 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', lfsm, 'linewidth', lw, 'Layer', 'Top', 'box', 'on');
set(gcf, 'unit', 'normalized', 'Position', fpos);