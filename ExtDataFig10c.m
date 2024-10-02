clearvars; close all; clc;
varphi = pi * 0.0; 

%% Parameter definition
t1 = 1; t2 = 1; t3 = 1;
Ux = (t2/2) * pauli(2) - (1i*t3/2) * pauli(3);
Uy = (t1/2) * pauli(1) - (1i*t3/2) * pauli(3);
m0 = 0.0; M = m0 * pauli(3); % Change as needed

kx_list = linspace(-1, 1, 101) * pi;
E = zeros(numel(kx_list), 2);

%% Solve bandstr
for k_index = 1 : numel(kx_list)
    kx = kx_list(k_index);
    ky = 3 * kx + varphi;
    H = exp(1i * kx) * Ux + exp(1i * ky) * Uy;
    H = H + H';
    H = H + M;
    E(k_index, :) = sort(eig(H));
end

%% Plotting
figure; lw = 3; ftsz = 30;
plot(kx_list/pi, E(:, 1), 'LineWidth', lw, 'Color', 'k'); hold on;
plot(kx_list/pi, E(:, 2), 'LineWidth', lw, 'Color', 'k'); hold off;

xlim([-1, 1]); xticks([-1, 0, 1]); 
ylim([-1, 1] * 4); yticks([-1, 0, 1] * 4); 
xlabel('kx / pi'); ylabel('E');

set(gca, 'fontname', 'Arial', 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', lw, 'Layer', 'Top', 'box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.3 0.6]);