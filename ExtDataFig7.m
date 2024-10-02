clearvars; close all; clc;
Ax = pi/2 * pauli(3); Ay = pi/2 * pauli(1);
A0 = 0.9 * pauli(2); % Change as needed

num_kx = 51; num_ky = 51;
kx_list = linspace(-pi, pi, num_kx);
ky_list = linspace(-pi, pi, num_ky);
E = zeros(num_kx, num_ky, 2);

for kx_index = 1 : num_kx
    for ky_index = 1 : num_ky
        kx = kx_list(kx_index); ky = ky_list(ky_index);
        H = cosm(kx * eye(2) - Ax) + cosm(ky * eye(2) - Ay) + A0;
        H = 2/3 * H;
        % Because we are using the three-roundtrip interleaving modulation
        % scheme, we need to rescale the Hamiltonian by 2/3 to make it
        % comparable to our results with two-roundtrip interleaving scheme.
        H = (H + H') / 2; % Force H to be Hermitian
        E(kx_index, ky_index, :) = sort(eig(H));
    end
end

figure; lw = 2; ftsz = 24;
[kx_plot, ky_plot] = meshgrid(kx_list, ky_list);
mesh(kx_plot/pi, ky_plot/pi, transpose(squeeze(E(:, :, 1))), ...
    'EdgeColor', 'k'); hold on;
mesh(kx_plot/pi, ky_plot/pi, transpose(squeeze(E(:, :, 2))), ...
    'EdgeColor', 'k'); hold off;

xlabel('kx / pi'); ylabel('ky / pi'); zlabel('E');
xticks([-1, 0, 1]); yticks([-1, 0, 1]); zticks([-1, 0, 1] * 2);
axis([-1, 1, -1, 1, -2, 2]); view([156, 10]); grid off;
set(gca, 'fontname', 'Arial', 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', lw, 'Layer', 'Top', 'box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.4 0.6]);