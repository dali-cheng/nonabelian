clearvars; close all; clc;

Ax = pi/2 * pauli(3); Ay = pi/2 * pauli(1);

num_kx = 51; num_ky = 51;
kx_list = linspace(-pi, pi, num_kx);
ky_list = linspace(-pi, pi, num_ky);
E = zeros(num_kx, num_ky, 2);

for kx_index = 1 : num_kx
    for ky_index = 1 : num_ky
        kx = kx_list(kx_index); ky = ky_list(ky_index);
        H = cosm(kx * eye(2) - Ax) + cosm(ky * eye(2) - Ay);
        E(kx_index, ky_index, :) = sort(eig(H));
    end
end

[kx_plot, ky_plot] = meshgrid(kx_list, ky_list);
figure; lw = 3; ftsz = 30;
mesh(kx_plot/pi, ky_plot/pi, transpose(squeeze(E(:, :, 1))), ...
     'EdgeColor', 'k'); hold on;
mesh(kx_plot/pi, ky_plot/pi, transpose(squeeze(E(:, :, 2))), ...
     'EdgeColor', 'k'); hold off;

xlabel('kx / pi'); xticks([-1, 0, 1]); 
ylabel('ky / pi'); yticks([-1, 0, 1]); 
zlabel('E'); zticks([-2, 0, 2]); 
axis([-1, 1, -1, 1, -2, 2]); view([156, 10]);
grid off;
set(gca, 'fontname', 'Arial', 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', lw, 'Layer', 'Top', 'box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.2 0.05 0.5 0.6]);