clearvars; close all; clc;
N = 50; % Number of lattice sites
Ax = (pi/2) * pauli(3); Ay = (pi/2) * pauli(1);
A0 = 0.1 * pauli('y') + 0.75 * pauli('z');

num_ky = 101;
ky_list = linspace(-pi, pi, num_ky);
E = zeros(4*N, num_ky);

for ky_index = 1 : num_ky
    ky = ky_list(ky_index);
    H = zeros(4*N, 4*N);
    
    % PBC along (11) direction
    for x = 1 : 2*N-1
        H(2*x-1 : 2*x, 2*x+1 : 2*x+2) = (expm(+1i * Ax) + expm(-1i * ky * eye(2) - 1i * Ay))/2;
        H(2*x+1 : 2*x+2, 2*x-1 : 2*x) = (expm(-1i * Ax) + expm(+1i * ky * eye(2) + 1i * Ay))/2;
    end 
    H(1 : 2, 4*N-1 : 4*N) = (expm(-1i * Ax) + expm(+1i * ky * eye(2) + 1i * Ay))/2;
    H(4*N-1 : 4*N, 1 : 2) = (expm(+1i * Ax) + expm(-1i * ky * eye(2) - 1i * Ay))/2;

    for x = 1 : 2*N
        H(2*x-1 : 2*x, 2*x-1 : 2*x) = ternary(x<=N, A0, -A0);
    end % Scalar potential

    H = 2/3 * H; % For consistency with two-roundtrip calculations
    H = (H + H')/2; % Forcing Hermiticity
    E(:, ky_index) = sort(eig(H, 'vector'));
end

%% Delete the edge states on the unwanted interface
% In this calculation we are effectively solving a system with two parallel
% interfaces. To show the band structure of a semi-infinite lattice, we
% delete the edge state that is localized on one of the interfaces.
E([99, 101], [1:13, 39:50]) = nan;
E([100, 102], [52:63, 89:101]) = nan;
E([99, 100], 14:38) = nan;
E([101, 102], 64:88) = nan;

%% Plotting bandstr
figure; lw = 2; sz = 6; ftsz = 24;
for E_index = 1 : 4*N
    scatter(ky_list/pi, E(E_index, :), sz, 'k', 'filled'); hold on;
end
xticks([-1, 0, 1] * 1.0); xlim([-1, 1] * 1.0);
yticks([-1, 0, 1] * 2); ylim([-1, 1] * 2);

xlabel('kp / (pi/sqrt(2))'); ylabel('E');

set(gca, 'fontname', 'Arial', 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', lw, 'Layer', 'Top', 'box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.05 0.05 0.3 0.8]);

%% Plot eigenstates
ky = 0.1 * pi; E_index = 101;
H = zeros(4*N, 4*N);
% PBC along (11) direction
for x = 1 : 2*N-1
    H(2*x-1 : 2*x, 2*x+1 : 2*x+2) = (expm(+1i * Ax) + expm(-1i * ky * eye(2) - 1i * Ay))/2;
    H(2*x+1 : 2*x+2, 2*x-1 : 2*x) = (expm(-1i * Ax) + expm(+1i * ky * eye(2) + 1i * Ay))/2;
end
H(1 : 2, 4*N-1 : 4*N) = (expm(-1i * Ax) + expm(+1i * ky * eye(2) + 1i * Ay))/2;
H(4*N-1 : 4*N, 1 : 2) = (expm(+1i * Ax) + expm(-1i * ky * eye(2) - 1i * Ay))/2;

for x = 1 : 2*N
    H(2*x-1 : 2*x, 2*x-1 : 2*x) = ternary(x<=N, A0, -A0);
end % Scalar potential

H = 2/3 * H; % For consistency with two-roundtrip calculations
H = (H + H')/2; % Forcing Hermiticity

[V, E] = eig(H, 'vector');
eigenstate_up = V(1 : 2 : end, E_index);
eigenstate_down = V(2 : 2 : end, E_index);

figure; mksz = 5;
stem(1:2*N, abs(eigenstate_up), "filled", LineWidth=lw, Color='k', MarkerSize=mksz); hold on;
stem(1:2*N, -abs(eigenstate_down), "filled", LineWidth=lw, Color='k', MarkerSize=mksz); hold off;
xticks([0, 0.5, 1] * 2*N); xlim([0, 2*N]);
yticks([-1, 0, 1]); ylim([-1, 1]); yticklabels([1, 0, 1]);

xlabel('Lattice site'); ylabel('|Eigenstate|');

set(gca, 'fontname', 'Arial', 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', 1, 'linewidth', lw, 'Layer', 'Top', 'box', 'on');
set(gcf, 'unit', 'normalized', 'Position', [0.05 0.05 0.6 0.5]);