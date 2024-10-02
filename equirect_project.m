function [xx, yy] = equirect_project(v, azi_ctr)
%EQUIRECT_PROJECT Perform and plot the equirectangular projection
%   Input v: An n*3 matrix, each row represents a datapoint on the 3D
%   sphere, and column 1 (2, 3) represents the x (y, z) coordinates.
%
%   Input azi_ctr: the center position of azimuthal angle for doing the 
%   projection. Default is 0.
%   
%   Output xx and yy: After the equirectangular projection, we use xx and
%   yy data to illustrate the datapoints in a 2D planar plot.

if nargin < 2, azi_ctr = 0; end

n = size(v, 1); 
xx = zeros(n, 1); yy = xx;

for loop_index = 1 : n
    P = v(loop_index, :); P = P / norm(P);
    % Calculate polar angle
    yy(loop_index) = acos(P(3));
    % Calculate azimuthal angle
    xx(loop_index) = atan2(P(2), P(1));
end

if azi_ctr > 0
    xx(xx < -pi + azi_ctr) = xx(xx < -pi + azi_ctr) + 2 * pi;
else
    xx(xx >  pi + azi_ctr) = xx(xx >  pi + azi_ctr) - 2 * pi;
end

% Plot the 2D planar projection
color_list = hsv(n); 
figure; lw = 2; sz = 30; ftsz = 16; lfsm = 1; ftnm = 'Arial'; fpos = [0.05 0.05 0.5 0.8];

for loop_index = 1 : n-1
    if abs(xx(loop_index) - xx(loop_index+1))/pi < 1
        line([xx(loop_index); xx(loop_index+1)]/pi, ...
             [yy(loop_index); yy(loop_index+1)]/pi, 'LineWidth', lw, 'Color', color_list(loop_index, :));
        hold on;
    else
        if xx(loop_index) < xx(loop_index+1)
            line([xx(loop_index); xx(loop_index+1) - 2*pi]/pi, ...
                 [yy(loop_index); yy(loop_index+1)]/pi, 'LineWidth', lw, 'Color', color_list(loop_index, :));
            hold on;
            line([xx(loop_index) + 2*pi; xx(loop_index+1)]/pi, ...
                 [yy(loop_index); yy(loop_index+1)]/pi, 'LineWidth', lw, 'Color', color_list(loop_index, :));
            hold on;
        else
            line([xx(loop_index); xx(loop_index+1) + 2*pi]/pi, ...
                 [yy(loop_index); yy(loop_index+1)]/pi, 'LineWidth', lw, 'Color', color_list(loop_index, :));
            hold on;
            line([xx(loop_index) - 2*pi; xx(loop_index+1)]/pi, ...
                 [yy(loop_index); yy(loop_index+1)]/pi, 'LineWidth', lw, 'Color', color_list(loop_index, :));
            hold on;
        end
    end
end
if abs(xx(1) - xx(n))/pi < 1
    line([xx(1); xx(n)]/pi, ...
         [yy(1); yy(n)]/pi, 'LineWidth', lw, 'Color', color_list(loop_index, :)); hold on;
end
scatter(xx/pi, yy/pi, sz, color_list, "filled"); hold off;

xlabel('Azimuthal angle / {\pi}'); ylabel('Polar angle / {\pi}');
axis equal; set(gca, 'YDir', 'reverse');
xlim([-1, 1] + azi_ctr/pi); ylim([0, 1]); yticks([0, 0.5, 1]);
set(gca, 'fontname', ftnm, 'fontsize', ftsz, 'fontweight', 'normal', ...
    'labelfontsizemultiplier', lfsm, 'linewidth', lw, 'Layer', 'Top', 'Box', 'on');
set(gcf, 'unit', 'normalized', 'Position', fpos);
end