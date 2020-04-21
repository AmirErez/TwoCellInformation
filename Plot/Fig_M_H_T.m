fulltab = readtable('../Data/randparams0_1-nc_3000-schlogl-collected.csv', 'delimiter',',');
ttl = 'Schlogl';
warning('Hack! g hardcoded to 1');

fulltab.g = ones(height(fulltab),1);

ncolors = 1024;
figure;
cmap = colormap(parula(ncolors));
close(gcf);

fulltab.Hx = (fulltab.h_x+fulltab.h_y).*fulltab.g+ fulltab.h_x.*fulltab.theta_y;
fulltab.Hy = (fulltab.h_x+fulltab.h_y).*fulltab.g+ fulltab.h_y.*fulltab.theta_x;
fulltab.H = (fulltab.Hx+fulltab.Hy)/2;
fulltab.Theta = fulltab.theta_x.*fulltab.theta_y + fulltab.g.*(fulltab.theta_x + fulltab.theta_y);
Hx_str = '$g(h_x+h_y)+ h_x \theta_y$';
Hy_str = '$g(h_x+h_y)+ h_y \theta_x$';
H_str = '$g(h_x+h_y) + \frac{1}{2}(h_x \theta_y + h_y\theta_x)$';
Theta_str = '$\theta_x\theta_y+ g(\theta_x+\theta_y)$';

% Plot mean_n vs. H and Theta
newfigure(4.5,3.75);
set(gca, 'FontSize', 18);
hold on
cols = zeros(height(fulltab), 3);
fulltab.renormTheta = fulltab.Theta-min(fulltab.Theta);
fulltab.renormTheta = fulltab.renormTheta / max(fulltab.renormTheta);
cols = 0*cols;
for ii=1:height(fulltab)
    cind = floor(fulltab.renormTheta(ii)*(ncolors-1))+1;
    cols(ii,:) = cmap(cind,:);
end
nc = 3000;
scatter(fulltab.H, (0.5*(fulltab.mean_n+fulltab.mean_m)-nc)/nc,10, cols, 'Filled');
xlabel('$H$','Interpreter', 'latex');
ylabel('$M$', 'Interpreter', 'Latex')

% title('Colored by Theta');
%H - T m - m^3/3 == 0
% Hs = linspace(-0.2, 0.2, 201);
% Thetas = [-0.2:0.1:0.2];
% vals = zeros(length(Hs), length(Thetas));
% for tt=1:length(Thetas)
%     for hh=1:length(Hs)
%         root = roots([-1/3, 0, -Thetas(tt), Hs(hh)]);
%         [~, mind] = min(abs(imag(root)));
%         vals(hh, tt) = root(mind);
%     end
%     plot(Hs, vals(:, tt), '-', 'DisplayName', ['$\Theta=' num2str(Thetas(tt)) '$']);
% end
colormap(cmap);
c = colorbar();
xlim([-0.2, 0.2]);
ylim([-1, 1]);
c.Ticks = [0, 0.5, 1];
c.TickLabels = {num2str(round(min(fulltab.Theta)*10)/10), '0', num2str(round(max(fulltab.Theta)*10)/10)};
cy = ylabel(c, '$T$', 'Interpreter', 'Latex');
cy.Position = [2.265          0.5            0];
print(gcf, 'SI-Fig6a.png', '-dpng', '-r600');

newfigure(4.5,3.75); 
set(gca, 'FontSize', 18);
hold on
cols = zeros(height(fulltab), 3);
fulltab.renormTheta = fulltab.Theta-min(fulltab.Theta);
fulltab.renormTheta = fulltab.renormTheta / max(fulltab.renormTheta);
cols = 0*cols;
fulltab.renormH = fulltab.H - min(fulltab.H);
fulltab.renormH = fulltab.renormH / max(fulltab.renormH);
for ii=1:height(fulltab)
    cind = floor(fulltab.renormH(ii)*(ncolors-1))+1;
    cols(ii,:) = cmap(cind,:);
end
% Thetas = linspace(-0.2, 0.2, 201);
% Hs = [-0.2:0.1:0.2];
% vals = zeros(length(Thetas), length(Hs));
% for hh=1:length(Hs)
%     for tt=1:length(Thetas)    
%         root = roots([-1/3, 0, -Thetas(tt), Hs(hh)]);
%         [~, mind] = min(abs(imag(root)));
%         vals(tt, hh) = root(mind);
%     end
%     plot(Thetas, vals(:, hh), '-', 'DisplayName', ['$H=' num2str(Hs(hh)) '$']);
% end
colormap(cmap);
c = colorbar();
c.Ticks = [0, 0.5, 1];
c.TickLabels = {num2str(round(min(fulltab.H)*10)/10), '0', num2str(round(max(fulltab.H)*10)/10)};
cy = ylabel(c, '$H$', 'Interpreter', 'Latex');
cy.Position = [2.265          0.5            0];
% legend();
scatter(fulltab.Theta, (0.5*(fulltab.mean_n+fulltab.mean_m)-nc)/nc, 10, cols, 'Filled');
ylabel('$M$', 'Interpreter', 'Latex')
xlabel('$T$', 'Interpreter', 'Latex');
xlim([-0.2, 0.2]);
ylim([-1, 1]);
print(gcf, 'SI-Fig6b.png', '-dpng', '-r600');
