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

%
newfigure(4.5,3.75);
set(gca, 'FontSize', 18);
hold on

cmap = colormap(jet(ncolors));
cols = zeros(height(fulltab), 3);
fulltab.H = 0.5*(fulltab.Hx+fulltab.Hy);
fulltab.renormH = fulltab.H - min(fulltab.H);
fulltab.renormH = fulltab.renormH / max(fulltab.renormH);
for ii=1:height(fulltab)
    cind = floor(fulltab.renormH(ii)*(ncolors-1))+1;
    cols(ii,:) = cmap(cind,:);
end


scatter(0.5*(fulltab.tau_x+fulltab.tau_y), fulltab.I, 8, cols, 'filled');
xlabel('$\tau$','Interpreter','Latex');
ylabel('$I$','Interpreter','Latex');
set(gca,'XScale', 'log');
% title('Colored by H');
xlim([1, 100]);
colormap(cmap);
c = colorbar();
c.Ticks = [0, 0.5, 1];
c.TickLabels = {num2str(round(min(fulltab.H)*10)/10), '0', num2str(round(max(fulltab.H)*10)/10)};
yl = ylabel(c,'$H$', 'Interpreter', 'Latex');
yl.Position = [2.2, 0.5, 0];


cols = zeros(height(fulltab), 3);
fulltab.renormTheta = fulltab.Theta-min(fulltab.Theta);
fulltab.renormTheta = fulltab.renormTheta / max(fulltab.renormTheta);
cols = 0*cols;
for ii=1:height(fulltab)
    cind = floor(fulltab.renormTheta(ii)*(ncolors-1))+1;
    cols(ii,:) = cmap(cind,:);
end
print(gcf, '-dpng', 'FreeLunchColH.png', '-r600');

%
newfigure(4.5,3.75);
set(gca, 'FontSize', 18);
hold on
scatter(0.5*(fulltab.tau_x+fulltab.tau_y), fulltab.I, 16, cols, 'filled');
xlabel('$\tau$','Interpreter','Latex');
ylabel('$I$','Interpreter','Latex');
% title('Colored by Theta');
set(gca,'XScale', 'log');
xlim([1,100]);

colormap(cmap);
c = colorbar();
c.Ticks = [0, 0.5, 1];
c.TickLabels = {num2str(round(min(fulltab.Theta)*10)/10), '0', num2str(round(max(fulltab.Theta)*10)/10)};
yl = ylabel(c,'$T$', 'Interpreter', 'Latex');
yl.Position = [2.2, 0.5, 0];
print(gcf, '-dpng', 'FreeLunchColT.png', '-r600');

