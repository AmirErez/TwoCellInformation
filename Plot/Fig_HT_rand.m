fulltab = readtable('../Data/randparams0_1-nc_3000-schlogl-long-collected.csv', 'delimiter',',');
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

%%

fulltab.renormI = (fulltab.I-min(fulltab.I));
fulltab.renormI = fulltab.renormI / max(fulltab.renormI);
fulltab = sortrows(fulltab, 'renormI');

newfigure(4.5,3.75);
set(gca, 'FontSize', 18);
hold on

ncolors = 1024;
cmap = colormap(parula(ncolors));
cols = zeros(height(fulltab), 3);
for ii=1:height(fulltab)
    cind = floor(fulltab.renormI(ii)*(ncolors-1))+1;
    cols(ii,:) = cmap(cind,:);
end
scatter(fulltab.Theta, fulltab.Hx, 50, cols, 'filled');
xlabel(Theta_str, 'Interpreter','Latex');
ylabel(Hx_str, 'Interpreter','Latex');

newfigure(4.5,3.75);
set(gca, 'FontSize', 18);
hold on
scatter(fulltab.Theta, fulltab.Hy, 50, cols, 'filled');
xlabel(Theta_str, 'Interpreter','Latex');
ylabel(Hy_str, 'Interpreter','Latex');

newfigure(4.5,3.75);
set(gca, 'FontSize', 18);
hold on
scatter(fulltab.Theta, 0.5*(fulltab.Hx+fulltab.Hy), 50, cols, 'filled');
xlabel(Theta_str, 'Interpreter','Latex');
ylabel(H_str, 'Interpreter','Latex');
