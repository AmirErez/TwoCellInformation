% Plots collected tau/theta/I as hx+hy=0
% fulltab = readtable('../Data/randparams0_1-nc_3000-schlogl-collected.csv', 'delimiter',',');
% ttl = 'Schlogl';
% outfile = 'ShannonI_rand_collective.png';
fulltab = readtable('../Data/randparams0_1-nc_3000-hill-collected.csv', 'delimiter',',');
ttl = 'Hill';
outfile = 'ShannonI_rand_collective_Hill.png';

hxs = unique(fulltab.h_x);
hys = unique(fulltab.h_y);
theta_xs = unique(fulltab.theta_x);
theta_ys = unique(fulltab.theta_y);

% Scatter plot according to Andrew's predictions

% Hx = [(hx + hy)*g + hx*thetay]
% Hy = [(hx + hy)*g + hy*thetax]
% Theta = thetax*thetay + g*(thetax + thetay)

g=1;
fulltab = fulltab(abs(fulltab.theta_x)<=0.1 & abs(fulltab.theta_y )<=0.1,:);
fulltab = fulltab(abs(fulltab.h_x)<=0.1 & abs(fulltab.h_x )<=0.1,:);
% fulltab = sortrows(fulltab,'I');

fulltab.Hx = (fulltab.h_x+fulltab.h_y)*g+ fulltab.h_x.*fulltab.theta_y;
fulltab.Hy = (fulltab.h_x+fulltab.h_y)*g+ fulltab.h_y.*fulltab.theta_x;
fulltab.Theta = fulltab.theta_x.*fulltab.theta_y + g*(fulltab.theta_x + fulltab.theta_y);
Hx_str = '$g(h_x+h_y)+ h_x \theta_y$';
Hy_str = '$g(h_x+h_y)+ h_y \theta_x$';
H_str = '$H=g(h_x+h_y) + \frac{1}{2}(h_x \theta_y + h_y\theta_x)$';
Theta_str = '$T=\theta_x\theta_y+ g(\theta_x+\theta_y)$';

fulltab.renormI = (fulltab.I-min(fulltab.I));
fulltab.renormI = fulltab.renormI / max(fulltab.renormI);

newfigure(4.5,3.75);
% Ivals = 0:0.5:2;
if(max(fulltab.I)>1.5)
    Ivals = [0, 0.5, 1, 1.5, max(fulltab.I)];
else
    Ivals = [0, 0.5, 1, 1.5];
end
ncolors = 1024;
Ivalsinds = linspace(min(Ivals), max(Ivals), ncolors);
cmap = colormap(hot(ncolors));
hold on
cols = zeros(height(fulltab), 3);
for ii=1:height(fulltab)
%     cind = floor(fulltab.renormI(ii)*(ncolors-1))+1;
    [mn, cind] = min(abs(fulltab.I(ii)-Ivalsinds));
    cols(ii,:) = cmap(cind,:);
end
% scatter(fulltab.Theta, fulltab.Hx, 50, cols, 'filled');
% xlabel(Theta_str, 'Interpreter','Latex');
% ylabel(Hx_str, 'Interpreter','Latex');
% 
% figure;
% scatter(fulltab.Theta, fulltab.Hy, 50, cols, 'filled');
% xlabel(Theta_str, 'Interpreter','Latex');
% ylabel(Hy_str, 'Interpreter','Latex');

scatter(0.5*(fulltab.Hx+fulltab.Hy), fulltab.Theta,50, cols, 'filled','s');
% ylim([-0.2, 0.2]);
% xlim([-0.2, 0.2]);
ylim([-0.1, 0.1]);
xlim([-0.1, 0.1]);
ylabel(Theta_str, 'Interpreter','Latex');
xlabel(H_str, 'Interpreter','Latex');
c = colorbar;
if (max(Ivals)==1.5)
    c.Ticks = Ivals/max(Ivals);
    c.TickLabels = Ivals;
else
    c.Ticks = Ivals(1:end-1)/max(Ivals);
    c.TickLabels = Ivals(1:end-1);    
end
% c.Ticks = [0, 0.5, 1];
% c.TickLabels = [min(fulltab.I), (min(fulltab.I)+max(fulltab.I))/2, max(fulltab.I)];
set(gca,'FontSize',18);
t = text(0.127, 0.107, '$I$', 'Interpreter', 'Latex', 'FontSize', 18);
t1 = text(-0.177, 0.107, '(c)', 'Interpreter', 'Latex', 'FontSize', 22);
print(gcf,'-dpng',outfile,'-r600');


