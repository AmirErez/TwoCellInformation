% plot I(\theta,h) for symmetric cells
% k_{1}^{1} = 1 for both cells
% g = 1 for both cells
% theta & h are varied for both cells, but are always equal,
% i.e. theta1 = theta2, h1 = h2 for every data point
% The color is the Shannon mutual information, in units of nats
% The bottom left corner is cutoff because of the restriction
% h + theta > 1/3

fulltab = readtable('../Data/scan_hx_hy_theta0_nc3000-collected.csv','delimiter',',');
fulltab = fulltab(fulltab.nc_x==3000, :);
thetazerotab = fulltab(fulltab.theta_x == 0 & fulltab.theta_y == 0,:);

h1_vary = unique(thetazerotab.h_x);
h2_vary = unique(thetazerotab.h_y);

info_vary = zeros(length(h2_vary), length(h1_vary));
for h2=1:length(h2_vary)
    for h1=1:length(h1_vary)
        g = thetazerotab(thetazerotab.h_x == h1_vary(h1) & thetazerotab.h_y == h2_vary(h2), :);
        if height(g) ~= 1
            disp(['Warning! g has ' num2str(height(g)) ' rows !']);
            continue
        end
        info_vary(h2, h1) = g.I;
    end
end

newfigure(4.5,3.75);
imagesc(h1_vary,h2_vary,info_vary);
colormap(hot(256));
colorbar;
xlabel('$h_X$','Interpreter','Latex');
ylabel('$h_Y$','Interpreter','Latex');
xlim([-0.3,0.3]);
ylim([-0.3,0.3]);
yticks([-0.3, 0, 0.3]);
xticks([-0.3, 0, 0.3]);
hold on
t = text(0.38, 0.33, '$I$', 'Interpreter', 'Latex', 'FontSize', 18);
t1 = text(-0.52, 0.32, '(b)', 'Interpreter', 'Latex', 'FontSize', 22);
set(gca,'FontSize',18);
set(gca,'YDir','Normal');
print(gcf,'-dpng','ShannonI_nonsymmetric_h.png','-r600');

%% Hill dynamics with constant gamma not g
fulltab = readtable('../Data/scan_hx_hy_theta0_nc3000-hill-collected.csv','delimiter',',');
fulltab = fulltab(fulltab.nc_x==3000, :);
thetazerotab = fulltab(fulltab.theta_x == 0 & fulltab.theta_y == 0,:);

h1_vary = unique(thetazerotab.h_x);
h2_vary = unique(thetazerotab.h_y);

info_vary = zeros(length(h2_vary), length(h1_vary));
for h2=1:length(h2_vary)
    for h1=1:length(h1_vary)
        g = thetazerotab(thetazerotab.h_x == h1_vary(h1) & thetazerotab.h_y == h2_vary(h2), :);
        if height(g) ~= 1
            disp(['Warning! g has ' num2str(height(g)) ' rows !']);
            continue
        end
        info_vary(h2, h1) = g.I;
        if h2_vary(h2)<-0.2 && h1_vary(h1)<-0.2
            info_vary(h2, h1) = 0;
        end
    end
end

newfigure(4.5,3.75);
imagesc(h1_vary,h2_vary,info_vary);
colormap(hot(256));
colorbar;
xlabel('$h_X$','Interpreter','Latex');
ylabel('$h_Y$','Interpreter','Latex');
xlim([-0.3,0.3]);
ylim([-0.3,0.3]);
yticks([-0.3, 0, 0.3]);
xticks([-0.3, 0, 0.3]);
set(gca,'FontSize',18);
set(gca,'YDir','Normal');
t = text(0.45, 0.3, '$I$', 'Interpreter', 'Latex', 'FontSize', 18);
t = text(-0.52, 0.32, '(b)', 'Interpreter', 'Latex', 'FontSize', 22);

print(gcf,'-dpng','ShannonI_nonsymmetric_h_Hill.png','-r600');
