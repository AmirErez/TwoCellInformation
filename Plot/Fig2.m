% Used to create Fig.2 in the manuscript

% figname = 'Fig2';
% newfigure(3.375,3.375);
% set(gca,'FontSize',9);
% nmap = colormap(parula(length(collected.ncs)+5));
% [ha, pos] = tight_subplot(2,2,[.05 .08],[.08 .01],[.08 .01]);

% manythetasfile = '../Data/ProductionManyThetas.csv';
many_thetas_file = '../Data/ProductionManyThetasParticularNcs.csv';
theta_zero_file = '../Data/ProductionTheta0Scaling.csv';

tabManyThetas = readtable(many_thetas_file);
tabManyThetas = sortrows(tabManyThetas,'nc');
tabThetaZero = readtable(theta_zero_file);
tabThetaZero = sortrows(tabThetaZero,'nc');


%% a - Cv
figname = 'Fig2';
newfigure(3.375,3.375/2.2);
set(gca,'FontSize',9);

subplot(1,2,1);
ncs = unique(tabManyThetas.nc);
nmap = colormap(parula(length(ncs)+4));
hold on
Cv_at_thetazero = NaN*zeros(length(ncs),1);
min_Cv = NaN*zeros(length(ncs),1);
for nn=1:length(ncs)
    f = find((tabManyThetas.g==1) .* (tabManyThetas.nc==ncs(nn)));
    subtab = tabManyThetas(f,:);
    subtab = sortrows(subtab,'theta');
    
    Cv = (1+subtab.theta(1:end-1)).*diff(subtab.Snm)./diff(subtab.theta);
    
    plot(subtab.theta(1:end-1), Cv,'.-','Color',nmap(nn,:),'DisplayName',[num2str(log10(ncs(nn)),3)]);
%     plot(subtab.theta(1:end-1)*sqrt(ncs(nn)), Cv/sqrt(ncs(nn)),'.-','Color',nmap(nn,:),'DisplayName',[num2str(log10(ncs(nn)),3)]);    
    %       plot(log(abs(subtab.theta(2:end))), log(abs(Cv)),'.-','Color',nmap(nn,:));


    f = find(subtab.theta==0);
    Cv_at_thetazero(nn) = Cv(f);
    min_Cv(nn) = min(Cv(subtab.theta(1:end-1)>-0.1));
end
xlabel('$\theta$','Interpreter','latex');
ylabel('$C_v$','Interpreter','latex');
ylim([-25,20]);
xlim([-0.02,0.06]);
text(-0.06,15,'(a)');

axes('Position',[.33 .77 .11 .2])
box on
plot(sqrt(ncs),-Cv_at_thetazero,'.-b');
% plot(sqrt(ncs.^(0.4)),-Cv_at_thetazero,'.-b');

set(gca,'FontSize',9);
hold on
% plot(sqrt(ncs),0.35*sqrt(ncs)-0.38,'-k','LineWidth',2);
% xlim(10.^[2.7,4]);
% xticks(10.^[3,4]);
ylim([5,30]);
xl = xlabel('$\sqrt{n_c}$','Interpreter','latex');
% xl = xlabel('$n_c^{0.4}$','Interpreter','latex');
xl.Position(2) = -3;
xl.Position(1) = xl.Position(1)-15;
xticks([0,100]);
ylabel('$|C_v|$','Interpreter','latex');
% yticks([5,15,35]);


%% Corr time (b)

% ax = subplot(1,2,2);
% ax.Position = [ax.Position(1), ax.Position(2)+ax.Position(4)/1.8, ax.Position(3), ax.Position(4)/2.1];
ax = axes('Parent', gcf);
ax.Position = [ax.Position(1)+ax.Position(3)/1.8, ax.Position(2)+ax.Position(4)/1.8, ax.Position(3)/2.1, ax.Position(4)/2.1];

set(gca,'FontSize',9);
hold on

% gs = unique(tabThetaZero.g);
gs = [0,1];
% gmap = colormap(copper(length(gs)));
gmap(1,:) = [0 0 0.7];
gmap(2,:) = [0.9 0 0];

slopes = NaN(length(gs),1);
slopes_err = zeros(length(gs),1);

local_slopes = cell(length(gs),1);

for gg=1:length(gs)    
    f = find((tabThetaZero.g==gs(gg)).*(tabThetaZero.theta==0));
    subtab = tabThetaZero(f,:);        
    ncs = subtab.nc;
    local_slopes{gg} = nan(length(ncs),3);
    
    loglog(subtab.nc,subtab.tau_n,'o','MarkerSize',2,'Color',gmap(gg,:),'MarkerFaceColor',gmap(gg,:),'DisplayName',['g=' num2str(gs(gg))]);
    hold on
    pl = loglog(subtab.nc,subtab.tau_m,'x','MarkerSize',2,'Color',gmap(gg,:),'MarkerFaceColor',gmap(gg,:),'DisplayName',['log_{10} g=' num2str(log10(gs(gg)),2)]);
    set(get(get( pl, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off' );

%     for nn=4:(length(ncs)-3)
      for nn=1:length(ncs)
        %     ft = fit(log(subtab.nc(flarge)),log(subtab.tau_n(flarge)/subtab.tau_n(1)),'poly1');
        nstart = nn-3;
        nend = nn+3;
        if(nstart<1), nstart=1; end
        if(nend>length(ncs)), nend=length(ncs); end        
        local_ncs = subtab.nc(nstart:nend);
        local_tau_n = subtab.tau_n(nstart:nend);
        local_tau_m = subtab.tau_m(nstart:nend);
        ft_n = fit(log(local_ncs),log(local_tau_n),'poly1');
        ci_n = confint(ft_n);
        ft_m = fit(log(local_ncs),log(local_tau_m),'poly1');
        ci_m = confint(ft_m);
        
        local_slopes{gg}(nn,1) = subtab.nc(nn);
        local_slopes{gg}(nn,2) = (ft_n.p1+ft_m.p1)/2;
        local_slopes{gg}(nn,3) = sqrt((ci_n(2,1)-ft_n.p1)^2+(ci_m(2,1)-ft_m.p1)^2)/sqrt(2);
    end
end
text(10^2.2,10^2.4,'(b)');
xl = xlabel('$n_c$','Interpreter','Latex')
xl.Position = [10^4.2, 15,-1];
yl = ylabel('$\tau$','Interpreter','Latex');
yl.Position = [90, 300,-1];
xlim(10.^[2,4]);
ylim([8, 10.^2.5]);
set(gca,'YTick',10.^[1,2]);
set(gca,'XTick',10.^[3,4]);
set(gca,'XScale','Log');
set(gca,'YScale','Log');
box('off');

[lb,licons,lplots,ltxt]= legend('show');
lb.Box = 'Off';
lb.Position = [0.61523        0.88       0.2284      0.10082];
licons(1).Position(1) = 0.5;
licons(2).Position(1) = 0.5;

%% (c) Exponent x in tau ~ n_c^x

ax = axes('Parent', gcf);
ax.Position = [ax.Position(1)+ax.Position(3)/1.8, ax.Position(2), ax.Position(3)/2.1, ax.Position(4)/2.1];
set(gca,'FontSize',9);
hold on
ax.Color = 'None';

for gg=1:length(gs)
%     errorbar(local_slopes{gg}(:,1),local_slopes{gg}(:,2),local_slopes{gg}(:,3),'.-','Color',gmap(gg,:))
    plot(local_slopes{gg}(:,1),local_slopes{gg}(:,2),'.-','Color',gmap(gg,:))
    maxslopes = local_slopes{gg}(:,2)+local_slopes{gg}(:,3);
    minslopes = local_slopes{gg}(:,2)-local_slopes{gg}(:,3);    
    h = area(local_slopes{gg}(:,1),[minslopes'; maxslopes'-minslopes']','LineStyle','None');
    h(1).FaceAlpha = 0;
    h(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
    h(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
    h(2).FaceAlpha = 0.3;
    h(2).FaceColor = gmap(gg,:);        
end

set(gca,'XScale','log');
text(10^2.2,1,'(c)');
ylim([0.4,1]);
xlim(10.^[2,4]);
xl = xlabel('$n_c$','Interpreter','Latex')
xl.Position = [10^4.2, 0.5,-1];
yl = ylabel('$x$','Interpreter','Latex');
yl.Position = [90,0.95,-1];
set(gca,'YTick',[0.5,0.75]);

print(gcf,'-dpng','Fig2_Composite.png','-r600');

