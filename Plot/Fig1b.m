% Generates Fig1b: mutual info examples

datadir = '../Data/ForFig1b';

datanames = {'out_nc_1000_theta_-0.1_g_0.mat','out_nc_1000_theta_0_g_0.mat',...
    'out_nc_1000_theta_0.1_g_0.mat',...
    'out_nc_1000_theta_-0.1.mat','out_nc_1000_theta_0.mat',...
    'out_nc_1000_theta_0.1.mat'};

xlabels = cell(6,1);
xlabels{4} = '$\theta = -0.1$';
xlabels{5} = '$\theta =  0$';
xlabels{6} = '$\theta =  0.1$';


%% Fig 1b
figname = 'JointXY';

newfigure(3.375,2.23);
%colormap(flipud(hot));

% R=[zeros(1,150), tanh((1:106)/45)]';
% G=[zeros(1,90), tanh((1:166)/60)]';
% B=tanh((1:256)/50)';
B=[zeros(1,50), (1:176)/176, ones(1,30)]';
R=[zeros(1,50), (1:206)/206]';
G=[(1:100)/100, ones(1,156)]';
colormap(flipud([R,G,B]));

[ha, pos] = tight_subplot(2,3,[.01 .005],[.11 .01],[.08 .08]);
for ff=1:6
    axes(ha(ff));
    loaded = load([datadir filesep datanames{ff}]);
    P = loaded.GillespieOut.Pnm / max(loaded.GillespieOut.Pnm(:));
%     P(P<1E-7) = 0;
    imagesc(log(P));
%     imagesc(P);
    set(gca,'YDir','Normal');
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
%     if(ff==1||ff==4)
%         ylabel('Y');
%         set(gca,'YTick',[0,1000,2000]);
%     end
    if(ff>3)
        l = xlabel(xlabels{ff},'Interpreter','Latex','FontSize',11);
%         l.Position = [l.Position(1), l.Position(2)+100, l.Position(3)];
%         set(gca,'XTick',[0,1000,2000]);
    end
    
    if(ff==1)
        l = ylabel('$g=0$','Interpreter','Latex');
    elseif(ff==4)
        l = ylabel('$g=1$','Interpreter','Latex');
        ar = annotation('arrow');
        ar.Position = [0.1    0.14    0.14    0];
        ar = annotation('arrow');
        ar.Position = [0.1    0.14    0    0.2];
        text(1400,250,'$X$','Interpreter','Latex');
        text(250,1350,'$Y$','Interpreter','Latex');
    end
      

    ylim([0,2200]);
    xlim([0,2200]);
%     if(ff<4)
%         xlim([0,700]);
%         ylim([0,700]);
%     else
%         ylim([0,2000]);
%         xlim([0,2000]);
%     end
end
c = colorbar();
c.Position = [0.9402    0.1435    0.0370    0.7917];
c.Ticks = [];

% set(gca,'FontSize',11);
% annotation('textbox',...
%     [0.14 0.01 0.3 0.08],...
%     'String',{'$\theta = 0.1$'},'Interpreter','Latex',...
%     'LineStyle','None','FontSize',11);
% text(0.1815,-0.1,'$\theta = 0.1$','Interpreter','Latex');

print(gcf,'-dpng',['Fig1b.png'], '-r600');


