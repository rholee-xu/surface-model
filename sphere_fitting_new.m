%% Sphere fitting NEW
load('Run_difgradients_database.mat','Run_CANew_Con10_MUeq','Run_CANew_Non10_5_MUeq','Run_CANew_Lin10_5_MUeq')
load('Run_difgradients_database.mat','Run_CHNew_Con10_MUeq','Run_CHNew_Non10_5_MUeq','Run_CHNew_Lin10_5_MUeq')
%% SIDE

side_strain_CA_Non10_5 = (max(Run_CANew_Non10_5_MUeq.rv_final128(2,:))-max(Run_CANew_Non10_5_MUeq.rv_init128(2,:)))/(max(Run_CANew_Non10_5_MUeq.rv_init128(2,:)));
Y_CA_Non10_5 = 2*max(Run_CANew_Non10_5_MUeq.rv_final128(2,:))/side_strain_CA_Non10_5;
err_CA_Non10_5 = (young_poisson(10,8)-Y_CA_Non10_5)/young_poisson(10,8)

side_strain_CA_Lin10_5 = (max(Run_CANew_Lin10_5_MUeq.rv_final128(2,:))-max(Run_CANew_Lin10_5_MUeq.rv_init128(2,:)))/(max(Run_CANew_Lin10_5_MUeq.rv_init128(2,:)));
Y_CA_Lin10_5 = 2*max(Run_CANew_Lin10_5_MUeq.rv_final128(2,:))/side_strain_CA_Lin10_5;
err_CA_Lin10_5 = (young_poisson(10,8)-Y_CA_Lin10_5)/young_poisson(10,8)

side_strain_CA_Con10 = (max(Run_CANew_Con10_MUeq.rv_final128(2,:))-max(Run_CANew_Con10_MUeq.rv_init128(2,:)))/(max(Run_CANew_Con10_MUeq.rv_init128(2,:)));
Y_CA_Con10 = 2*max(Run_CANew_Con10_MUeq.rv_final128(2,:))/side_strain_CA_Con10;
err_CA_Con10 = (young_poisson(10,8)-Y_CA_Con10)/young_poisson(10,8)

side_strain_CH_Non10_5 = (max(Run_CHNew_Non10_5_MUeq.rv_final128(2,:))-max(Run_CHNew_Non10_5_MUeq.rv_init128(2,:)))/(max(Run_CHNew_Non10_5_MUeq.rv_init128(2,:)));
Y_CH_Non10_5 = 2*max(Run_CHNew_Non10_5_MUeq.rv_final128(2,:))/side_strain_CH_Non10_5;
err_CH_Non10_5 = (young_poisson(10,8)-Y_CH_Non10_5)/young_poisson(10,8)

side_strain_CH_Lin10_5 = (max(Run_CHNew_Lin10_5_MUeq.rv_final128(2,:))-max(Run_CHNew_Lin10_5_MUeq.rv_init128(2,:)))/(max(Run_CHNew_Lin10_5_MUeq.rv_init128(2,:)));
Y_CH_Lin10_5 = 2*max(Run_CHNew_Lin10_5_MUeq.rv_final128(2,:))/side_strain_CH_Lin10_5;
err_CH_Lin10_5 = (young_poisson(10,8)-Y_CH_Lin10_5)/young_poisson(10,8)

side_strain_CH_Con10 = (max(Run_CHNew_Con10_MUeq.rv_final128(2,:))-max(Run_CHNew_Con10_MUeq.rv_init128(2,:)))/(max(Run_CHNew_Con10_MUeq.rv_init128(2,:)));
Y_CH_Con10 = 2*max(Run_CHNew_Con10_MUeq.rv_final128(2,:))/side_strain_CH_Con10;
err_CH_Con10 = (young_poisson(10,8)-Y_CH_Con10)/young_poisson(10,8)
%%
n = 1;
for i = 10:10:60
    [r,r_after,~,y_tip]=K_radius_tip(Run_CANew_Non10_5_MUeq.rv_final128,Run_CANew_Non10_5_MUeq.rv_init128,Run_CANew_Non10_5_MUeq.synTensC128,Run_CANew_Non10_5_MUeq.synTensM128,i);
    rad_CA_Non10_5(n) = r;
    rad_after_CA_Non10_5(n) = r_after;
    Ytip_CA_Non10_5(n) = y_tip;
    % hold on;scatter(i,y_tip,'filled')
    n = n+1;
end
% yline(young_poisson(Run_CANew_Non10_5.K128(end),Run_CANew_Non10_5.MU128(end)))
n = 1;
for i = 10:10:60
    [r,r_after,~,y_tip]=K_radius_tip(Run_CANew_Lin10_5_MUeq.rv_final128,Run_CANew_Lin10_5_MUeq.rv_init128,Run_CANew_Lin10_5_MUeq.synTensC128,Run_CANew_Lin10_5_MUeq.synTensM128,i);
    rad_CA_Lin10_5(n) = r;
    rad_after_CA_Lin10_5(n) = r_after;
    Ytip_CA_Lin10_5(n) = y_tip;
    % hold on;scatter(i,y_tip,'filled')
    n = n+1;
end
% yline(young_poisson(Run_CANew_Lin10_5.K128(end),Run_CANew_Lin10_5.MU128(end)))
n = 1;
for i = 10:10:60
    [r,r_after,~,y_tip]=K_radius_tip(Run_CANew_Con10_MUeq.rv_final128,Run_CANew_Con10_MUeq.rv_init128,Run_CANew_Con10_MUeq.synTensC128,Run_CANew_Con10_MUeq.synTensM128,i);
    rad_CA_Con10(n) = r;
    rad_after_CA_Con10(n) = r_after;
    Ytip_CA_Con10(n) = y_tip;
    % hold on;scatter(i,y_tip,'filled')
    n = n+1;
end
% yline(young_poisson(Run_CANew_Con10.K128(end),Run_CANew_Con10.MU128(end)))
%%
n = 1;
for i = 10:10:60
    [r,r_after,~,y_tip]=K_radius_tip(Run_CHNew_Non10_5_MUeq.rv_final128,Run_CHNew_Non10_5_MUeq.rv_init128,Run_CHNew_Non10_5_MUeq.synTensC128,Run_CHNew_Non10_5_MUeq.synTensM128,i);
    rad_CH_Non10_5(n) = r;
    rad_after_CH_Non10_5(n) = r_after;
    Ytip_CH_Non10_5(n) = y_tip;
    % hold on;scatter(i,y_tip,'filled')
    n = n+1;
end
% yline(young_poisson(Run_CHNew_Non10_5.K128(end),Run_CHNew_Non10_5.MU128(end)))

n = 1;
for i = 10:10:60
    [r,r_after,~,y_tip]=K_radius_tip(Run_CHNew_Lin10_5_MUeq.rv_final128,Run_CHNew_Lin10_5_MUeq.rv_init128,Run_CHNew_Lin10_5_MUeq.synTensC128,Run_CHNew_Lin10_5_MUeq.synTensM128,i);
    rad_CH_Lin10_5(n) = r;
    rad_after_CH_Lin10_5(n) = r_after;
    Ytip_CH_Lin10_5(n) = y_tip;
    % hold on;scatter(i,y_tip,'filled')
    n = n+1;
end
% yline(young_poisson(Run_CHNew_Lin10_5.K128(end),Run_CHNew_Lin10_5.MU128(end)))

n = 1;
for i = 10:10:60
    [r,r_after,~,y_tip]=K_radius_tip(Run_CHNew_Con10_MUeq.rv_final128,Run_CHNew_Con10_MUeq.rv_init128,Run_CHNew_Con10_MUeq.synTensC128,Run_CHNew_Con10_MUeq.synTensM128,i);
    rad_CH_Con10(n) = r;
    rad_after_CH_Con10(n) = r_after;
    Ytip_CH_Con10(n) = y_tip;
    % hold on;scatter(i,y_tip,'filled')
    n = n+1;
end
% yline(young_poisson(Run_CHNew_Con10.K128(end),Run_CHNew_Con10.MU128(end)))
%%
Run = Run_CANew_Lin10_5_MUeq;
n = 1;
for i = 10:10:60
    for j = 1:10

        scale = (1+0.4*rand(1,1));
        rv_final128 = Run.rv_final128.*scale;
        rv_init128 = Run.rv_init128.*scale;
        [~,~,~,~,~,y_tip]=K_radius_tip(rv_final128,rv_init128,Run.synTensC128,Run.synTensM128,i);
        avgY_LinCANew(j,n) = y_tip;
    end
    n = n+1;
end

Run = Run_CANew_Non10_5_MUeq;
n = 1;
for i = 10:10:60
    for j = 1:10

        scale = (1+0.4*rand(1,1));
        rv_final128 = Run.rv_final128.*scale;
        rv_init128 = Run.rv_init128.*scale;
        [~,~,~,~,~,y_tip]=K_radius_tip(rv_final128,rv_init128,Run.synTensC128,Run.synTensM128,i);
        avgY_NonCANew(j,n) = y_tip;
    end
    n = n+1;
end

Run = Run_CANew_Con10_MUeq;
n = 1;
for i = 10:10:60
    for j = 1:10

        scale = (1+0.4*rand(1,1));
        rv_final128 = Run.rv_final128.*scale;
        rv_init128 = Run.rv_init128.*scale;
        [~,~,~,~,~,y_tip]=K_radius_tip(rv_final128,rv_init128,Run.synTensC128,Run.synTensM128,i);
        avgY_ConCANew(j,n) = y_tip;
    end
    n = n+1;
end
%%
Run = Run_CHNew_Lin10_5_MUeq;
n = 1;
for i = 10:10:60
    for j = 1:10

        scale = (1+0.4*rand(1,1));
        rv_final128 = Run.rv_final128.*scale;
        rv_init128 = Run.rv_init128.*scale;
        [~,~,~,~,~,y_tip]=K_radius_tip(rv_final128,rv_init128,Run.synTensC128,Run.synTensM128,i);
        avgY_LinCHNew(j,n) = y_tip;
    end
    n = n+1;
end

Run = Run_CHNew_Non10_5_MUeq;
n = 1;
for i = 10:10:60
    for j = 1:10

        scale = (1+0.4*rand(1,1));
        rv_final128 = Run.rv_final128.*scale;
        rv_init128 = Run.rv_init128.*scale;
        [~,~,~,~,~,y_tip]=K_radius_tip(rv_final128,rv_init128,Run.synTensC128,Run.synTensM128,i);
        avgY_NonCHNew(j,n) = y_tip;
    end
    n = n+1;
end

Run = Run_CHNew_Con10_MUeq;
n = 1;
for i = 10:10:60
    for j = 1:10

        scale = (1+0.4*rand(1,1));
        rv_final128 = Run.rv_final128.*scale;
        rv_init128 = Run.rv_init128.*scale;
        [~,~,~,~,~,y_tip]=K_radius_tip(rv_final128,rv_init128,Run.synTensC128,Run.synTensM128,i);
        avgY_ConCHNew(j,n) = y_tip;
    end
    n = n+1;
end

%% sphere fitting radius results
tiledlayout(2,3);
nexttile;hold on;scatter(10:10:60,rad_CA_Con10(1,1:6),'filled');%title('Constant','FontSize',13);
hold on;scatter(10:10:60,rad_after_CA_Con10(1,1:6),'filled');
% e=errorbar(10:10:60,mean(rad_Con_CANew(:,1:6)),std(rad_Con_CANew(:,1:6)),'LineStyle','none');e.LineWidth = 2;
% e=errorbar(10:10:60,mean(radafter_Con_CANew(:,1:6)),std(radafter_Con_CANew(:,1:6)),'LineStyle','none');e.LineWidth = 2;
xlabel('Degree fitting','FontSize',12);ylabel('Sphere Radius','FontSize',12);set(gca,'LineWidth',1,'FontSize',13,'FontWeight','bold');
xlim([10,65]);ylim([0.45,0.75]);xticks([0,10,20,30,40,50,60]);yticks([0.45:0.1:0.75]);legend('Radius turgid (R_B)','Radius unturgid (R_A)','Fontweight','normal','Location','northoutside')

nexttile;hold on;scatter(10:10:60,rad_CA_Non10_5(1,1:6),'filled');title('Nonlinear','FontSize',13);
hold on;scatter(10:10:60,rad_after_CA_Non10_5(1,1:6),'filled');
xlabel('Degree fitting','FontSize',12);ylabel('Sphere Radius','FontSize',12);set(gca,'LineWidth',1,'FontSize',13,'FontWeight','bold');
xlim([10,65]);ylim([0.45,0.75]);xticks([0,10,20,30,40,50,60]);yticks([0.45:0.1:0.75]);%legend('Radius turgid','Radius unturgid','Fontweight','normal')

nexttile;hold on;scatter(10:10:60,rad_CA_Lin10_5(1,1:6),'filled');title('Linear','FontSize',13);
hold on;scatter(10:10:60,rad_after_CA_Lin10_5(1,1:6),'filled');
xlabel('Degree fitting','FontSize',12);ylabel('Sphere Radius','FontSize',12);set(gca,'LineWidth',1,'FontSize',13,'FontWeight','bold');
xlim([10,65]);ylim([0.45,0.75]);xticks([0,10,20,30,40,50,60]);yticks([0.45:0.1:0.75]);%legend('Radius turgid','Radius unturgid','Fontweight','normal')


nexttile;
hold on;scatter(10:10:60,rad_CH_Con10(1,1:6)','filled');title('Constant','FontSize',13);
hold on;scatter(10:10:60,rad_after_CH_Con10(1,1:6)','filled');
xlabel('Degree fitting','FontSize',12);ylabel('Sphere Radius','FontSize',12);set(gca,'LineWidth',1,'FontSize',13,'FontWeight','bold');
legend('Radius turgid (R_B)','Radius unturgid (R_A)','Fontweight','normal','Location','northeast')
xlim([10,65]);ylim([0.6,0.9]);xticks([0,10,20,30,40,50,60]);yticks([0.6,0.7,0.8,0.9,1]);%legend('Radius turgid','Radius unturgid','Fontweight','normal')


nexttile;hold on;scatter(10:10:60,rad_CH_Non10_5(1,1:6),'filled');title('Nonlinear','FontSize',13);
hold on;scatter(10:10:60,rad_after_CH_Non10_5(1,1:6),'filled');
xlabel('Degree fitting','FontSize',12);ylabel('Sphere Radius','FontSize',12);set(gca,'LineWidth',1,'FontSize',13,'FontWeight','bold')
xlim([10,65]);ylim([0.6,0.9]);xticks([0,10,20,30,40,50,60]);yticks([0.6,0.7,0.8,0.9,1]);%legend('Radius turgid','Radius unturgid','Fontweight','normal')

nexttile;hold on;scatter(10:10:60,rad_CH_Lin10_5(1,1:6),'filled');title('Linear','FontSize',13);
hold on;scatter(10:10:60,rad_after_CH_Lin10_5(1,1:6),'filled');
xlabel('Degree fitting','FontSize',12);ylabel('Sphere Radius','FontSize',12);set(gca,'LineWidth',1,'FontSize',13,'FontWeight','bold')
xlim([10,65]);ylim([0.6,0.9]);xticks([0,10,20,30,40,50,60]);yticks([0.6,0.7,0.8,0.9,1]);%legend('Radius turgid','Radius unturgid','Fontweight','normal')

ax = axes('Position',[.24 .62 .1 .1]);
K_radius_tip(Run_CANew_Con10_MUeq.rv_final128,Run_CANew_Con10_MUeq.rv_init128,Run_CANew_Con10_MUeq.synTensC128,Run_CANew_Con10_MUeq.synTensM128,60,1);
grid on;box off;xticks([]);yticks([]);xlabel([]);ylabel([]);ax.XColor = 'w';ax.YColor= 'w';title('60 degrees')

ax2 = axes('Position',[.1 .62 .1 .1]);
K_radius_tip(Run_CANew_Con10_MUeq.rv_final128,Run_CANew_Con10_MUeq.rv_init128,Run_CANew_Con10_MUeq.synTensC128,Run_CANew_Con10_MUeq.synTensM128,20,1);
grid on;box off;xticks([]);yticks([]);xlabel([]);ylabel([]);ax2.XColor = 'w';ax2.YColor= 'w';title('20 degrees')

ax3 = axes('Position',[.53 .62 .1 .1]);
K_radius_tip(Run_CANew_Non10_5_MUeq.rv_final128,Run_CANew_Non10_5_MUeq.rv_init128,Run_CANew_Non10_5_MUeq.synTensC128,Run_CANew_Non10_5_MUeq.synTensM128,60,1);
grid on;box off;xticks([]);yticks([]);xlabel([]);ylabel([]);ax3.XColor = 'w';ax3.YColor= 'w';title('60 degrees')

ax4 = axes('Position',[.39 .62 .1 .1]);
K_radius_tip(Run_CANew_Non10_5_MUeq.rv_final128,Run_CANew_Non10_5_MUeq.rv_init128,Run_CANew_Non10_5_MUeq.synTensC128,Run_CANew_Lin10_5_MUeq.synTensM128,20,1);
grid on;box off;xticks([]);yticks([]);xlabel([]);ylabel([]);ax4.XColor = 'w';ax4.YColor= 'w';title('20 degrees')

ax5 = axes('Position',[.82 .62 .1 .1]);
K_radius_tip(Run_CANew_Lin10_5_MUeq.rv_final128,Run_CANew_Lin10_5_MUeq.rv_init128,Run_CANew_Lin10_5_MUeq.synTensC128,Run_CANew_Non10_5_MUeq.synTensM128,60,1);
grid on;box off;xticks([]);yticks([]);xlabel([]);ylabel([]);ax5.XColor = 'w';ax5.YColor= 'w';title('60 degrees')

ax6 = axes('Position',[.68 .62 .1 .1]);
K_radius_tip(Run_CANew_Lin10_5_MUeq.rv_final128,Run_CANew_Lin10_5_MUeq.rv_init128,Run_CANew_Lin10_5_MUeq.synTensC128,Run_CANew_Lin10_5_MUeq.synTensM128,20,1);
grid on;box off;xticks([]);yticks([]);xlabel([]);ylabel([]);ax6.XColor = 'w';ax6.YColor= 'w';title('20 degrees')


ax = axes('Position',[.24 .145 .1 .1]);
K_radius_tip(Run_CHNew_Con10_MUeq.rv_final128,Run_CHNew_Con10_MUeq.rv_init128,Run_CHNew_Con10_MUeq.synTensC128,Run_CHNew_Con10_MUeq.synTensM128,60,1);
grid on;box off;xticks([]);yticks([]);xlabel([]);ylabel([]);ax.XColor = 'w';ax.YColor= 'w';title('60 degrees')

ax2 = axes('Position',[.1 .145 .1 .1]);
K_radius_tip(Run_CHNew_Con10_MUeq.rv_final128,Run_CHNew_Con10_MUeq.rv_init128,Run_CHNew_Con10_MUeq.synTensC128,Run_CHNew_Con10_MUeq.synTensM128,20,1);
grid on;box off;xticks([]);yticks([]);xlabel([]);ylabel([]);ax2.XColor = 'w';ax2.YColor= 'w';title('20 degrees')

ax3 = axes('Position',[.53 .145 .1 .1]);
K_radius_tip(Run_CHNew_Non10_5_MUeq.rv_final128,Run_CHNew_Non10_5_MUeq.rv_init128,Run_CHNew_Non10_5_MUeq.synTensC128,Run_CHNew_Non10_5_MUeq.synTensM128,60,1);
grid on;box off;xticks([]);yticks([]);xlabel([]);ylabel([]);ax3.XColor = 'w';ax3.YColor= 'w';title('60 degrees')

ax4 = axes('Position',[.39 .145 .1 .1]);
K_radius_tip(Run_CHNew_Non10_5_MUeq.rv_final128,Run_CHNew_Non10_5_MUeq.rv_init128,Run_CHNew_Non10_5_MUeq.synTensC128,Run_CHNew_Lin10_5_MUeq.synTensM128,20,1);
grid on;box off;xticks([]);yticks([]);xlabel([]);ylabel([]);ax4.XColor = 'w';ax4.YColor= 'w';title('20 degrees')

ax5 = axes('Position',[.82 .145 .1 .1]);
K_radius_tip(Run_CHNew_Lin10_5_MUeq.rv_final128,Run_CHNew_Lin10_5_MUeq.rv_init128,Run_CHNew_Lin10_5_MUeq.synTensC128,Run_CHNew_Non10_5_MUeq.synTensM128,60,1);
grid on;box off;xticks([]);yticks([]);xlabel([]);ylabel([]);ax5.XColor = 'w';ax5.YColor= 'w';title('60 degrees')

ax6 = axes('Position',[.68 .145 .1 .1]);
K_radius_tip(Run_CHNew_Lin10_5_MUeq.rv_final128,Run_CHNew_Lin10_5_MUeq.rv_init128,Run_CHNew_Lin10_5_MUeq.synTensC128,Run_CHNew_Lin10_5_MUeq.synTensM128,20,1);
grid on;box off;xticks([]);yticks([]);xlabel([]);ylabel([]);ax6.XColor = 'w';ax6.YColor= 'w';title('20 degrees')

set(gcf,'color','w')
set(gcf, 'Position', [1 1 1000 500]);
%%
% young_poisson(Run_CANew_Con10_MUeq.K128(end),Run_CANew_Con10_MUeq.MU128(end))
%% sphere fitting modulus results
figure()
tiledlayout(2,3);


nexttile;hold on; yline(young_poisson(Run_CANew_Con10_MUeq.K128(end),Run_CANew_Con10_MUeq.MU128(end))/2,'LineWidth',1);scatter(10:10:60,mean(avgY_ConCANew)/2,30,'filled','MarkerFaceColor',[0.9290 0.6940 0.1250]);title('Constant','FontSize',13);
e=errorbar(10:10:60,mean(avgY_ConCANew)/2,std(avgY_ConCANew)/2,'LineStyle','none');e.LineWidth = 2;e.Color = [0.9290 0.6940 0.1250];e.CapSize = 15;
xlabel('Degree fitting','FontSize',12);ylabel('Surface Modulus','FontSize',12);set(gca,'LineWidth',1,'FontSize',13,'FontWeight','bold');xlim([10,65]);ylim([4,12])
xticks([0,10,20,30,40,50,60]);
nexttile;hold on; yline(young_poisson(Run_CANew_Non10_5_MUeq.K128(end),Run_CANew_Non10_5_MUeq.MU128(end))/2,'LineWidth',1);scatter(10:10:60,mean(avgY_NonCANew)/2,30,'filled','MarkerFaceColor',[0.9290 0.6940 0.1250]);title('Nonlinear','FontSize',13);
e=errorbar(10:10:60,mean(avgY_NonCANew)/2,std(avgY_NonCANew)/2,'LineStyle','none');e.LineWidth = 2;e.Color = [0.9290 0.6940 0.1250];e.CapSize = 15;
xlabel('Degree fitting','FontSize',12);ylabel('Surface Modulus','FontSize',12);set(gca,'LineWidth',1,'FontSize',13,'FontWeight','bold');xlim([10,65]);ylim([2,7])
xticks([0,10,20,30,40,50,60]);
nexttile;hold on; yline(young_poisson(Run_CANew_Lin10_5_MUeq.K128(end),Run_CANew_Lin10_5_MUeq.MU128(end))/2,'LineWidth',1);scatter(10:10:60,mean(avgY_LinCANew)/2,30,'filled','MarkerFaceColor',[0.9290 0.6940 0.1250]);title('Linear','FontSize',13);
e=errorbar(10:10:60,mean(avgY_LinCANew)/2,std(avgY_LinCANew)/2,'LineStyle','none');e.LineWidth = 2;e.Color = [0.9290 0.6940 0.1250];e.CapSize = 15;
xlabel('Degree fitting','FontSize',12);ylabel('Surface Modulus','FontSize',12);set(gca,'LineWidth',1,'FontSize',13,'FontWeight','bold');xlim([10,65]);ylim([2,7])
xticks([0,10,20,30,40,50,60]);

nexttile;
hold on; yline(young_poisson(Run_CHNew_Con10_MUeq.K128(end),Run_CHNew_Con10_MUeq.MU128(end))/2,'LineWidth',1);scatter(10:10:60,mean(avgY_ConCHNew)/2,30,'filled','MarkerFaceColor',[0.9290 0.6940 0.1250]);title('Constant','FontSize',13);
e=errorbar(10:10:60,mean(avgY_ConCHNew)/2,std(avgY_ConCHNew)/2,'LineStyle','none');e.LineWidth = 2;e.Color = [0.9290 0.6940 0.1250];e.CapSize = 15;
xlabel('Degree fitting','FontSize',12);ylabel('Surface Modulus','FontSize',12);set(gca,'LineWidth',1,'FontSize',13,'FontWeight','bold');
xlim([10,65]);xticks([0,10,20,30,40,50,60]);


nexttile;hold on; yline(young_poisson(Run_CHNew_Non10_5_MUeq.K128(end),Run_CHNew_Non10_5_MUeq.MU128(end))/2,'LineWidth',1);scatter(10:10:60,mean(avgY_NonCHNew)/2,30,'filled','MarkerFaceColor',[0.9290 0.6940 0.1250]);title('Nonlinear','FontSize',13);
e=errorbar(10:10:60,mean(avgY_NonCHNew)/2,std(avgY_NonCHNew)/2,'LineStyle','none');e.LineWidth = 2;e.Color = [0.9290 0.6940 0.1250];e.CapSize = 15;
xlabel('Degree fitting','FontSize',12);ylabel('Surface Modulus','FontSize',12);set(gca,'LineWidth',1,'FontSize',13,'FontWeight','bold');
xlim([10,65]);xticks([0,10,20,30,40,50,60]);

nexttile;hold on; yline(young_poisson(Run_CHNew_Lin10_5_MUeq.K128(end),Run_CHNew_Lin10_5_MUeq.MU128(end))/2,'LineWidth',1);scatter(10:10:60,mean(avgY_LinCHNew)/2,30,'filled','MarkerFaceColor',[0.9290 0.6940 0.1250]);title('Linear','FontSize',13);
e=errorbar(10:10:60,mean(avgY_LinCHNew)/2,std(avgY_LinCHNew)/2,'LineStyle','none');e.LineWidth = 2;e.Color = [0.9290 0.6940 0.1250];e.CapSize = 15;
xlabel('Degree fitting','FontSize',12);ylabel('Surface Modulus','FontSize',12);set(gca,'LineWidth',1,'FontSize',13,'FontWeight','bold');
xlim([10,65]);xticks([0,10,20,30,40,50,60]);

axes('Position',[.24 .145 .07 .07]);
% box on
hold on;yline(young_poisson(Run_CHNew_Con10_MUeq.K128(end),Run_CHNew_Con10_MUeq.MU128(end))/2,'LineWidth',1)
scatter(50:10:60,mean(avgY_ConCHNew(:,[5,6]))/2,30,'filled','MarkerFaceColor',[0.9290 0.6940 0.1250]);
xlim([45,62]);ylim([5,inf]);
e=errorbar(50:10:60,mean(avgY_ConCHNew(:,5:6))/2,std(avgY_ConCHNew(:,5:6))/2,'LineStyle','none');e.LineWidth = 2;e.Color = [0.9290 0.6940 0.1250];e.CapSize = 15;

ax3 = axes('Position',[.532 .145 .07 .07]);
% box on
hold on;yline(young_poisson(Run_CHNew_Non10_5_MUeq.K128(end),Run_CHNew_Non10_5_MUeq.MU128(end))/2,'LineWidth',1)
scatter(50:10:60,mean(avgY_NonCHNew(:,[5,6]))/2,30,'filled','MarkerFaceColor',[0.9290 0.6940 0.1250]);
xlim([45,62]);ylim([5,inf]);
e=errorbar(50:10:60,mean(avgY_NonCHNew(:,5:6))/2,std(avgY_NonCHNew(:,5:6))/2,'LineStyle','none');e.LineWidth = 2;e.Color = [0.9290 0.6940 0.1250];e.CapSize = 15;
hold off

ax4 = axes('Position',[.825 .35 .07 .07]);
% box on
hold on;yline(young_poisson(Run_CHNew_Lin10_5_MUeq.K128(end),Run_CHNew_Lin10_5_MUeq.MU128(end))/2,'LineWidth',1)
scatter(50:10:60,mean(avgY_LinCHNew(:,[5,6]))/2,30,'filled','MarkerFaceColor',[0.9290 0.6940 0.1250]);
xlim([45,62]);ylim([5,inf]);
e=errorbar(50:10:60,mean(avgY_LinCHNew(:,5:6))/2,std(avgY_LinCHNew(:,5:6))/2,'LineStyle','none');e.LineWidth = 2;e.Color = [0.9290 0.6940 0.1250];e.CapSize = 15;
hold off

% axes('Position',[.8 .6 .08 .08])
% % box on
% hold on;yline(10)
% scatter(50:10:60,avgY_NonCHNew(:,[5,6]),'filled');xlim([45,62]);ylim([5,inf]);hold off

set(gcf,'color','w')
set(gcf, 'Position', [1 1 1000 500]);
%%
tiledlayout(2,3,'TileSpacing','tight')
nexttile;hold on;
plot(Run_CANew_Con10_MUeq.rv_final128(1,:),Run_CANew_Con10_MUeq.rv_final128(2,:),'LineWidth',2)
plot(Run_CANew_Con10_MUeq.rv_init128(1,:),Run_CANew_Con10_MUeq.rv_init128(2,:),'LineWidth',2);axis equal;xlim([0,4.3]);ylim([0,1.2])

nexttile;hold on;
plot(Run_CANew_Non10_5_MUeq.rv_final128(1,:),Run_CANew_Non10_5_MUeq.rv_final128(2,:),'LineWidth',2)
plot(Run_CANew_Non10_5_MUeq.rv_init128(1,:),Run_CANew_Non10_5_MUeq.rv_init128(2,:),'LineWidth',2);axis equal;xlim([0,4.3]);ylim([0,1.2])

nexttile;hold on;
plot(Run_CANew_Lin10_5_MUeq.rv_final128(1,:),Run_CANew_Lin10_5_MUeq.rv_final128(2,:),'LineWidth',2)
plot(Run_CANew_Lin10_5_MUeq.rv_init128(1,:),Run_CANew_Lin10_5_MUeq.rv_init128(2,:),'LineWidth',2);axis equal;xlim([0,4.3]);ylim([0,1.2])

nexttile;hold on;
plot(Run_CHNew_Con10_MUeq.rv_final128(1,:),Run_CHNew_Con10_MUeq.rv_final128(2,:),'LineWidth',2)
plot(Run_CHNew_Con10_MUeq.rv_init128(1,:),Run_CHNew_Con10_MUeq.rv_init128(2,:),'LineWidth',2);axis equal;xlim([0,4.3]);ylim([0,1.2])

nexttile;hold on;
plot(Run_CHNew_Non10_5_MUeq.rv_final128(1,:),Run_CHNew_Non10_5_MUeq.rv_final128(2,:),'LineWidth',2)
plot(Run_CHNew_Non10_5_MUeq.rv_init128(1,:),Run_CHNew_Non10_5_MUeq.rv_init128(2,:),'LineWidth',2);axis equal;xlim([0,4.3]);ylim([0,1.2])

nexttile;hold on;
plot(Run_CHNew_Lin10_5_MUeq.rv_final128(1,:),Run_CHNew_Lin10_5_MUeq.rv_final128(2,:),'LineWidth',2)
plot(Run_CHNew_Lin10_5_MUeq.rv_init128(1,:),Run_CHNew_Lin10_5_MUeq.rv_init128(2,:),'LineWidth',2);axis equal;xlim([0,4.3]);ylim([0,1.2])

set(gcf,'Color','w','Position',[1,1,1000,500])
%%
function [young,poisson] = young_poisson(K,MU,circ_stretch,mer_stretch)

% poisson = (2*K-2*MU)./(2*K+2*MU);
% young = 2*K.*(1-poisson);

if nargin>2
poisson = (mer_stretch-1)./(circ_stretch-1);
young = 2*K.*(1-poisson);
else
poisson = (2*K-2*MU)./(2*K+2*MU);
young = 2*K.*(1-poisson);
end
end