
function Origin = make_plot_Figure_5to9_vari(SEIS,DEIS,change)

Z_SEIS = option_input(SEIS,change);
Z_DEIS = option_input(DEIS,change);

colors_1 = {'g','r','b','m','k',...
            'g','r','b','m','k'};

colors_2 = {'--g','--r','--b','--m','--k',...
            '--g','--r','--b','--m','--k'};

%% OCV
switch change
    case 'cs_bar_neg'
figure(3)
    case 'cs_bar_pos'
figure(4)
    case 'ce'
figure(5)
    case 'phis_neg'
figure(6)
    case 'phis_pos'
figure(7)
    case 'phie'
figure(8)
    case 'j_neg'
figure(9)
    case 'j_pos'
figure(10)
end
subplot(2,3,1)
% draw a 3D line graph using plot3
for m = 1:length(SEIS{1,1}.soc_thr) % soc
    plot3(SEIS{1,1}.soc_thr(m)*ones(length(Z_SEIS{1,1}{m}),1),...
          Z_SEIS{1,1}{m}(:,1),...
          Z_SEIS{1,1}{m}(:,2),colors_2{m},'LineWidth',2);
    hold on
end
xlabel('SOC',              'fontname','Times'); % X axis label
ylabel('{\itx} / {\itL}_4','fontname','Times'); % Y axis label
xlim([0 1])
switch change
    case 'cs_bar_neg'
Fig = 'Fig .5';
zlabel('{\itc}_{ss1} / {\itc}_{s1,max}','fontname','Times'); % Z axis label
ylim([0 0.456])
zlim([0 1])
    case 'cs_bar_pos'
Fig = 'Fig .5';
zlabel('{\itc}_{ss2} / {\itc}_{s2,max}','fontname','Times'); % Z axis label
ylim([0.5855 1])
zlim([0 1])
    case 'ce'
Fig = 'Fig .6';
zlabel('{\itc}_e / {\itc}_{e,0}','fontname','Times'); % Z axis label
ylim([0 1])
zlim([0.5 1.5])
    case 'phis_neg'
Fig = 'Fig .7';
zlabel('{\itφ}_{s1} (V)','fontname','Times'); % Z axis label
ylim([0 0.456])
zlim([0 0.4])
    case 'phis_pos'
Fig = 'Fig .7';
zlabel('{\itφ}_{s2} (V)','fontname','Times'); % Z axis label
ylim([0.5855 1])
zlim([3.6 4.4])
    case 'phie'
Fig = 'Fig .8';
zlabel('{\itφ}_e (V)','fontname','Times'); % Z axis label
ylim([0 1])
zlim([-0.1 0])
    case 'j_neg'
Fig = 'Fig .9';
zlabel('{\itj}_1 × {\itF} / {\iti}','fontname','Times'); % Z axis label
ylim([0 0.456])
zlim([0 0.025])
    case 'j_pos'
Fig = 'Fig .9';
zlabel('{\itj}_2 × {\itF} / {\iti}','fontname','Times'); % Z axis label
ylim([0.5855 1])
zlim([-0.025 0])
set(gca,'FontSize',18,'zDir','reverse')
end
set(gca,'FontSize',18,'yDir','reverse')
grid on; % show grid
switch change
    case {'cs_bar_neg','ce','phis_neg','phie','j_neg'}
        title([Fig,'(a) OCV'],'fontsize',18,'fontweight','bold')
    case {'cs_bar_pos','phis_pos','j_pos'}
        title([Fig,'(e) OCV'],'fontsize',18,'fontweight','bold')
end

%% different C-rates
for k = 1:size(DEIS,1) % C
switch change
    case {'cs_bar_neg','ce','phis_neg','phie','j_neg'}
        if k == 1
            Label = '(b) ';  n_k = 2;
        elseif k == 2
            Label = '(c) ';  n_k = 4;
        elseif k == 3
            Label = '(d) ';  n_k = 5;
        end
    case {'cs_bar_pos','phis_pos','j_pos'}
        if k == 1
            Label = '(f) ';  n_k = 2;
        elseif k == 2
            Label = '(g) ';  n_k = 4;
        elseif k == 3
            Label = '(h) ';  n_k = 5;
        end
end

subplot(2,3,n_k)
% draw a 3D line graph using plot3
for m = 1:length(SEIS{1,1}.soc_thr) % soc
    plot3(SEIS{1,1}.soc_thr(m)*ones(length(Z_SEIS{1,1}{m}),1),...
          Z_SEIS{1,1}{m}(:,1),...
          Z_SEIS{1,1}{m}(:,2),colors_2{m},'LineWidth',2);
    hold on
end
for m = 1:length(SEIS{1,1}.soc_thr) % soc
    plot3(SEIS{1,1}.soc_thr(m)*ones(length(Z_DEIS{k,1}{m}),1),...
          Z_DEIS{k,1}{m}(:,1),...
          Z_DEIS{k,1}{m}(:,2),colors_1{m},'LineWidth',2);
end
xlabel('SOC',              'fontname','Times'); % X axis label
ylabel('{\itx} / {\itL}_4','fontname','Times'); % Y axis label
xlim([0 1])
switch change
    case 'cs_bar_neg'
Fig = 'Fig .5';
zlabel('{\itc}_{ss1} / {\itc}_{s1,max}','fontname','Times'); % Z axis label
ylim([0 0.456])
zlim([0 1])
    case 'cs_bar_pos'
Fig = 'Fig .5';
zlabel('{\itc}_{ss2} / {\itc}_{s2,max}','fontname','Times'); % Z axis label
ylim([0.5855 1])
zlim([0 1])
    case 'ce'
Fig = 'Fig .6';
zlabel('{\itc}_e / {\itc}_{e,0}','fontname','Times'); % Z axis label
ylim([0 1])
zlim([0.5 1.5])
    case 'phis_neg'
Fig = 'Fig .7';
zlabel('{\itφ}_{s1} (V)','fontname','Times'); % Z axis label
ylim([0 0.456])
zlim([0 0.4])
    case 'phis_pos'
Fig = 'Fig .7';
zlabel('{\itφ}_{s2} (V)','fontname','Times'); % Z axis label
ylim([0.5855 1])
zlim([3.6 4.4])
    case 'phie'
Fig = 'Fig .8';
zlabel('{\itφ}_e (V)','fontname','Times'); % Z axis label
ylim([0 1])
zlim([-0.1 0])
    case 'j_neg'
Fig = 'Fig .9';
zlabel('{\itj}_1 × {\itF} / {\iti}','fontname','Times'); % Z axis label
ylim([0 0.456])
zlim([0 0.025])
    case 'j_pos'
Fig = 'Fig .9';
zlabel('{\itj}_2 × {\itF} / {\iti}','fontname','Times'); % Z axis label
ylim([0.5855 1])
zlim([-0.025 0])
set(gca,'FontSize',18,'zDir','reverse')
end
set(gca,'FontSize',18,'yDir','reverse')
grid on; % show grid
title([Fig,Label,num2str(DEIS{k,1}.Crates),'C'],'fontsize',18,'fontweight','bold')

end

%% output origin
Origin{1,1} = ...
            [Z_SEIS{1,1}{1}(:,1) SEIS{1,1}.soc_thr(1)*ones(length(Z_SEIS{1,1}{1}),1) Z_SEIS{1,1}{1}(:,2)...
             Z_SEIS{1,1}{2}(:,1) SEIS{1,1}.soc_thr(2)*ones(length(Z_SEIS{1,1}{2}),1) Z_SEIS{1,1}{2}(:,2)...
             Z_SEIS{1,1}{3}(:,1) SEIS{1,1}.soc_thr(3)*ones(length(Z_SEIS{1,1}{3}),1) Z_SEIS{1,1}{3}(:,2)...
             Z_SEIS{1,1}{4}(:,1) SEIS{1,1}.soc_thr(4)*ones(length(Z_SEIS{1,1}{4}),1) Z_SEIS{1,1}{4}(:,2)...
             Z_SEIS{1,1}{5}(:,1) SEIS{1,1}.soc_thr(5)*ones(length(Z_SEIS{1,1}{5}),1) Z_SEIS{1,1}{5}(:,2)...
             Z_SEIS{1,1}{6}(:,1) SEIS{1,1}.soc_thr(5)*ones(length(Z_SEIS{1,1}{6}),1) Z_SEIS{1,1}{6}(:,2)...
             Z_SEIS{1,1}{7}(:,1) SEIS{1,1}.soc_thr(6)*ones(length(Z_SEIS{1,1}{7}),1) Z_SEIS{1,1}{7}(:,2)...
             Z_SEIS{1,1}{8}(:,1) SEIS{1,1}.soc_thr(7)*ones(length(Z_SEIS{1,1}{8}),1) Z_SEIS{1,1}{8}(:,2)...
             Z_SEIS{1,1}{9}(:,1) SEIS{1,1}.soc_thr(8)*ones(length(Z_SEIS{1,1}{9}),1) Z_SEIS{1,1}{9}(:,2)...
             Z_SEIS{1,1}{10}(:,1) SEIS{1,1}.soc_thr(9)*ones(length(Z_SEIS{1,1}{10}),1) Z_SEIS{1,1}{10}(:,2)];

for k = 1:size(DEIS,1) % C
Origin{k+1,1} = ...
            [Z_DEIS{k,1}{1}(:,1) DEIS{k,1}.soc_thr(1)*ones(length(Z_DEIS{k,1}{1}),1) Z_DEIS{k,1}{1}(:,2)...
             Z_DEIS{k,1}{2}(:,1) DEIS{k,1}.soc_thr(2)*ones(length(Z_DEIS{k,1}{2}),1) Z_DEIS{k,1}{2}(:,2)...
             Z_DEIS{k,1}{3}(:,1) DEIS{k,1}.soc_thr(3)*ones(length(Z_DEIS{k,1}{3}),1) Z_DEIS{k,1}{3}(:,2)...
             Z_DEIS{k,1}{4}(:,1) DEIS{k,1}.soc_thr(4)*ones(length(Z_DEIS{k,1}{4}),1) Z_DEIS{k,1}{4}(:,2)...
             Z_DEIS{k,1}{5}(:,1) DEIS{k,1}.soc_thr(5)*ones(length(Z_DEIS{k,1}{5}),1) Z_DEIS{k,1}{5}(:,2)...
             Z_DEIS{k,1}{6}(:,1) DEIS{k,1}.soc_thr(6)*ones(length(Z_DEIS{k,1}{6}),1) Z_DEIS{k,1}{6}(:,2)...
             Z_DEIS{k,1}{7}(:,1) DEIS{k,1}.soc_thr(7)*ones(length(Z_DEIS{k,1}{7}),1) Z_DEIS{k,1}{7}(:,2)...
             Z_DEIS{k,1}{8}(:,1) DEIS{k,1}.soc_thr(8)*ones(length(Z_DEIS{k,1}{8}),1) Z_DEIS{k,1}{8}(:,2)...
             Z_DEIS{k,1}{9}(:,1) DEIS{k,1}.soc_thr(9)*ones(length(Z_DEIS{k,1}{9}),1) Z_DEIS{k,1}{9}(:,2)...
             Z_DEIS{k,1}{10}(:,1) DEIS{k,1}.soc_thr(10)*ones(length(Z_DEIS{k,1}{10}),1) Z_DEIS{k,1}{10}(:,2)];
end

end
