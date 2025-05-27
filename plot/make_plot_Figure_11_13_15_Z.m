
function Origin = make_plot_Figure_11_13_15_Z(SEIS,DEIS,change)

Z_SEIS = option_input(SEIS,change);
Z_DEIS = option_input(DEIS,change);

colors_1 = {'g','r','b','m','k',...
            'g','r','b','m','k'};

colors_2 = {'og','or','ob','om','ok',...
            'og','or','ob','om','ok'};

%% OCV
switch change
    case 'Z_neg'
figure(11)
    case 'Z_pos'
figure(12)
    case 'Z_sep'
figure(13)
end
subplot(2,3,1)
% draw a 3D line graph using plot3
for m = 1:length(SEIS{1,1}.soc_thr) % soc
    plot3(Z_SEIS{1,1}{m}(:,1)*1e3,...
          SEIS{1,1}.soc_thr(m)*ones(length(Z_SEIS{1,1}{m}),1),...
          Z_SEIS{1,1}{m}(:,2)*1e3,colors_1{m},'LineWidth',2);
    hold on
end
ylabel('SOC');           % Y axis label
ylim([0 1])
switch change
    case 'Z_neg'
Fig = 'Fig .11';
xlabel('{\itZ}_1`   (mΩ·m^2)','fontname','Times'); % X axis label
zlabel('{\itZ}_1``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([0 1])
zlim([-1 0])
    case 'Z_pos'
Fig = 'Fig .13';
xlabel('{\itZ}_2`   (mΩ·m^2)','fontname','Times'); % X axis label
zlabel('{\itZ}_2``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([0 1])
zlim([-1 0])
    case 'Z_sep'
Fig = 'Fig .15';
xlabel('{\itZ}_3`   (mΩ·m^2)','fontname','Times'); % X axis label
zlabel('{\itZ}_3``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([0.02 0.1])
zlim([-0.08 0])
end
set(gca,'FontSize',18,'zDir','reverse')
grid on; % show grid
title([Fig,'(a) OCV'],'fontsize',18,'fontweight','bold')

%% different C-rates
for k = 1:size(DEIS,1) % C
        if k == 1
            Label = '(b) ';  n_k = 2;
        elseif k == 2
            Label = '(c) ';  n_k = 4;
        elseif k == 3
            Label = '(d) ';  n_k = 5;
        end

subplot(2,3,n_k)
% draw a 3D line graph using plot3
for m = 1:length(SEIS{1,1}.soc_thr) % soc
    plot3(Z_SEIS{1,1}{m}(:,1)*1e3,...
          SEIS{1,1}.soc_thr(m)*ones(length(Z_SEIS{1,1}{m}),1),...
          Z_SEIS{1,1}{m}(:,2)*1e3,colors_1{m},'LineWidth',2);
    hold on
end
for m = 1:length(SEIS{1,1}.soc_thr) % soc
    plot3(Z_DEIS{k,1}{m}(:,1)*1e3,...
          SEIS{1,1}.soc_thr(m)*ones(length(Z_DEIS{k,1}{m}),1),...
          Z_DEIS{k,1}{m}(:,2)*1e3,colors_2{m},'markersize',4);
end
ylabel('SOC');           % Y axis label
ylim([0 1])
switch change
    case 'Z_neg'
Fig = 'Fig .11';
xlabel('{\itZ}_1`   (mΩ·m^2)','fontname','Times'); % X axis label
zlabel('{\itZ}_1``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([0 1])
zlim([-1 0])
    case 'Z_pos'
Fig = 'Fig .13';
xlabel('{\itZ}_2`   (mΩ·m^2)','fontname','Times'); % X axis label
zlabel('{\itZ}_2``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([0 1])
zlim([-1 0])
    case 'Z_sep'
Fig = 'Fig .15';
xlabel('{\itZ}_3`   (mΩ·m^2)','fontname','Times'); % X axis label
zlabel('{\itZ}_3``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([0.02 0.1])
zlim([-0.08 0])
end
set(gca,'FontSize',18,'zDir','reverse')
grid on; % show grid
title([Fig,Label,num2str(DEIS{k,1}.Crates),'C'],'fontsize',18,'fontweight','bold')

end

%% output origin
Origin{1,1} = ...
            [SEIS{1,1}.soc_thr(1)*ones(length(Z_SEIS{1,1}{1}),1) Z_SEIS{1,1}{1}*1e3...
             SEIS{1,1}.soc_thr(2)*ones(length(Z_SEIS{1,1}{2}),1) Z_SEIS{1,1}{2}*1e3...
             SEIS{1,1}.soc_thr(3)*ones(length(Z_SEIS{1,1}{3}),1) Z_SEIS{1,1}{3}*1e3...
             SEIS{1,1}.soc_thr(4)*ones(length(Z_SEIS{1,1}{4}),1) Z_SEIS{1,1}{4}*1e3...
             SEIS{1,1}.soc_thr(5)*ones(length(Z_SEIS{1,1}{5}),1) Z_SEIS{1,1}{5}*1e3...
             SEIS{1,1}.soc_thr(6)*ones(length(Z_SEIS{1,1}{6}),1) Z_SEIS{1,1}{6}*1e3...
             SEIS{1,1}.soc_thr(7)*ones(length(Z_SEIS{1,1}{7}),1) Z_SEIS{1,1}{7}*1e3...
             SEIS{1,1}.soc_thr(8)*ones(length(Z_SEIS{1,1}{8}),1) Z_SEIS{1,1}{8}*1e3...
             SEIS{1,1}.soc_thr(9)*ones(length(Z_SEIS{1,1}{9}),1) Z_SEIS{1,1}{9}*1e3...
             SEIS{1,1}.soc_thr(10)*ones(length(Z_SEIS{1,1}{10}),1) Z_SEIS{1,1}{10}*1e3];

for k = 1:size(DEIS,1) % C
Origin{k+1,1} = ...
            [DEIS{k,1}.soc_thr(1)*ones(length(Z_DEIS{k,1}{1}),1) Z_DEIS{k,1}{1}*1e3...
             DEIS{k,1}.soc_thr(2)*ones(length(Z_DEIS{k,1}{2}),1) Z_DEIS{k,1}{2}*1e3...
             DEIS{k,1}.soc_thr(3)*ones(length(Z_DEIS{k,1}{3}),1) Z_DEIS{k,1}{3}*1e3...
             DEIS{k,1}.soc_thr(4)*ones(length(Z_DEIS{k,1}{4}),1) Z_DEIS{k,1}{4}*1e3...
             DEIS{k,1}.soc_thr(5)*ones(length(Z_DEIS{k,1}{5}),1) Z_DEIS{k,1}{5}*1e3...
             DEIS{k,1}.soc_thr(6)*ones(length(Z_DEIS{k,1}{6}),1) Z_DEIS{k,1}{6}*1e3...
             DEIS{k,1}.soc_thr(7)*ones(length(Z_DEIS{k,1}{7}),1) Z_DEIS{k,1}{7}*1e3...
             DEIS{k,1}.soc_thr(8)*ones(length(Z_DEIS{k,1}{8}),1) Z_DEIS{k,1}{8}*1e3...
             DEIS{k,1}.soc_thr(9)*ones(length(Z_DEIS{k,1}{9}),1) Z_DEIS{k,1}{9}*1e3...
             DEIS{k,1}.soc_thr(10)*ones(length(Z_DEIS{k,1}{10}),1) Z_DEIS{k,1}{10}*1e3];
end

end
