
function Origin = make_plot_Figure_18_state(SEIS,DEIS,Multiple)

%% plot Figure 18: cs, ce, φs, φe, j
Origin.cs_bar_neg = make_plot_vari(SEIS,DEIS,Multiple,6,'cs_bar_neg');
Origin.cs_bar_pos = make_plot_vari(SEIS,DEIS,Multiple,6,'cs_bar_pos');
Origin.ce         = make_plot_vari(SEIS,DEIS,Multiple,6,'ce');
Origin.phie       = make_plot_vari(SEIS,DEIS,Multiple,6,'phie');
Origin.phis_neg   = make_plot_vari(SEIS,DEIS,Multiple,6,'phis_neg');
Origin.phis_pos   = make_plot_vari(SEIS,DEIS,Multiple,6,'phis_pos');
Origin.j_neg      = make_plot_vari(SEIS,DEIS,Multiple,6,'j_neg');
Origin.j_pos      = make_plot_vari(SEIS,DEIS,Multiple,6,'j_pos');
Origin.electrode = [Origin.cs_bar_neg Origin.cs_bar_pos...
                    Origin.phis_neg   Origin.phis_pos...
                    Origin.j_neg      Origin.j_pos];
Origin.cell = [Origin.ce Origin.phie];

end

function out = make_plot_vari(SEIS,DEIS,Multiple,m,change)

Z_SEIS = option_input(SEIS,change);
Z_DEIS = option_input(DEIS,change);

% color
colors_1 = {'g','k','b','m','r',...
            'g','k','b','m','r'};

colors_2 = {'--g','--k','--b','--m','--r',...
            '--g','--k','--b','--m','--r'};

%% 4C
switch change
    case {'cs_bar_neg'}
        figure(2);        subplot(2,3,1);        Label = 'Fig. 18(a) {\itc}_{ss1}';
    case {'cs_bar_pos'}
        figure(2);        subplot(2,3,2);        Label = 'Fig. 18(b) {\itc}_{ss2}';
    case {'ce'}
        figure(2);        subplot(2,3,4);        Label = 'Fig. 18(c) {\itc}_e';
    case {'phie'}
        figure(2);        subplot(2,3,5);        Label = 'Fig. 18(d) {\itφ}_e';
    case {'phis_neg'}
        figure(3);        subplot(2,3,1);        Label = 'Fig. 18(e) {\itφ}_{s1}';
    case {'phis_pos'}
        figure(3);        subplot(2,3,2);        Label = 'Fig. 18(f) {\itφ}_{s2}';
    case {'j_neg'}
        figure(3);        subplot(2,3,4);        Label = 'Fig. 18(g) {\itj}_1';
    case {'j_pos'}
        figure(3);        subplot(2,3,5);        Label = 'Fig. 18(h) {\itj}_2';
end

for n = [1 3 5] % Multiple
    plot(Z_DEIS{1,n}{m}(:,1),Z_DEIS{1,n}{m}(:,2),colors_1{n},'LineWidth',2);
    hold on
end
plot(Z_SEIS{1,5}{m}(:,1),Z_SEIS{1,5}{m}(:,2),colors_2{2},'LineWidth',2);
% Add a legend.
h = legend(['{\itD}_e: ',num2str(Multiple(1))],...
           ['{\itD}_e: ',num2str(Multiple(3))],...
           ['{\itD}_e: ',num2str(Multiple(5))],'OCV','Location','southwest','FontSize',12);
xlabel('{\itx} / {\itL}_4','fontname','Times'); % Y axis label
switch change
    case 'cs_bar_neg'
ylabel('{\itc}_{ss1} / {\itc}_{s1,max}','fontname','Times'); % Z axis label
xlim([0 0.456])
ylim([0 1])
    case 'cs_bar_pos'
ylabel('{\itc}_{ss2} / {\itc}_{s2,max}','fontname','Times'); % Z axis label
xlim([0.5855 1])
ylim([0 1])
    case 'ce'
ylabel('{\itc}_e / {\itc}_{e,0}','fontname','Times'); % Z axis label
xlim([0 1])
ylim([0 2])
    case 'phis_neg'
ylabel('{\itφ}_{s1} (V)','fontname','Times'); % Z axis label
xlim([0 0.456])
ylim([0.0 0.3])
    case 'phis_pos'
ylabel('{\itφ}_{s2} (V)','fontname','Times'); % Z axis label
xlim([0.5855 1])
ylim([3.7 4.0])
    case 'phie'
ylabel('{\itφ}_e (V)','fontname','Times'); % Z axis label
xlim([0 1])
ylim([-0.2 0.02])
    case 'j_neg'
ylabel('{\itj}_1 × {\itF} / {\iti}','fontname','Times'); % Z axis label
xlim([0 0.456])
ylim([-0.0025 0.025])
    case 'j_pos'
ylabel('{\itj}_2 × {\itF} / {\iti}','fontname','Times'); % Z axis label
xlim([0.5855 1])
ylim([-0.025 0.0025])
set(gca,'FontSize',18,'yDir','reverse')
end
set(gca,'FontSize',18)
grid on; % show grid
title(Label,'fontsize',18,'fontweight','bold')

%% output origin
out = [Z_DEIS{1,1}{m}...
       Z_DEIS{1,2}{m}(:,2)...
       Z_DEIS{1,3}{m}(:,2)...
       Z_DEIS{1,4}{m}(:,2)...
       Z_DEIS{1,5}{m}(:,2)...
       Z_SEIS{1,5}{m}(:,2)];

end
