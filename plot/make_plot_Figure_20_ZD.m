
function Origin = make_plot_Figure_20_ZD(SEIS,DEIS,Multiple)

Origin.ZDs_neg  = make_plot_ZD(SEIS,DEIS,Multiple,6,'ZDs_neg');
Origin.ZDe_neg  = make_plot_ZD(SEIS,DEIS,Multiple,6,'ZDe_neg');
Origin.ZDs_pos  = make_plot_ZD(SEIS,DEIS,Multiple,6,'ZDs_pos');
Origin.ZDe_pos  = make_plot_ZD(SEIS,DEIS,Multiple,6,'ZDe_pos');
Origin.ZD_tot   = [Origin.ZDs_neg Origin.ZDe_neg Origin.ZDs_pos Origin.ZDe_pos];

f_ZD = SEIS{1,1}.f(1:length(SEIS{1,1}.Z{1}.Ds_neg));

Origin.ZDs_neg_re  = make_plot_ZD_re(SEIS,DEIS,Multiple,f_ZD,6,'ZDs_neg');
Origin.ZDe_neg_re  = make_plot_ZD_re(SEIS,DEIS,Multiple,f_ZD,6,'ZDe_neg');
Origin.ZDs_pos_re  = make_plot_ZD_re(SEIS,DEIS,Multiple,f_ZD,6,'ZDs_pos');
Origin.ZDe_pos_re  = make_plot_ZD_re(SEIS,DEIS,Multiple,f_ZD,6,'ZDe_pos');
Origin.ZD_re_tot   = [Origin.ZDs_neg_re Origin.ZDe_neg_re Origin.ZDs_pos_re Origin.ZDe_pos_re];

Origin.ZDs_neg_im  = make_plot_ZD_im(SEIS,DEIS,Multiple,f_ZD,6,'ZDs_neg');
Origin.ZDe_neg_im  = make_plot_ZD_im(SEIS,DEIS,Multiple,f_ZD,6,'ZDe_neg');
Origin.ZDs_pos_im  = make_plot_ZD_im(SEIS,DEIS,Multiple,f_ZD,6,'ZDs_pos');
Origin.ZDe_pos_im  = make_plot_ZD_im(SEIS,DEIS,Multiple,f_ZD,6,'ZDe_pos');
Origin.ZD_im_tot   = [Origin.ZDs_neg_im Origin.ZDe_neg_im Origin.ZDs_pos_im Origin.ZDe_pos_im];

function out = make_plot_ZD(SEIS,DEIS,Multiple,m,change)

Z_SEIS = option_input(SEIS,change);
Z_DEIS = option_input(DEIS,change);

% color
colors_1 = {'g','k','b','m','r',...
            'g','k','b','m','r'};

colors_2 = {'og','ok','ob','om','or',...
            'og','ok','ob','om','or'};

%% 4C
switch change
    case {'ZDs_neg'}
        figure(6);        subplot(2,3,1);        Label = 'Fig .20(a) {\itZ}_{Ds1}';
    case {'ZDe_neg'}
        figure(6);        subplot(2,3,4);        Label = 'Fig .20(c) {\itZ}_{De1}';
    case {'ZDs_pos'}
        figure(7);        subplot(2,3,1);        Label = 'Fig .20(e) {\itZ}_{Ds2}';
    case {'ZDe_pos'}
        figure(7);        subplot(2,3,4);        Label = 'Fig .20(g) {\itZ}_{De2}';
end

for n = [1 3 5] % Multiple
    plot(Z_SEIS{1,n}{m}(:,1)*1e3,Z_SEIS{1,n}{m}(:,2)*1e3,colors_1{n},'LineWidth',2);
    hold on
end
for n = [1 3 5] % Multiple
    plot(Z_DEIS{1,n}{m}(:,1)*1e3,Z_DEIS{1,n}{m}(:,2)*1e3,colors_2{n},'markersize',4);
end
% Add a legend.
h = legend(['{\itD}_e: ',num2str(Multiple(1))],...
           ['{\itD}_e: ',num2str(Multiple(3))],...
           ['{\itD}_e: ',num2str(Multiple(5))],'DEIS','Location','northwest','FontSize',12);
switch change
    case 'ZDs_neg'
xlabel('{\itZ}_{Ds1}`   (mΩ·m^2)','fontname','Times'); % X axis label
ylabel('{\itZ}_{Ds1}``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([-0.03 0.3])
ylim([-0.3 0.03])
    case 'ZDe_neg'
xlabel('{\itZ}_{De1}`   (mΩ·m^2)','fontname','Times'); % X axis label
ylabel('{\itZ}_{De1}``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([-0.03 0.3])
ylim([-0.3 0.03])
    case 'ZDs_pos'
xlabel('{\itZ}_{Ds2}`   (mΩ·m^2)','fontname','Times'); % X axis label
ylabel('{\itZ}_{Ds2}``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([-0.06 0.6])
ylim([-0.6 0.06])
    case 'ZDe_pos'
xlabel('{\itZ}_{De2}`   (mΩ·m^2)','fontname','Times'); % X axis label
ylabel('{\itZ}_{De2}``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([-0.03 0.3])
ylim([-0.3 0.03])
end
% Set grid.
set(gca,'FontSize',18,'yDir','reverse')
grid on; % show grid
title(Label,'fontsize',18,'fontweight','bold')

%% output origin
out{2,1} = [Z_SEIS{1,1}{m}*1e3...
            Z_SEIS{1,2}{m}*1e3...
            Z_SEIS{1,3}{m}*1e3...
            Z_SEIS{1,4}{m}*1e3...
            Z_SEIS{1,5}{m}*1e3];
    
out{1,1} = [Z_DEIS{1,1}{m}*1e3...
            Z_DEIS{1,2}{m}*1e3...
            Z_DEIS{1,3}{m}*1e3...
            Z_DEIS{1,4}{m}*1e3...
            Z_DEIS{1,5}{m}*1e3];

out = [out{1,1} out{2,1}];

end

function out = make_plot_ZD_re(SEIS,DEIS,Multiple,f_ZD,m,change)

Z_SEIS = option_input(SEIS,change);
Z_DEIS = option_input(DEIS,change);

% color
colors_1 = {'g','k','b','m','r',...
            'g','k','b','m','r'};

colors_2 = {'og','ok','ob','om','or',...
            'og','ok','ob','om','or'};

%% 4C
switch change
    case {'ZDs_neg'}
        figure(6);        subplot(2,3,2);        Label = 'Fig .20(b top) {\itZ}_{Ds1}`';
    case {'ZDe_neg'}
        figure(6);        subplot(2,3,5);        Label = 'Fig .20(d top) {\itZ}_{De1}`';
    case {'ZDs_pos'}
        figure(7);        subplot(2,3,2);        Label = 'Fig .20(f top) {\itZ}_{Ds2}`';
    case {'ZDe_pos'}
        figure(7);        subplot(2,3,5);        Label = 'Fig .20(h top) {\itZ}_{De2}`';
end

for n = [1 3 5] % Multiple
    semilogx(f_ZD,Z_SEIS{1,n}{m}(:,1)*1e3,colors_1{n},'LineWidth',2);
    hold on
end
for n = [1 3 5] % Multiple
    semilogx(f_ZD,Z_DEIS{1,n}{m}(:,1)*1e3,colors_2{n},'markersize',4);
end
% Add a legend.
h = legend(['{\itD}_e: ',num2str(Multiple(1))],...
           ['{\itD}_e: ',num2str(Multiple(3))],...
           ['{\itD}_e: ',num2str(Multiple(5))],'DEIS','Location','northwest','FontSize',12);
switch change
    case 'ZDs_neg'
xlabel('{\itf} (Hz)','fontname','Times'); % X axis label
ylabel('{\itZ}_{Ds1}`  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([10^-3 10^0])
ylim([0 0.3])
    case 'ZDe_neg'
xlabel('{\itf} (Hz)','fontname','Times'); % X axis label
ylabel('{\itZ}_{De1}`  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([10^-3 10^0])
ylim([0 0.3])
    case 'ZDs_pos'
xlabel('{\itf} (Hz)','fontname','Times'); % X axis label
ylabel('{\itZ}_{Ds2}`  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([10^-3 10^0])
ylim([0 0.6])
    case 'ZDe_pos'
xlabel('{\itf} (Hz)','fontname','Times'); % X axis label
ylabel('{\itZ}_{De2}`  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([10^-3 10^0])
ylim([0 0.3])
end
% Set grid.
set(gca,'FontSize',18)
grid on; % show grid
% The logspace function is used here to generate logarithmic scale values
xticks(logspace(-3, 0, 4));
title(Label,'fontsize',18,'fontweight','bold')

%% output origin
out{2,1} = [f_ZD,Z_SEIS{1,1}{m}(:,1)*1e3...
                 Z_SEIS{1,2}{m}(:,1)*1e3...
                 Z_SEIS{1,3}{m}(:,1)*1e3...
                 Z_SEIS{1,4}{m}(:,1)*1e3...
                 Z_SEIS{1,5}{m}(:,1)*1e3];
    
out{1,1} = [f_ZD,Z_DEIS{1,1}{m}(:,1)*1e3...
                 Z_DEIS{1,2}{m}(:,1)*1e3...
                 Z_DEIS{1,3}{m}(:,1)*1e3...
                 Z_DEIS{1,4}{m}(:,1)*1e3...
                 Z_DEIS{1,5}{m}(:,1)*1e3];

out = [out{1,1} out{2,1}];

end

function out = make_plot_ZD_im(SEIS,DEIS,Multiple,f_ZD,m,change)

Z_SEIS = option_input(SEIS,change);
Z_DEIS = option_input(DEIS,change);

% color
colors_1 = {'g','k','b','m','r',...
            'g','k','b','m','r'};

colors_2 = {'og','ok','ob','om','or',...
            'og','ok','ob','om','or'};

%% 4C
switch change
    case {'ZDs_neg'}
        figure(6);        subplot(2,3,3);        Label = 'Fig .20(b bottom) {\itZ}_{Ds1}``';
    case {'ZDe_neg'}
        figure(6);        subplot(2,3,6);        Label = 'Fig .20(d bottom) {\itZ}_{De1}``';
    case {'ZDs_pos'}
        figure(7);        subplot(2,3,3);        Label = 'Fig .20(f bottom) {\itZ}_{Ds2}``';
    case {'ZDe_pos'}
        figure(7);        subplot(2,3,6);        Label = 'Fig .20(h bottom) {\itZ}_{De2}``';
end

for n = [1 3 5] % Multiple
    semilogx(f_ZD,Z_SEIS{1,n}{m}(:,2)*1e3,colors_1{n},'LineWidth',2);
    hold on
end
for n = [1 3 5] % Multiple
    semilogx(f_ZD,Z_DEIS{1,n}{m}(:,2)*1e3,colors_2{n},'markersize',4);
end
% Add a legend.
h = legend(['{\itD}_e: ',num2str(Multiple(1))],...
           ['{\itD}_e: ',num2str(Multiple(3))],...
           ['{\itD}_e: ',num2str(Multiple(5))],'DEIS','Location','northwest','FontSize',12);
switch change
    case 'ZDs_neg'
xlabel('{\itf} (Hz)','fontname','Times'); % X axis label
ylabel('{\itZ}_{Ds1}``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([10^-3 10^0])
ylim([-0.3 0])
    case 'ZDe_neg'
xlabel('{\itf} (Hz)','fontname','Times'); % X axis label
ylabel('{\itZ}_{De1}``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([10^-3 10^0])
ylim([-0.3 0])
    case 'ZDs_pos'
xlabel('{\itf} (Hz)','fontname','Times'); % X axis label
ylabel('{\itZ}_{Ds2}``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([10^-3 10^0])
ylim([-0.6 0])
    case 'ZDe_pos'
xlabel('{\itf} (Hz)','fontname','Times'); % X axis label
ylabel('{\itZ}_{De2}``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([10^-3 10^0])
ylim([-0.3 0])
end
% Set grid.
set(gca,'FontSize',18,'yDir','reverse')
grid on; % show grid
% The logspace function is used here to generate logarithmic scale values
xticks(logspace(-3, 0, 4));
title(Label,'fontsize',18,'fontweight','bold')

%% output origin
out{2,1} = [f_ZD,Z_SEIS{1,1}{m}(:,2)*1e3...
                 Z_SEIS{1,2}{m}(:,2)*1e3...
                 Z_SEIS{1,3}{m}(:,2)*1e3...
                 Z_SEIS{1,4}{m}(:,2)*1e3...
                 Z_SEIS{1,5}{m}(:,2)*1e3];
    
out{1,1} = [f_ZD,Z_DEIS{1,1}{m}(:,2)*1e3...
                 Z_DEIS{1,2}{m}(:,2)*1e3...
                 Z_DEIS{1,3}{m}(:,2)*1e3...
                 Z_DEIS{1,4}{m}(:,2)*1e3...
                 Z_DEIS{1,5}{m}(:,2)*1e3];

out = [out{1,1} out{2,1}];

end

end