
function Origin = make_plot_Figure_19_Z(SEIS,DEIS,Multiple)

%% plot Figure 19: Z1/Z2/Z3/Z4 and ΔZ1/ΔZ2/ΔZ3/ΔZ4
Origin.Z_neg  = make_plot_Z(SEIS,DEIS,Multiple,6,'Z_neg');
Origin.Z_pos  = make_plot_Z(SEIS,DEIS,Multiple,6,'Z_pos');
Origin.Z_sep  = make_plot_Z(SEIS,DEIS,Multiple,6,'Z_sep');
Origin.Z_cell = make_plot_Z(SEIS,DEIS,Multiple,6,'Z_cell');
Origin.Z_tot  = [Origin.Z_neg Origin.Z_pos Origin.Z_sep Origin.Z_cell];

Origin.delta_Mag_neg  = make_plot_diff(SEIS,DEIS,Multiple,6,'delta_Mag_neg');
Origin.delta_Mag_pos  = make_plot_diff(SEIS,DEIS,Multiple,6,'delta_Mag_pos');
Origin.delta_Mag_sep  = make_plot_diff(SEIS,DEIS,Multiple,6,'delta_Mag_sep');
Origin.delta_Mag_cell = make_plot_diff(SEIS,DEIS,Multiple,6,'delta_Mag_cell');
Origin.delta_Mag_tot  = [Origin.delta_Mag_neg Origin.delta_Mag_pos Origin.delta_Mag_sep Origin.delta_Mag_cell];
end

function out = make_plot_Z(SEIS,DEIS,Multiple,m,change)

Z_SEIS = option_input(SEIS,change);
Z_DEIS = option_input(DEIS,change);

% color
colors_1 = {'g','k','b','m','r',...
            'g','k','b','m','r'};

colors_2 = {'og','ok','ob','om','or',...
            'og','ok','ob','om','or'};

%% 4C
switch change
    case {'Z_neg'}
        figure(4);        subplot(2,3,1);        Label = 'Fig .19(a) {\itZ}_1';
    case {'Z_pos'}
        figure(4);        subplot(2,3,4);        Label = 'Fig .19(c) {\itZ}_2';
    case {'Z_sep'}
        figure(5);        subplot(2,3,1);        Label = 'Fig .19(e) {\itZ}_3';
    case {'Z_cell'}
        figure(5);        subplot(2,3,4);        Label = 'Fig .19(g) {\itZ}_4';
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
    case 'Z_neg'
xlabel('{\itZ}_1`   (mΩ·m^2)','fontname','Times'); % X axis label
ylabel('{\itZ}_1``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([-0.1 1])
ylim([-0.5 0.05])
    case 'Z_pos'
xlabel('{\itZ}_2`   (mΩ·m^2)','fontname','Times'); % X axis label
ylabel('{\itZ}_2``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([-0.1 1])
ylim([-0.5 0.05])
    case 'Z_sep'
xlabel('{\itZ}_3`   (mΩ·m^2)','fontname','Times'); % X axis label
ylabel('{\itZ}_3``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([-0.02 0.2])
ylim([-0.1 0.01])
    case 'Z_cell'
xlabel('{\itZ}_4`   (mΩ·m^2)','fontname','Times'); % X axis label
ylabel('{\itZ}_4``  (mΩ·m^2)','fontname','Times'); % Z axis label
xlim([-0.2 2])
ylim([-1 0.1])
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

function out = make_plot_diff(SEIS,DEIS,Multiple,m,change)

for k = 1:size(DEIS,1)                       % C
    for n = 1:size(DEIS,2)                   % multiple
        switch change
            case 'delta_Mag_neg'
                delta_Z{k,n}{m} = (abs(DEIS{k,n}.Z{m}.DFN_neg)...
                                 - abs(SEIS{1,n}.Z{m}.DFN_neg))...
                                ./ abs(SEIS{1,n}.Z{m}.DFN_neg) * 100;
            case 'delta_Mag_pos'
                delta_Z{k,n}{m} = (abs(DEIS{k,n}.Z{m}.DFN_pos)...
                                 - abs(SEIS{1,n}.Z{m}.DFN_pos))...
                                ./ abs(SEIS{1,n}.Z{m}.DFN_pos) * 100;
            case 'delta_Mag_sep'
                delta_Z{k,n}{m} = (abs(DEIS{k,n}.Z{m}.DFN_sep)...
                                 - abs(SEIS{1,n}.Z{m}.DFN_sep))...
                                ./ abs(SEIS{1,n}.Z{m}.DFN_sep) * 100;
            case 'delta_Mag_cell'
                delta_Z{k,n}{m} = (abs(DEIS{k,n}.Z{m}.DFN_cell)...
                                 - abs(SEIS{1,n}.Z{m}.DFN_cell))...
                                ./ abs(SEIS{1,n}.Z{m}.DFN_cell) * 100;
        end
    end
end

% color
colors_1 = {'g','k','b','m','r',...
            'g','k','b','m','r'};

%% different C-rates
switch change
    case {'delta_Mag_neg'}
        figure(4);        subplot(2,3,2);        Label = 'Fig .19(b) Δ|{\itZ}_1|';
    case {'delta_Mag_pos'}
        figure(4);        subplot(2,3,5);        Label = 'Fig .19(d) Δ|{\itZ}_2|';
    case {'delta_Mag_sep'}
        figure(5);        subplot(2,3,2);        Label = 'Fig .19(f) Δ|{\itZ}_3|';
    case {'delta_Mag_cell'}
        figure(5);        subplot(2,3,5);        Label = 'Fig .19(h) Δ|{\itZ}_4|';
end

for n = [1 3 5] % multiple
    semilogx(SEIS{1,1}.f,delta_Z{1,n}{6},colors_1{n},'LineWidth',2)
    hold on
end
% Add a legend.
h = legend(['{\itD}_e: ',num2str(Multiple(1))],...
           ['{\itD}_e: ',num2str(Multiple(3))],...
           ['{\itD}_e: ',num2str(Multiple(5))],'Location','northeast');
% Set the properties of the coordinate axis.
xlabel('{\itf} (Hz)','fontsize',15,'fontname','Times')
switch change
    case 'delta_Mag_neg'
ylabel('Δ|{\itZ}_1| (%)','fontsize',15,'fontname','Times')
axis([10^-3 10^4  -10 30])
    case 'delta_Mag_pos'
ylabel('Δ|{\itZ}_2| (%)','fontsize',15,'fontname','Times')
axis([10^-3 10^4  -20 60])
    case 'delta_Mag_sep'
ylabel('Δ|{\itZ}_3| (%)','fontsize',15,'fontname','Times')
axis([10^-3 10^4 -10 10])
    case 'delta_Mag_cell'
ylabel('Δ|{\itZ}_4| (%)','fontsize',15,'fontname','Times')
axis([10^-3 10^4  -10 30])
end
axis([10^-3 10^4  -20 60])
% Set grid.
set(gca,'FontSize',18,'xgrid','on')
set(gca,'FontSize',18,'ygrid','on')
% The logspace function is used here to generate logarithmic scale values
xticks(logspace(-3, 3, 4));
title(Label,'fontsize',18,'fontweight','bold')

%% output origin
out = [SEIS{1,1}.f delta_Z{1,1}{6}...
                   delta_Z{1,2}{6}...
                   delta_Z{1,3}{6}...
                   delta_Z{1,4}{6}...
                   delta_Z{1,5}{6}];

end
