
function Origin = make_plot_Figure_12_14_16_diffZ(SEIS,DEIS,Crates,diff)

% color
colors{1} = {'c','c--','c:','c-.'};
colors{2} = {'r','r--','r:','r-.'};
colors{3} = {'b','b--','b:','b-.'};
colors{4} = {'m','m--','m:','m-.'};
colors{5} = {'g','g--','g:','g-.'};
colors{6} = {'k','k--','k:','k-.'};

for k = 1:size(DEIS,1)                   % C
    for m = 1:length(DEIS{k,1}.soc_thr)  % soc
        switch diff
            case 'delta_Mag_neg'
                delta_Z{k,1}{m} = (abs(DEIS{k,1}.Z{m}.DFN_neg)...
                                 - abs(SEIS{1,1}.Z{m}.DFN_neg))...
                                ./ abs(SEIS{1,1}.Z{m}.DFN_neg) * 100;
            case 'delta_Mag_pos'
                delta_Z{k,1}{m} = (abs(DEIS{k,1}.Z{m}.DFN_pos)...
                                 - abs(SEIS{1,1}.Z{m}.DFN_pos))...
                                ./ abs(SEIS{1,1}.Z{m}.DFN_pos) * 100;
            case 'delta_Mag_sep'
                delta_Z{k,1}{m} = (abs(DEIS{k,1}.Z{m}.DFN_sep)...
                                 - abs(SEIS{1,1}.Z{m}.DFN_sep))...
                                ./ abs(SEIS{1,1}.Z{m}.DFN_sep) * 100;
        end
    end
end

for m = 1:length(SEIS{1,1}.soc_thr)  % soc
    make_plot_fre_diff(SEIS{1,1}.f,...
                       SEIS{1,1}.soc_thr,...
                       delta_Z,Crates,colors,diff,m);
end

%% output origin
for k = 1:size(DEIS,1)                  % C
    for m = 1:length(DEIS{k,1}.soc_thr) % soc
        Origin.delta_Z_soc{m}(:,2*k-1:2*k) = [SEIS{1,1}.f delta_Z{k,1}{m}];
    end
end

Origin.delta_Z = [Origin.delta_Z_soc{1} Origin.delta_Z_soc{2} Origin.delta_Z_soc{3}...
                  Origin.delta_Z_soc{4} Origin.delta_Z_soc{5} Origin.delta_Z_soc{6}...
                  Origin.delta_Z_soc{7} Origin.delta_Z_soc{8} Origin.delta_Z_soc{9}...
                  Origin.delta_Z_soc{10}];

end

function make_plot_fre_diff(f,soc_thr,delta_Z,Crates,colors,diff,m)

%% plot
switch diff
    case 'delta_Mag_neg'
k_figure = 13+0;
    case 'delta_Mag_pos'
k_figure = 13+2;
    case 'delta_Mag_sep'
k_figure = 13+4;
end
if m == 1
    figure(k_figure+1);  subplot(3,3,1);  Label = '(a) SOC=';
elseif m == 2
    figure(k_figure+1);  subplot(3,3,4);  Label = '(b) SOC=';
elseif m == 3
    figure(k_figure+1);  subplot(3,3,7);  Label = '(c) SOC=';
elseif m == 4
    figure(k_figure+2);  subplot(3,3,1);  Label = '(d) SOC=';
elseif m == 5
    figure(k_figure+2);  subplot(3,3,4);  Label = '(e) SOC=';
elseif m == 6
    figure(k_figure+1);  subplot(3,3,2);  Label = '(f) SOC=';
elseif m == 7
    figure(k_figure+1);  subplot(3,3,5);  Label = '(g) SOC=';
elseif m == 8
    figure(k_figure+1);  subplot(3,3,8);  Label = '(h) SOC=';
elseif m == 9
    figure(k_figure+2);  subplot(3,3,2);  Label = '(i) SOC=';
elseif m == 10
    figure(k_figure+2);  subplot(3,3,5);  Label = '(j) SOC=';
end
for k = 1:size(delta_Z,1)  % C
    semilogx(f,delta_Z{k,1}{m},colors{k}{1},'LineWidth',2)
    hold on
end
% Add a legend.
h = legend([num2str(Crates(1)),' C'],...
           [num2str(Crates(2)),' C'],...
           [num2str(Crates(3)),' C'],'Location','northeast');
% Set the properties of the coordinate axis.
xlabel('{\itf} (Hz)','fontsize',15,'fontname','Times')
switch diff
    case 'delta_Mag_neg'
Fig = 'Fig .12';
ylabel('Δ|{\itZ}_1| (%)','fontsize',15,'fontname','Times')
if m == 1 || m == 2 || m == 3 || m == 4 || m == 5
    axis([10^-3 10^4 -10 10])
else 
    axis([10^-3 10^4  -4 40])
end
    case 'delta_Mag_pos'
Fig = 'Fig .14';
ylabel('Δ|{\itZ}_2| (%)','fontsize',15,'fontname','Times')
if m == 1 || m == 2 || m == 3 || m == 4 || m == 5
    axis([10^-3 10^4 -20 20])
else 
    axis([10^-3 10^4  -4 40])
end
    case 'delta_Mag_sep'
Fig = 'Fig .16';
ylabel('Δ|{\itZ}_3| (%)','fontsize',15,'fontname','Times')
axis([10^-3 10^4 -5 5])
end
% Set grid.
set(gca,'FontSize',18,'xgrid','on')
set(gca,'FontSize',18,'ygrid','on')
% The logspace function is used here to generate logarithmic scale values
xticks(logspace(-3, 3, 4));
title([Fig,Label,num2str(soc_thr(m))],'fontsize',18,'fontweight','bold')

end
