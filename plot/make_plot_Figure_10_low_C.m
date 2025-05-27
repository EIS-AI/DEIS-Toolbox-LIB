
function Origin = make_plot_Figure_10_low_C(SEIS,volt,DEIS,Crates)

% color
colors{1} = {'c','c--','c:','c-.'};
colors{2} = {'r','r--','r:','r-.'};
colors{3} = {'b','b--','b:','b-.'};
colors{4} = {'m','m--','m:','m-.'};
colors{5} = {'g','g--','g:','g-.'};
colors{6} = {'k','k--','k:','k-.'};

% voltage relative difference
for k = 1:size(DEIS,1)     % C
    for n = 1:size(DEIS,2) % multiple
        % stoichiometric coefficient for porous electrode [-]
        w{k,n} = volt{k,n}.p.NP.s0_neg   + volt{k,n}.soc * ...
                (volt{k,n}.p.NP.s100_neg - volt{k,n}.p.NP.s0_neg);
        z{k,n} = volt{k,n}.p.NP.s0_pos   + volt{k,n}.soc * ...
                (volt{k,n}.p.NP.s100_pos - volt{k,n}.p.NP.s0_pos);
        % the equilibrium potential of the full cell [V]
        volt{k,n}.Eq_cell = volt{k,n}.p.U_pos(z{k,n})...
                          - volt{k,n}.p.U_neg(w{k,n});
        % relative difference [-]
        delta{k,n}.v_cell = (volt{k,n}.v_cell - volt{k,n}.Eq_cell)...
                                             ./ volt{k,n}.Eq_cell * 100;
    end
end

% impedance magnitude relative difference
for k = 1:size(DEIS,1)                       % C
    for n = 1:size(DEIS,2)                   % multiple
        for m = 1:length(DEIS{k,n}.soc_thr)  % soc
        DEIS{k,n}.delta{m}.DFN_cell = (abs(DEIS{k,n}.Z{m}.DFN_cell)...
                                     - abs(SEIS{1,n}.Z{m}.DFN_cell))...
                                    ./ abs(SEIS{1,n}.Z{m}.DFN_cell) * 100;
        end
    end
end

%% cell voltage
figure(1)
subplot(2,3,1)
for k = 1:size(Crates,1)  % C
    plot(volt{k,1}.soc,...
         volt{k,1}.v_cell,colors{k+2}{1},'LineWidth',2)
    hold on
end
plot(volt{1,1}.soc,...
     volt{1,1}.Eq_cell,colors{6}{1},'LineWidth',2)
% Add a legend.
h = legend([num2str(Crates(1)),' C'],...
           [num2str(Crates(2)),' C'],...
           [num2str(Crates(3)),' C'],'OCV','Location','southwest');
% Set the properties of the coordinate axis.
xlabel('SOC','fontsize',24,'fontname','Times')
ylabel('{\itv}_4 (V)','fontsize',24,'fontname','Times')
axis([-0.1 1.1 2.7 4.2])
yticks([2.8:0.4:4.0])
% Flip the X axis up and down.
set(gca, 'XDir','reverse')
% Set grid.
set(gca,'FontSize',18,'xgrid','on')
set(gca,'FontSize',18,'ygrid','on')
title('Fig .10(a) {\itv}_4','fontsize',18,'fontweight','bold')

%% cell voltage relative difference
subplot(2,3,2)
for k = 1:size(Crates,1)  % C
    plot(volt{k,1}.soc,...
         delta{k,1}.v_cell,colors{k+2}{1},'LineWidth',2)
    hold on
end
% Add a legend.
h = legend([num2str(Crates(1)),' C'],...
           [num2str(Crates(2)),' C'],...
           [num2str(Crates(3)),' C'],'Location','northwest');
% Set the properties of the coordinate axis.
xlabel('SOC','fontsize',24,'fontname','Times')
ylabel('Δ{\itv}_4 (%)','fontsize',24,'fontname','Times')
axis([-0.1 1.1 -15 1.5])
% Flip the axis up and down.
set(gca, 'XDir','reverse')
set(gca, 'YDir','reverse')
% Set grid.
set(gca,'FontSize',18,'xgrid','on')
set(gca,'FontSize',18,'ygrid','on')
title('Fig .10(b) Δ{\itv}_4','fontsize',18,'fontweight','bold')

%% cell impedance
mm = 8;  % soc
subplot(2,3,4)
for k = 1:size(Crates,1)  % C
    plot(real(DEIS{k,1}.Z{mm}.DFN_cell)*1e3,...
         imag(DEIS{k,1}.Z{mm}.DFN_cell)*1e3,colors{k+2}{1},'LineWidth',2)
    hold on
end
plot(real(SEIS{1,1}.Z{mm}.DFN_cell)*1e3,...
     imag(SEIS{1,1}.Z{mm}.DFN_cell)*1e3,colors{6}{1},'LineWidth',2)
% Add a legend.
h = legend([num2str(Crates(1)),' C'],...
           [num2str(Crates(2)),' C'],...
           [num2str(Crates(3)),' C'],'OCV','Location','northwest');
% Set the properties of the coordinate axis.
xlabel('{\itZ}_4`  (mΩ·m^2)','fontsize',15,'fontname','Times')
ylabel('{\itZ}_4`` (mΩ·m^2)','fontsize',15,'fontname','Times')
axis([-0.06 1.2 -1.2 0.06])
xticks([ 0.0:0.3:1.2])
yticks([-1.2:0.3:0.0])
% Flip the Y axis up and down.
set(gca, 'YDir','reverse')
% Set grid.
set(gca,'FontSize',18,'xgrid','on')
set(gca,'FontSize',18,'ygrid','on')
title(['Fig .10(c) SOC=',num2str(SEIS{1,1}.soc_thr(mm))],'fontsize',18,'fontweight','bold')

%% cell impedance relative difference
subplot(2,3,5)
for k = 1:size(Crates,1)  % C
    semilogx(DEIS{k,1}.f,...
             DEIS{k,1}.delta{mm}.DFN_cell,colors{k+2}{1},'LineWidth',2)
    hold on
end
% Add a legend.
h = legend([num2str(Crates(1)),' C'],...
           [num2str(Crates(2)),' C'],...
           [num2str(Crates(3)),' C'],'Location','northeast');
% Set the properties of the coordinate axis.
xlabel('{\itf} (Hz)','fontsize',15,'fontname','Times')
ylabel('Δ|{\itZ}_4| (%)','fontsize',15,'fontname','Times')
axis([10^-3 10^4 -0.5 5])
% Set grid.
set(gca,'FontSize',18,'xgrid','on')
set(gca,'FontSize',18,'ygrid','on')
% The logspace function is used here to generate logarithmic scale values
xticks(logspace(-3, 3, 4));
title(['Fig .10(d) SOC=',num2str(SEIS{1,1}.soc_thr(mm))],'fontsize',18,'fontweight','bold')

%% output origin
Origin.v_soc{1,1} = [volt{1,1}.soc'...
                     volt{1,1}.Eq_cell'];
Origin.Z_cell = [real(SEIS{1,1}.Z{mm}.DFN_cell)...
                 imag(SEIS{1,1}.Z{mm}.DFN_cell)]*1e3;

for k = 1:size(Crates,1)  % C
    Origin.v_soc{1,k+1} = [volt{k,1}.soc' volt{k,1}.Eq_cell'];
    Origin.v_soc{2,k+1} = [volt{k,1}.soc' delta{k,1}.v_cell'];

    Origin.Z_cell(:,2*k+1:2*k+2) = [real(DEIS{k,1}.Z{mm}.DFN_cell)...
                                    imag(DEIS{k,1}.Z{mm}.DFN_cell)]*1e3;
    Origin.Z_cell(:,2*k+7:2*k+8) = [DEIS{k,1}.f...
                                    DEIS{k,1}.delta{mm}.DFN_cell];
end

end
