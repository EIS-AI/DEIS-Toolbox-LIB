
function Origin = make_plot_Figure_4_v(volt,Crates)

% color
colors{1} = {'c','c--','c:','c-.'};
colors{2} = {'r','r--','r:','r-.'};
colors{3} = {'b','b--','b:','b-.'};
colors{4} = {'m','m--','m:','m-.'};
colors{5} = {'g','g--','g:','g-.'};
colors{6} = {'k','k--','k:','k-.'};

% voltage absolute difference
for k = 1:size(volt,1)     % C
    for n = 1:size(volt,2) % multiple
        % stoichiometric coefficient for porous electrode [-]
        w{k,n} = volt{k,n}.p.NP.s0_neg   + volt{k,n}.soc * ...
                (volt{k,n}.p.NP.s100_neg - volt{k,n}.p.NP.s0_neg);
        z{k,n} = volt{k,n}.p.NP.s0_pos   + volt{k,n}.soc * ...
                (volt{k,n}.p.NP.s100_pos - volt{k,n}.p.NP.s0_pos);
        % the equilibrium potential of the full cell and its components [V]
        volt{k,n}.Eq_neg  = - volt{k,n}.p.U_neg(w{k,n});
        volt{k,n}.Eq_pos  =   volt{k,n}.p.U_pos(z{k,n});
        volt{k,n}.Eq_sep  = zeros(1,length(volt{k,n}.soc));
        volt{k,n}.Eq_cell =   volt{k,n}.p.U_pos(z{k,n})...
                            - volt{k,n}.p.U_neg(w{k,n});
        % relative difference of the full cell and its components [-]
        delta{k,n}.v_neg  = (volt{k,n}.v_neg  - volt{k,n}.Eq_neg);
        delta{k,n}.v_pos  = (volt{k,n}.v_pos  - volt{k,n}.Eq_pos);
        delta{k,n}.v_sep  = (volt{k,n}.v_sep  - volt{k,n}.Eq_sep);
        delta{k,n}.v_cell = (volt{k,n}.v_cell - volt{k,n}.Eq_cell);
    end
end

%% plot v and Δv
make_plot_v_soc(volt,delta,Crates,colors,'neg');
make_plot_v_soc(volt,delta,Crates,colors,'pos');
make_plot_v_soc(volt,delta,Crates,colors,'sep');
make_plot_v_soc(volt,delta,Crates,colors,'cell');

%% output origin
Origin.v{1,1} = [volt{1,1}.soc'...
                 volt{1,1}.Eq_neg'...
                 volt{1,1}.Eq_pos'...
                 volt{1,1}.Eq_sep'...
                 volt{1,1}.Eq_cell'];

for k = 1:size(Crates,1)  % C
    Origin.v{k+1,1} = [volt{k,1}.soc'...
                       volt{k,1}.v_neg'...
                       volt{k,1}.v_pos'...
                       volt{k,1}.v_sep'...
                       volt{k,1}.v_cell'];
    
    Origin.delta_v{k+1,1} = [volt{k,1}.soc'...
                             delta{k,1}.v_neg'...
                             delta{k,1}.v_pos'...
                             delta{k,1}.v_sep'...
                             delta{k,1}.v_cell'];
end

end

function make_plot_v_soc(volt,delta,Crates,colors,local)

switch local
    case 'neg'
        Ueq_soc = volt{1,1}.Eq_neg;
    case 'pos'
        Ueq_soc = volt{1,1}.Eq_pos;
    case 'sep'
        Ueq_soc = volt{1,1}.Eq_sep;
    case 'cell'
        Ueq_soc = volt{1,1}.Eq_cell;
end
for k = 1:size(Crates,1)  % C
    switch local
        case 'neg'
                  v{k,1} =  volt{k,1}.v_neg;
            delta_v{k,1} = delta{k,1}.v_neg;
        case 'pos'
                  v{k,1} =  volt{k,1}.v_pos;
            delta_v{k,1} = delta{k,1}.v_pos;
        case 'sep'
                  v{k,1} =  volt{k,1}.v_sep;
            delta_v{k,1} = delta{k,1}.v_sep;
        case 'cell'
                  v{k,1} =  volt{k,1}.v_cell;
            delta_v{k,1} = delta{k,1}.v_cell;
    end
end

%% v
switch local
    case 'neg'
figure(1);  subplot(2,3,1);  Label = 'Fig .4(a) {\itv}_1';
    case 'pos'
figure(1);  subplot(2,3,4);  Label = 'Fig .4(c) {\itv}_2';
    case 'sep'
figure(2);  subplot(2,3,1);  Label = 'Fig .4(e) {\itv}_3';
    case 'cell'
figure(2);  subplot(2,3,4);  Label = 'Fig .4(g) {\itv}_4';
end
plot(volt{1,1}.soc,Ueq_soc,colors{6}{1},'LineWidth',2)
hold on
for k = 1:size(Crates,1)  % C
    plot(volt{k,1}.soc,v{k,1},colors{k}{1},'LineWidth',2)
end
% Add a legend.
h = legend('OCV',[num2str(Crates(1)),' C'],...
                 [num2str(Crates(2)),' C'],...
                 [num2str(Crates(3)),' C'],'Location','southwest');
% Set the properties of the coordinate axis.
xlabel('SOC','fontsize',24,'fontname','Times')
switch local
    case 'neg'
ylabel('{\itv}_1 (V)','fontsize',24,'fontname','Times')
axis([-0.1 1.1 -1.5 0.5])
    case 'pos'
ylabel('{\itv}_2 (V)','fontsize',24,'fontname','Times')
axis([-0.1 1.1 3.5 4.5])
    case 'sep'
ylabel('{\itv}_3 (V)','fontsize',24,'fontname','Times')
axis([-0.1 1.1 -0.05 0.05])
    case 'cell'
ylabel('{\itv}_4 (V)','fontsize',24,'fontname','Times')
axis([-0.1 1.1 2.5 4.5])
end
% Flip the X axis up and down.
set(gca, 'XDir','reverse')
% Set grid.
set(gca,'FontSize',18,'xgrid','on')
set(gca,'FontSize',18,'ygrid','on')
title(Label,'fontsize',18,'fontweight','bold')

%% Δv
switch local
    case 'neg'
figure(1);  subplot(2,3,2);  Label = 'Fig .4(b) Δ{\itv}_1';
    case 'pos'
figure(1);  subplot(2,3,5);  Label = 'Fig .4(d) Δ{\itv}_2';
    case 'sep'
figure(2);  subplot(2,3,2);  Label = 'Fig .4(f) Δ{\itv}_3';
    case 'cell'
figure(2);  subplot(2,3,5);  Label = 'Fig .4(h) Δ{\itv}_4';
end
for k = 1:size(Crates,1)  % C
    plot(volt{k,1}.soc,delta_v{k,1},colors{k}{1},'LineWidth',2)
    hold on
end
% Add a legend.
h = legend([num2str(Crates(1)),' C'],...
           [num2str(Crates(2)),' C'],...
           [num2str(Crates(3)),' C'],'Location','northwest');
% Set the properties of the coordinate axis.
xlabel('SOC','fontsize',24,'fontname','Times')
switch local
    case 'neg'
ylabel('Δ{\itv}_1 (V)','fontsize',24,'fontname','Times')
axis([-0.1 1.1 -1 0.1])
    case 'pos'
ylabel('Δ{\itv}_2 (V)','fontsize',24,'fontname','Times')
axis([-0.1 1.1 -0.2 0.02])
    case 'sep'
ylabel('Δ{\itv}_3 (V)','fontsize',24,'fontname','Times')
axis([-0.1 1.1 -0.02 0.002])
    case 'cell'
ylabel('Δ{\itv}_4 (V)','fontsize',24,'fontname','Times')
axis([-0.1 1.1 -1 0.1])
end
% Flip the axis up and down.
set(gca,'XDir','reverse')
set(gca,'YDir','reverse')
% Set grid.
set(gca,'FontSize',18,'xgrid','on')
set(gca,'FontSize',18,'ygrid','on')
title(Label,'fontsize',18,'fontweight','bold')

end
