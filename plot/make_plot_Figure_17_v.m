
function Origin = make_plot_Figure_17_v(volt,Multiple)

%% equilibrium potential: Eq1/Eq2/Eq3/Eq4
for k = 1:size(volt,1)                 % C
    for n = 1:size(volt,2)             % multiple
    % stoichiometric coefficient for porous electrode [-]
    w{k,n} = volt{k,n}.p.NP.s0_neg + volt{k,n}.soc * ...
            (volt{k,n}.p.NP.s100_neg - volt{k,n}.p.NP.s0_neg);
    z{k,n} = volt{k,n}.p.NP.s0_pos + volt{k,n}.soc * ...
            (volt{k,n}.p.NP.s100_pos - volt{k,n}.p.NP.s0_pos);
    % the equilibrium potential of the full cell and its components [V]
    volt{k,n}.Eq_neg  = - volt{k,n}.p.U_neg(w{k,n});
    volt{k,n}.Eq_pos  = volt{k,n}.p.U_pos(z{k,n});
    volt{k,n}.Eq_sep  = zeros(1,length(volt{k,n}.soc));
    volt{k,n}.Eq_cell = volt{k,n}.p.U_pos(z{k,n}) - volt{k,n}.p.U_neg(w{k,n});
    end
end

%% plot Figure 17: v1/v2/v3/v4
make_plot_v(volt,Multiple,'v_neg');
make_plot_v(volt,Multiple,'v_pos');
make_plot_v(volt,Multiple,'v_sep');
make_plot_v(volt,Multiple,'v_cell');

%% output origin
Origin{2,5} = [volt{1,5}.soc' volt{1,5}.Eq_neg'...
                              volt{1,5}.Eq_pos'...
                              volt{1,5}.Eq_sep'...
                              volt{1,5}.Eq_cell'];
for n = 1:size(volt,2) % multiple
    Origin{1,n} = [volt{1,n}.soc' volt{1,n}.v_neg'...
                                  volt{1,n}.v_pos'...
                                  volt{1,n}.v_sep'...
                                  volt{1,n}.v_cell'];
end

end

function make_plot_v(volt,Multiple,change)

% color
colors_1 = {'g','k','b','m','r',...
            'g','k','b','m','r'};

colors_2 = {'--g','--k','--b','--m','--r',...
            '--g','--k','--b','--m','--r'};

switch change
    case 'v_neg'
figure(1);     subplot(2,3,1)
for n = [1 3 5] % Multiple
    plot(volt{1,n}.soc,volt{1,n}.v_neg,colors_1{n},'LineWidth',2)
    hold on
end
plot(volt{1,end}.soc,volt{1,end}.Eq_neg,colors_2{2},'LineWidth',2)
    case 'v_pos'
figure(1);     subplot(2,3,2)
for n = [1 3 5] % Multiple
    plot(volt{1,n}.soc,volt{1,n}.v_pos,colors_1{n},'LineWidth',2)
    hold on
end
plot(volt{1,end}.soc,volt{1,end}.Eq_pos,colors_2{2},'LineWidth',2)
    case 'v_sep'
figure(1);     subplot(2,3,4)
for n = [1 3 5] % Multiple
    plot(volt{1,n}.soc,volt{1,n}.v_sep,colors_1{n},'LineWidth',2)
    hold on
end
plot(volt{1,end}.soc,volt{1,end}.Eq_sep,colors_2{2},'LineWidth',2)
    case 'v_cell'
figure(1);     subplot(2,3,5)
for n = [1 3 5] % Multiple
    plot(volt{1,n}.soc,volt{1,n}.v_cell,colors_1{n},'LineWidth',2)
    hold on
end
plot(volt{1,end}.soc,volt{1,end}.Eq_cell,colors_2{2},'LineWidth',2)
end
% Add a legend.
h = legend(['{\itD}_e: ',num2str(Multiple(1))],...
           ['{\itD}_e: ',num2str(Multiple(3))],...
           ['{\itD}_e: ',num2str(Multiple(5))],'OCV','Location','southwest','FontSize',12);
% Set the properties of the coordinate axis.
xlabel('SOC','fontsize',24,'fontname','Times')
switch change
    case 'v_neg'
ylabel('{\itv}_1 (V)','fontsize',24,'fontname','Times')
axis([-0.1 1.1 -1 0])
title('Fig. 17(a) {\itv}_1','fontsize',18,'fontweight','bold')
    case 'v_pos'
ylabel('{\itv}_2 (V)','fontsize',24,'fontname','Times')
axis([-0.1 1.1 3.5 4.5])
title('Fig. 17(b) {\itv}_2','fontsize',18,'fontweight','bold')
    case 'v_sep'
ylabel('{\itv}_3 (V)','fontsize',24,'fontname','Times')
axis([-0.1 1.1 -0.05 0.05])
title('Fig. 17(c) {\itv}_3','fontsize',18,'fontweight','bold')
    case 'v_cell'
ylabel('{\itv}_4 (V)','fontsize',24,'fontname','Times')
axis([-0.1 1.1 2.5 4.5])
title('Fig. 17(d) {\itv}_4','fontsize',18,'fontweight','bold')
end
% Flip the X axis up and down.
set(gca, 'XDir','reverse')
% Set grid.
set(gca,'FontSize',18)
grid on; % show grid

end
