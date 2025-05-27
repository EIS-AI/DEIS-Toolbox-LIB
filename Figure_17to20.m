%% Script description
% Script used to obtain time-frequency domain datas shown in Figures 17-20 in Ref. [1].

% Author: Yuxuan Bai; Qiu-An Huang
% Affliction: Shanghai University, China
% Email address: yuxuan_bai@shu.edu.cn; hqahqahqa@163.com
% Created March 22, 2025

% References:
% [1] Bai et al., A unified framework to decouple charge transfer reactions
% in the Li-ion battery under dynamic conditions,
% Journal of Energy Chemistry, 2025

%% clear past variables and plots
clc; clear; close all

% add path
addpath('Functions')
addpath('Functions\Frequency functions')
addpath('plot')

mode = 1;   % Choice of charge and discharge mode.   Choices: 1: discharge, 2: charge.

if mode == 1
    soc_init = 1;
elseif mode == 2
    soc_init = 0;
end

% change = 'Ds';
% change = 'rs';
% change = 'k';
% change = 'L';
% change = 'epse';
change = 'De';
% change = 'kappa';
% change = 't_plus';

result = 'constant rate';
% result = 'sensitivity';

switch result
    case 'constant rate'
        Crates(:,1) = 4;
        if mode == 1
            soc_thr = 0.98:-0.1:0.08;
        elseif mode == 2
            soc_thr = 0.08: 0.1:0.98;
        end
    case 'sensitivity'
        Crates(:,1) = 0.25:0.25:4.5;
        if mode == 1
            soc_thr = 0.98:-0.05:0.08;
        elseif mode == 2
            soc_thr = 0.08: 0.05:0.98;
        end
end

f(:,1) = logspace(-3,4,100);

%% settings for parametric studies
switch result
    case 'constant rate'
        switch change
            case {'Ds','k','De'}
                Multiple = [0.5 0.7 1 1.3 1.8];  % Ds, k, De
            case  {'rs','L','kappa'}
                Multiple = 0.5:0.5:1.5;  % rs, L, κ
            case {'epse','t_plus'}
                Multiple = 0.8:0.2:1.2;  % εe, t0+
        end
    case 'sensitivity'
        Multiple = logspace(-1,0,15); Multiple = Multiple(6:10) / Multiple(8);
end

%% Cyclic calculations at different rates
p = parameters_init;            % the parameters used in simulation
p.gridsize = [13 10 13 41 41];  % compact specification of grid parameters

switch change
    case 'Ds'
        Ds_neg = p.Ds_neg * Multiple;  % settings for parametric studies
        Ds_pos = p.Ds_pos * Multiple;  % settings for parametric studies
    case 'rs'
        rs_neg = p.rs_neg * Multiple;  % settings for parametric studies
        rs_pos = p.rs_pos * Multiple;  % settings for parametric studies
    case 'k'
        k_neg = p.k_neg * Multiple;    % settings for parametric studies
        k_pos = p.k_pos * Multiple;    % settings for parametric studies
    case 'L'
        L_neg = p.L_neg * Multiple;    % settings for parametric studies
        L_pos = p.L_pos * Multiple;    % settings for parametric studies
        L_sep = p.L_sep * Multiple;    % settings for parametric studies
    case 'epse'
        epse_neg = p.epse_neg * Multiple;   % settings for parametric studies
        epse_pos = p.epse_pos * Multiple;   % settings for parametric studies
        epse_sep = p.epse_sep * Multiple;   % settings for parametric studies
    case 'De'
        De = @(c,T) p.De(c,T);          % settings for parametric studies
    case 'kappa'
        kappa = @(c,T) p.kappa(c,T);    % settings for parametric studies
    case 't_plus'
        t_plus = p.t_plus * Multiple;   % settings for parametric studies
end

for n = 1:length(Multiple)
    switch change
        case 'Ds'
            p.Ds_neg = Ds_neg(n);    % settings for parametric studies
            p.Ds_pos = Ds_pos(n);    % settings for parametric studies
        case 'rs'
            p.rs_neg = rs_neg(n);    % settings for parametric studies
            p.rs_pos = rs_pos(n);    % settings for parametric studies
        case 'k'
            p.k_neg = k_neg(n);      % settings for parametric studies
            p.k_pos = k_pos(n);      % settings for parametric studies
        case 'L'
            p.L_neg = L_neg(n);      % settings for parametric studies
            p.L_pos = L_pos(n);      % settings for parametric studies
            p.L_sep = L_sep(n);      % settings for parametric studies
        case 'epse'
            p.epse_neg = epse_neg(n);       % settings for parametric studies
            p.epse_pos = epse_pos(n);       % settings for parametric studies
            p.epse_sep = epse_sep(n);       % settings for parametric studies
            epss_neg(n) = (1-p.epse_neg)*0.4824/(0.4824+0.0326);  % settings for parametric studies
            p.epsf_neg = 1-p.epse_neg-epss_neg(n);                % settings for parametric studies
            epss_pos(n) = (1-p.epse_pos)*0.59/(0.59+0.025);       % settings for parametric studies
            p.epsf_pos = 1-p.epse_pos-epss_pos(n);                % settings for parametric studies
        case 'De'
            p.De = @(c,T) Multiple(n) * De(c,T);        % settings for parametric studies
        case 'kappa'
            p.kappa = @(c,T) Multiple(n) * kappa(c,T);  % settings for parametric studies
        case 't_plus'
            p.t_plus = t_plus(n);                       % settings for parametric studies
    end

    p.NP = computer_cap(p,mode);
    % calculate the static electrochemical impedance spectroscope
    out.SEIS{1,n} = frequency_input(f,p,[],[],0,soc_thr,[]);
    for k = 1:length(Crates)
        p.dt = 1/Crates(k);
        switch result
            case 'sensitivity'
                switch change
                    case 'epse'
                        if Crates > 4
                            p.dt = 1/Crates(k)/4;
                        end
                end
        end
        % calculate the galvanostatic charge-discharge
        [out.volt{k,n},vari{k,n}] = DFN_time(p,soc_init,Crates(k),mode);
        % calculate the dynamic electrochemical impedance spectroscope
        out.DEIS{k,n} = frequency_input(f,p,out.volt{k,n},vari{k,n},Crates(k),soc_thr,mode);
        out.define{k,n} = [result,': at rate = ',num2str(Crates(k)),' C and para = ',num2str(Multiple(n)),' ',change];
    end
    dispstat(sprintf('At Multiple : %2.2f, finished the simulation \n',Multiple(n)),'keepthis','timestamp');
end

%% save and plot
out.Crates = Crates;
out.soc_thr = soc_thr;
out.Multiple = Multiple;
out.change = change;

switch result
    case 'constant rate'
        switch change
            case 'Ds'
                save('discharge_17to20_C4_1_Ds','out')
            case 'rs'
                save('discharge_17to20_C4_2_rs','out')
            case 'k'
                save('discharge_17to20_C4_3_k','out')
            case 'L'
                save('discharge_17to20_C4_4_L','out')
            case 'epse'
                save('discharge_17to20_C4_5_epse','out')
            case 'De'
                save('discharge_17to20_C4_6_De','out')
            case 'kappa'
                save('discharge_17to20_C4_7_kappa','out')
            case 't_plus'
                save('discharge_17to20_C4_8_t_plus','out')
        end
        Origin.Fig_17 = make_plot_Figure_17_v(out.volt,out.Multiple);
        Origin.Fig_18 = make_plot_Figure_18_state(out.SEIS,out.DEIS,out.Multiple);
        Origin.Fig_19 = make_plot_Figure_19_Z(out.SEIS,out.DEIS,out.Multiple);
        Origin.Fig_20 = make_plot_Figure_20_ZD(out.SEIS,out.DEIS,out.Multiple);
    case 'sensitivity'
        switch change
            case 'Ds'
                save('discharge_21to22_sens_1_Ds','out')
            case 'rs'
                save('discharge_21to22_sens_2_rs','out')
            case 'k'
                save('discharge_21to22_sens_3_k','out')
            case 'L'
                save('discharge_21to22_sens_4_L','out')
            case 'epse'
                save('discharge_21to22_sens_5_epse','out')
            case 'De'
                save('discharge_21to22_sens_6_De','out')
            case 'kappa'
                save('discharge_21to22_sens_7_kappa','out')
            case 't_plus'
                save('discharge_21to22_sens_8_t_plus','out')
        end
        Origin.Fig_21to22 = make_plot_Figure_21to22_sensitivity(out.SEIS,out.DEIS,out.soc_thr,out.Crates,para);
end
