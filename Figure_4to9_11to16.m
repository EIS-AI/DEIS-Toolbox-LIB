%% Script description
% Script used to obtain time-frequency domain datas shown in Figures 4-9 and 11-16 in Ref. [1].

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
    soc_thr = 0.98:-0.1:0.08;
elseif mode == 2
    soc_init = 0;
    soc_thr = 0.08: 0.1:0.98;
end

% Rate = 'low_C';
Rate = 'high_C';
switch Rate
    case 'low_C'
        Crates(:,1) = [0.01 0.1 1];
    case 'high_C'
        Crates(:,1) = [1 2 4];
end

f(:,1) = logspace(-3,4,100);

%% Simulate static electrochemical impedance spectroscope
p = parameters_init;            % the parameters used in simulation
p.gridsize = [13 10 13 41 41];  % compact specification of grid parameters
p.NP = computer_cap(p,mode);

% calculate the static electrochemical impedance spectroscope
out.SEIS{1,1} = frequency_input(f,p,[],[],0,soc_thr,[]);

%% Cyclic calculations at different rates.
for k = 1:length(Crates)
    p.dt = 1/Crates(k);
    % calculate the galvanostatic charge-discharge
    [out.volt{k,1},out.vari{k,1}] = DFN_time(p,soc_init,Crates(k),mode);
    % calculate the dynamic electrochemical impedance spectroscope
    out.DEIS{k,1} = frequency_input(f,p,out.volt{k,1},out.vari{k,1},Crates(k),soc_thr,mode);
end

%% save and plot
out.Crates = Crates;

switch Rate
    case 'low_C'
save('discharge_10_low_C','out')
Origin.Fig_10 = make_plot_Figure_10_low_C(out.SEIS,out.volt,out.DEIS,out.Crates);
    case 'high_C'
save('discharge_4to9_11to16_high_C','out')
Origin.Fig_4 = make_plot_Figure_4_v(out.volt,out.Crates);
% internal state variable
Origin.Fig_5_left  = make_plot_Figure_5to9_vari(out.SEIS,out.DEIS,'cs_bar_neg');
Origin.Fig_5_right = make_plot_Figure_5to9_vari(out.SEIS,out.DEIS,'cs_bar_pos');
Origin.Fig_6       = make_plot_Figure_5to9_vari(out.SEIS,out.DEIS,'ce');
Origin.Fig_7_left  = make_plot_Figure_5to9_vari(out.SEIS,out.DEIS,'phis_neg');
Origin.Fig_7_right = make_plot_Figure_5to9_vari(out.SEIS,out.DEIS,'phis_pos');
Origin.Fig_8       = make_plot_Figure_5to9_vari(out.SEIS,out.DEIS,'phie');
Origin.Fig_9_left  = make_plot_Figure_5to9_vari(out.SEIS,out.DEIS,'j_neg');
Origin.Fig_9_right = make_plot_Figure_5to9_vari(out.SEIS,out.DEIS,'j_pos');
% cell components
Origin.Fig_11 = make_plot_Figure_11_13_15_Z(out.SEIS,out.DEIS,'Z_neg');
Origin.Fig_13 = make_plot_Figure_11_13_15_Z(out.SEIS,out.DEIS,'Z_pos');
Origin.Fig_15 = make_plot_Figure_11_13_15_Z(out.SEIS,out.DEIS,'Z_sep');
% relative difference
Origin.Fig_12 = make_plot_Figure_12_14_16_diffZ(out.SEIS,out.DEIS,out.Crates,'delta_Mag_neg');
Origin.Fig_14 = make_plot_Figure_12_14_16_diffZ(out.SEIS,out.DEIS,out.Crates,'delta_Mag_pos');
Origin.Fig_16 = make_plot_Figure_12_14_16_diffZ(out.SEIS,out.DEIS,out.Crates,'delta_Mag_sep');
end
