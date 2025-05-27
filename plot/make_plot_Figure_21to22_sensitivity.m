
function Origin = make_plot_Figure_21to22_sensitivity(SEIS,DEIS,soc_thr,Crates,para)

Crate = [0; Crates];

%% static/dynamic EIS
for n = 1:size(SEIS,2)                    % Multiple
    for m = 1:length(SEIS{1,n}.soc_thr)   % soc
        EIS{1,m}.ZDs_pos(:,n) = SEIS{1,n}.Z{m}.Ds_pos;
        EIS{1,m}.ZDe_pos(:,n) = SEIS{1,n}.Z{m}.De_pos;
        for k = 1:size(DEIS,1)            % C
            EIS{k+1,m}.ZDs_pos(:,n) = DEIS{k,n}.Z{m}.Ds_pos;
            EIS{k+1,m}.ZDe_pos(:,n) = DEIS{k,n}.Z{m}.De_pos;
        end
    end
end

%% sensitivity calculate
for k = 1:size(EIS,1)      % C
    for m = 1:size(EIS,2)  % soc
        SD.ZDs_pos(k,m) = sensitivity_calculate(EIS{k,m}.ZDs_pos) * 1e2;
        SD.ZDe_pos(k,m) = sensitivity_calculate(EIS{k,m}.ZDe_pos) * 1e2;
    end
end

%% plot
make_plot_SD(soc_thr,Crate,SD.ZDs_pos,'ZDs_pos',para)
make_plot_SD(soc_thr,Crate,SD.ZDe_pos,'ZDe_pos',para)

%% output origin
Origin.SD_ZDs_pos = SD.ZDs_pos';
Origin.SD_ZDe_pos = SD.ZDe_pos';

Origin.max_min(1,:) = [min(min(SD.ZDs_pos)) max(max(SD.ZDs_pos))];
Origin.max_min(2,:) = [min(min(SD.ZDe_pos)) max(max(SD.ZDe_pos))];

end

function SD = sensitivity_calculate(Z)

% This function calculates the relative sensitivity of the parameter X to the diffusion processes.

% Inputs:
%    - Z is a MATLAB struct containing the diffusion impedance and the interface impedance for Li-ion batteries,
%      such as the solid diffusion impedance Z_Ds(s) or the electrolyte diffusion impedance Z_De(s).

% Outputs:
%    - SD is the relative sensitivity of the solid/electrolyte diffusion impedance to parameter set X.

% Note: Change range of parameter X is set as:
% Multiple = logspace(-1,0,15); Multiple = Multiple(6:10) / Multiple(8);
% i.e., Multiple(:,1) = [0.7197, 0.8483, 1.0000, 1.1788, 1.3895];

N = length(Z(1,:));     % N = 5 in the present study
M = length(Z(:,1));     % M = 43 for diffusion impedance

% the average value of Z_fi at frequency f over i = 1~N times simulating results.
% find the average of each row.
Z_av = mean(Z(1:M,:) , 2);                    % Equation (87a)

% the standard derivation of Z_fi over M sampling frequency points.
Z_M = abs(Z(1:M,:) - Z_av);
% find the sum of each row.
Z_S = (sum(Z_M .^ 2 , 2) / N) .^ 0.5;         % Equation (87b)

% the relative sensitivity of the parameter X to the solid/electrolyte diffusion processes.
SD = (sum(Z_S) / M) / (sum(abs(Z_av)) / M);   % Equation (87c)

end

function make_plot_SD(soc_thr,Crate,SD,change,para)

switch para
    case 'Ds'
        Label = '(a)'; Para = 'to {\itD}_{s2}';
    case 'rs'
        Label = '(b)'; Para = 'to {\itr}_{s2}';
    case 'k'
        Label = '(c)'; Para = 'to {\itk}_2';
    case 'L'
        Label = '(d)'; Para = 'to {\itL}_2';
    case 'epse'
        Label = '(e)'; Para = 'to {\itε}_{e2}';
    case 'De'
        Label = '(f)'; Para = 'to {\itD}_e';
    case 'kappa'
        Label = '(g)'; Para = 'to {\itκ}_{s2}';
    case 't_plus'
        Label = '(h)'; Para = 'to {\itt}^0_+';
end

switch change
    case 'ZDs_pos'
        figure(1);  subplot(2,3,1);  Fig = 'Fig .21';  Pur = ' {\itZ}_{Ds2} ';
    case 'ZDe_pos'
        figure(1);  subplot(2,3,2);  Fig = 'Fig .22';  Pur = ' {\itZ}_{De2} ';
end
% create a data matrix
[X, Y] = meshgrid(soc_thr,Crate);
% draw surface
surf(X, Y, SD);
% add labels and titles
xlabel('SOC');
ylabel('C');
zlabel('SZ (%)');
% xlim([0 1])
xticks([0:0.2:1])
% ylim([0 5])
yticks([0: 1 :5]);
% zlim([0 50])
switch change
    case 'ZDs_pos'
        switch para
            case {'Ds','k','De','kappa','t_plus'}
                zlim([0  20])
            case {'rs','L'}
                zlim([0  40])
            case 'epse'
                zlim([0 160])
        end
    case 'ZDe_pos'
        switch para
            case {'Ds','rs','k','De','kappa'}
                zlim([0 20])
            case {'L','t_plus'}
                zlim([0 40])
            case 'epse'
                zlim([0 80])
        end
end
% flip the axis up and down.
set(gca, 'XDir','reverse')
% set grid.
set(gca,'FontSize',18,'xgrid','on')
set(gca,'FontSize',18,'ygrid','on')
title([Fig,Label,Pur,Para],'fontsize',18,'fontweight','bold')

end
