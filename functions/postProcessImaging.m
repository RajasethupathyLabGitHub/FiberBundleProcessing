function [rois,Cn,C,dF_F,dF_F0,T] = postProcessImaging(CNMFE_file,Fs)
% CNMFE_file = path to final output of CNMF-E
% Fs = imaging sample rate

CNMFE_load = load(CNMFE_file);
neuron = CNMFE_load.neuron;

roitemp = reshape(full(neuron.A),size(neuron.Cn,1),size(neuron.Cn,2),size(neuron.A,2));
rois = permute(roitemp,[3,1,2]);

Cn = neuron.Cn;
C = neuron.C_raw;

dF_F = DF_F_runningW(C,20*Fs);

 t_threshold = 3; % standard deviations above noise
 mindur = 1*Fs; % minimum duration for a transient

dF_F0 = zeros(size(dF_F));
T = zeros(size(dF_F));
for i_c = 1:size(dF_F,1)
    [Ti,~,T_idxs] = get_transients(dF_F(i_c,:),C(i_c,:),t_threshold,mindur,10);
    T(i_c,:) = Ti;

    T_idxs = cell2mat(T_idxs);
    tmp = zeros(size(dF_F(i_c,:)));
    tmp(T_idxs) = dF_F(i_c,T_idxs);
    dF_F0(i_c,:) = tmp;
end