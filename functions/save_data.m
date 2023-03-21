clear;
animal = 'b790lf';
load(fullfile('../../../data/processed/',animal,'processed_data'));
sessions = {'preexposure1','T1','T3','T5','R1','R8','R15','R21','RR1','RR2'};
regions = {'ACC','HPC','AM'};
% C = cell(numel(regions),numel(sessions));
% dF_F = cell(numel(regions),numel(sessions));
% T = cell(numel(regions),numel(sessions));
% behavior = cell(1,numel(sessions));
%Fs = 34;
%%
for i_s = 6%1:numel(sessions)
    session = sessions{i_s};
    disp(session);
    if ~isfile( fullfile('../../../data/processed/',animal,session,'imaging.mat')) ||...
        ~isfile(fullfile('../../../data/processed/',animal,session,'behavior.mat'))
        continue;
    end
    load(fullfile('../../../data/processed/',animal,session,'imaging.mat'));
    load(fullfile('../../../data/processed/',animal,session,'behavior.mat'));
    behavior{i_s} = b_tbl;
    for i_r = 1:numel(regions)
        if isempty(neurons{i_r});
            continue;
        end
        C_raw = neurons{i_r}.C_raw;
        noise_est = zeros(size(C_raw,1),1);
        for i_c = 1:numel(noise_est)
            noise_est(i_c) = GetSn(C_raw(i_c,:));
        end
        Creg = scale_C(neurons{i_r}.C_raw,noise_est,8);

        C{i_r,i_s} = Creg;
        dF_Freg = DF_F_runningW(Creg,20*Fs);
        dF_F{i_r,i_s} = dF_Freg;
        %%
        nC = size(dF_Freg,1);
        Treg = zeros(size(Creg));
        
        t_threshold = 3;
        mindur =1*Fs;
        
        for i_c = 110:nC
            dFi = dF_Freg(i_c,:);
            Ci = Creg(i_c,:);
            [Ti,T_Ci] = get_transients2(dFi,Ci,t_threshold,mindur,round(0.25*Fs));
            Treg(i_c,:) = Ti;
%             figure(2); clf;
%             plot(dFi); hold on;
%             plot(Treg,'linewidth',2);
            figure(1);clf; plot(dFi); hold on; plot(Ti)
                        figure(2);clf; plot(Ci); hold on; plot(T_Ci)

            pause(0.1);

            
        end
        T{i_r,i_s} = Treg;
    end
end
%%
save(fullfile('../../../data/processed/',animal,'processed_data'),'animal','sessions','regions','behavior','C','dF_F','T','Fs','t_threshold');