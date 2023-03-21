function file_path_out = run_CNMFE_editsources(subFolderPath,fileName,opts)
%%% adapted from demo_large_data_1p by:
%%% Andrew Toader
%%% The Rockefeller University

load(fullfile(subFolderPath,[fileName,'results.mat']));

thresh = opts.neuron_pixel_thresh;
neuron = remove_rough_sources(neuron,thresh);
    
Coor = neuron.show_contours(0.6);

%% manual add missing cellsn
if(1)
coords_to_add = [];
title('Pick missing cells, then press any key')
while(1)
    [x,y,button] = ginput(1);
    x = round(x);
    y = round(y);
    if(button>10)
        close(gcf);
       break
    else
    coords_to_add = [coords_to_add;[x,y]];
    plot(x,y,'r*','markersize',10);
    end
end
%%K
    ssub = 1;           % spatial downsampling factor
neuron.updateParams('ssub', ssub)
if numel(coords_to_add)>0
[neuron,center_res, Cn_res, PNR_res] =init_addPoints(neuron,[], save_initialization, use_parallel, min_corr_res, min_pnr_res, seed_method_res,coords_to_add);
end
Coor = neuron.show_contours(0.6);
title('Final initialization');

%% udpate spatial&temporal components, delete false positives and merge neurons
% update spatial
% if update_sn
%     neuron.update_spatial_parallel(use_parallel, true);
%     udpate_sn = false;
% else
%     neuron.update_spatial_parallel(use_parallel);
% end
% merge neurons based on correlations 
neuron.merge_high_corr(show_merge, merge_thr_spatial);

for m=1:2
    % update temporal
    neuron.update_temporal_parallel(use_parallel);
    
    % delete bad neurons
    neuron.remove_false_positives();
    
    % merge neurons based on temporal correlation + distances 
    neuron.merge_neurons_dist_corr(show_merge);
end
end
%%
neuron.options.spatial_algorithm = 'nnls';
if with_manual_intervention
    show_merge = true;
    %dmin_only = 3;
    neuron.merge_close_neighbors(true, dmin_only);
    
    neuron.update_background_parallel(use_parallel);
    neuron.update_temporal_parallel(use_parallel);
    Coor = neuron.show_contours(0.6);
    title('Post-merge');

    neuron.orderROIs('snr');   % order neurons in different ways {'snr', 'decay_time', 'mean', 'circularity'}
    [neuron,~] = viewNeurons_fullvideo(neuron,[], neuron.C_raw,[],nam);
    
    % merge closeby neurons
    neuron.merge_close_neighbors(true, dmin_only);
    
    % delete neurons
    tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
    ids = find(tags>0); 
    if ~isempty(ids)
        [neuron,~] = viewNeurons_fullvideo(neuron,[], neuron.C_raw,[],nam);
    end
end

%%
%% run more iterations
neuron.update_background_parallel(use_parallel);
neuron.update_spatial_parallel(use_parallel);
neuron.update_temporal_parallel(use_parallel);

K = size(neuron.A,2);
tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
neuron.remove_false_positives();
neuron.merge_neurons_dist_corr(show_merge);
neuron.merge_high_corr(show_merge, merge_thr_spatial);

if K~=size(neuron.A,2)
    %neuron.update_spatial_parallel(use_parallel);
    neuron.update_temporal_parallel(use_parallel);
    neuron.remove_false_positives();
end
    neuron.merge_close_neighbors(true, dmin_only);


%% save the workspace for future analysis
neuron.orderROIs('snr');
%%
fprintf('------------- SAVE THE WHOLE WORKSPACE ----------\n\n');
            
            neuron.compress_results();
            file_path_out = fullfile(neuron.P.log_folder,  ['final_',strrep(get_date(), ' ', '_'), '.mat']);
            log_file = neuron.P.log_file; 
            if exist('original_logfolder', 'var')
               neuron.P.log_folder = original_log_folder;  
            end
            if exist('original_logfile', 'var')
                neuron.P.log_file = original_logfile; 
            end
            save(file_path_out, 'neuron', 'save_*', 'show_*', 'use_parallel', 'with_*', '-v7.3');
            try
                fp = fopen(log_file, 'a');
                fprintf(fp, '\n--------%s--------\n[%s]\bSave the current workspace into file \n\t%s\n\n', get_date(), get_minute(), file_path_out);
                fprintf('The current workspace has been saved into file \n\t%s\n\n', file_path_out);
                fp.close();
            end
            
            if exist('zip_file_path', 'var') && ~isempty(zip_file_path)
               [zip_dir, zip_name, ~] = fileparts(get_fullname(zip_file_path)); 
               zip([zip_name, '.zip'], {file_path_out, log_file}, zip_dir); 
               fprintf('The results and the log files were compresed into file \n%s\n', zip_file_path); 
            end



%[neuron,cnmfe_path] = save_whole_workspace(neuron);

%% show neuron contours
Coor = neuron.show_contours(0.2);

%% save neurons shapes
neuron.save_neurons();

