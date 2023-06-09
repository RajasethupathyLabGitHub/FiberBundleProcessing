function run_CNMFE_preprocess(subFolderPath,fileName,opts)
Fs = opts.Fs;
tsub = opts.tsub;
ssub = opts.ssub;
do_residual = opts.do_residual;

with_dendrites = opts.with_dendrites;
    %%
    neuron = Sources2D();
    nam = get_fullname(fullfile(subFolderPath,[fileName,'.mat']));          % this demo data is very small, here we just use it as an example
    nam = neuron.select_data(nam);  %if nam is [], then select data interactively
    
    m = matfile(nam);
    % older versions of mot correct code didnt save a "Ysiz" variable to
    % the .mat file.
    if ~any(strcmp('Ysiz',who(m)))
        m.Properties.Writable = true;
        m.Yr = [];
        m.Ysiz = size(m,'Y');
        m.Properties.Writable = false;
    end
    
    %% parameters
    % -------------------------    COMPUTATION    -------------------------  %
    pars_envs = struct('memory_size_to_use', 64, ...   % GB, memory space you allow to use in MATLAB
        'memory_size_per_patch',4, ...   % GB, space for loading data within one patch
        'patch_dims', [25,25]);  %GB, patch size
%     pars_envs = struct('memory_size_to_use', 4, ...   % GB, memory space you allow to use in MATLAB
%         'memory_size_per_patch', 0.1, ...   % GB, space for loading data within one patch
%         'patch_dims', [25, 25]);  %GB, patch size
    % -------------------------      SPATIAL      -------------------------  %
    gSig = opts.gSig;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
    gSiz = opts.gSiz;          % pixel, neuron diameter
    %ssub = 2;           % spatial downsampling factor
    %with_dendrites = false;   % with dendrites or not
    if with_dendrites
        % determine the search locations by dilating the current neuron shapes
        updateA_search_method = 'dilate';  %#ok<UNRCH>
        updateA_bSiz = 5;
        updateA_dist = neuron.options.dist;
    else
        % determine the search locations by selecting a round area
        updateA_search_method = 'ellipse'; %#ok<UNRCH>
        updateA_dist = 5;
        updateA_bSiz = neuron.options.dist;
    end
    spatial_constraints = struct('connected', true, 'circular', true);  % you can include following constraints: 'circular'
    spatial_algorithm = 'hals_thresh';

    % -------------------------      TEMPORAL     -------------------------  %
    %Fs = 25;%34;             % frame rate
    %tsub = 2;%5;           % temporal downsampling factor
    deconv_flag = true;     % run deconvolution or not 
    deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
        'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
        'smin', -5, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
        'optimize_pars', true, ...  % optimize AR coefficients
        'optimize_b', true, ...% optimize the baseline);
        'max_tau', 200);    % maximum decay time (unit: frame);

    nk = 3; %2             % detrending the slow fluctuation. usually 1 is fine (no detrending)
    % when changed, try some integers smaller than total_frame/(Fs*30)
    detrend_method = 'spline';  % compute the local minimum as an estimation of trend.

    % -------------------------     BACKGROUND    -------------------------  %
    bg_model = 'ring';  % model of the background {'ring', 'svd'(default), 'nmf'}
    nb = opts.nb;             % number of background sources for each patch (only be used in SVD and NMF model)
    ring_radius = opts.ring_radius;  % when the ring model used, it is the radius of the ring used in the background model.
    %otherwise, it's just the width of the overlapping area
    num_neighbors = []; % number of neighbors for each neuron
    bg_ssub = 2;        % downsample background for a faster speed 

    % -------------------------      MERGING      -------------------------  %
    show_merge = false;  % if true, manually verify the merging step
    merge_thr = opts.merge_thr;  % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
    method_dist = 'mean';   % method for computing neuron distances {'mean', 'max'}
    dmin = opts.dmin;       % minimum distances between two neurons. it is used together with merge_thr
    dmin_only = opts.dmin_only;  % merge neurons if their distances are smaller than dmin_only.
    merge_thr_spatial = [0.6, 0.3, -inf];  % merge components with highly correlated spatial shapes (corr=0.8) and small temporal correlations (corr=0.1)

    % -------------------------  INITIALIZATION   -------------------------  %
    K = [];             % maximum number of neurons per patch. when K=[], take as many as possible.
    min_corr = opts.min_corr;%0.65;     % minimum local correlation for a seeding pixel
    min_pnr = opts.min_pnr;%5.5;       % minimum peak-to-noise ratio for a seeding pixel
    min_pixel = (gSiz*0.75)^2;      % minimum number of nonzero pixels for each neuron
    bd = 0;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
    frame_range = [];   % when [], uses all frames
    save_initialization = false;    % save the initialization procedure as a video.
    use_parallel = true;    % use parallel computation for parallel computing
    show_init = true;   % show initialization results
    choose_params = false; % manually choose parameters
    center_psf = true;  % set the value as true when the background fluctuation is large (usually 1p data)
    % set the value as false when the background fluctuation is small (2p)

    % -------------------------  Residual   -------------------------  %
    min_corr_res = opts.min_corr_res;
    min_pnr_res = opts.min_pnr_res;
    seed_method_res = 'auto';  % method for initializing neurons from the residual
    update_sn = true;

    % ----------------------  WITH MANUAL INTERVENTION  --------------------  %
    with_manual_intervention = true;

    % -------------------------  FINAL RESULTS   -------------------------  %
    save_demixed = true;    % save the demixed file or not
    kt = 3;                 % frame intervals

    % -------------------------    UPDATE ALL    -------------------------  %
    neuron.updateParams('gSig', gSig, ...       % -------- spatial --------
        'gSiz', gSiz, ...
        'ring_radius', ring_radius, ...
        'ssub', ssub, ...
        'search_method', updateA_search_method, ...
        'bSiz', updateA_bSiz, ...
        'dist', updateA_bSiz, ...
        'spatial_constraints', spatial_constraints, ...
        'spatial_algorithm', spatial_algorithm, ...
        'tsub', tsub, ...                       % -------- temporal --------
        'deconv_flag', deconv_flag, ...
        'deconv_options', deconv_options, ...
        'nk', nk, ...
        'detrend_method', detrend_method, ...
        'background_model', bg_model, ...       % -------- background --------
        'nb', nb, ...
        'ring_radius', ring_radius, ...
        'num_neighbors', num_neighbors, ...
        'bg_ssub', bg_ssub, ...
        'merge_thr', merge_thr, ...             % -------- merging ---------
        'dmin', dmin, ...
        'method_dist', method_dist, ...
        'min_corr', min_corr, ...               % ----- initialization -----
        'min_pnr', min_pnr, ...
        'min_pixel', min_pixel, ...
        'bd', bd, ...
        'center_psf', center_psf);
    neuron.Fs = Fs;

    %% distribute data and be ready to run source extraction
    neuron.getReady(pars_envs);

    %% initialize neurons from the video data within a selected temporal range
    if choose_params
        % change parameters for optimized initialization
        [gSig, gSiz, ring_radius, min_corr, min_pnr] = neuron.set_parameters();
    end

    [center, Cn, PNR] = neuron.initComponents_parallel(K, frame_range, save_initialization, use_parallel);
    neuron.compactSpatial();
    if show_init
        figure();
        ax_init= axes();
        imagesc(Cn, [0, 1]); colormap gray;
        hold on;
        plot(center(:, 2), center(:, 1), '.r', 'markersize', 10);
    end
    n_init = numel(neuron.ids);

    %% estimate the background components
    neuron.update_background_parallel(use_parallel);
    neuron_init = neuron.copy();

    %%  merge neurons and update spatial/temporal components
    neuron.merge_neurons_dist_corr(show_merge);
    neuron.merge_high_corr(show_merge, merge_thr_spatial);
    n_init_trim = numel(neuron.ids);

    %% pick neurons from the residual
    if do_residual 
        [center_res, Cn_res, PNR_res] =neuron.initComponents_residual_parallel([], save_initialization, use_parallel, min_corr_res, min_pnr_res, seed_method_res);
        if show_init
            axes(ax_init);
            plot(center_res(:, 2), center_res(:, 1), '.g', 'markersize', 10);
        end
        neuron_init_res = neuron.copy();
    end
    n_res = numel(neuron.ids - n_init);

    %% udpate spatial&temporal components, delete false positives and merge neurons
    if update_sn
        neuron.update_spatial_parallel(use_parallel, true);
        udpate_sn = false;
    else
        neuron.update_spatial_parallel(use_parallel);
    end
    % merge neurons based on correlations 
    neuron.merge_high_corr(false, merge_thr_spatial);

    for m=1:2
        % update temporal
        neuron.update_temporal_parallel(use_parallel);
        % delete bad neurons
        neuron.remove_false_positives();
        % merge neurons based on temporal correlation + distances 
        neuron.merge_neurons_dist_corr(false);
    end

    
    Coor = neuron.show_contours(0.6);

    sf = fileparts(nam);
        save(fullfile(sf,[fileName,'results']),'-v7.3');

end
