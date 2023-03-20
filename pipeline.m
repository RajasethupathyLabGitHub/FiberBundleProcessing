
%%% Andrew Toader
%%% The Rockefeller University

%%% Run this after FOV is cropped and saved as a .tif to extract sources.
%%% You also should run any sort of preprocessing on your behavior data and
%%% align it to image triggers prior to running this code
%%%
%%% It has 4 major steps:
%%%     1) Convert .tif to .mat and run motion correction (if needed)
%%%     2) Initial source extraction via CNMF-E
%%%     3) improved manual verification of CNMF-E results
%%%     4) Save final data with only necessary variables
%%% These steps are particularly valuable 

%%% Requirements: A decent amount of memory, depending on your video size
%%% (min ~32GB for running provided sample dat              a)
%%% Capability of running multiple parallel threads will also help (current
%%% settings divide videos in 3x3=9 patches, to take advantaged of the 12
%%% cores available.

%%% This pipeline will assume the following file structure
%%% Subject name
%%%     - Session name
%%%         - Brain region name (assumes multiple regions)
%%%             -> In this folder will be:
%%%                 1) [region_name].mat, the raw cropped andmotion-corrected
%%%                     video use for source extracion 
%%%                  (sized: height * width * #of frames, in a .mat file)
%%%                 2) folder of CNMF-E source extraction results. If inputting the
%%%                 above file into CNMF-E, the default way this will be
%%%                 saved is as folder titled:
%%%                 [region_name]_source_extraction (if using the Sources2D
%%%                 object)
%%
%%%                 Then, depending on settings, this folder will typically have the following
%%%                 file structure:
%%%                     -frames_xxxx
%%%                         -LOGS_dd_mmm_yy_xx_xx
%%%                             In this final "LOGS" folder is where
%%%                             intermediate results (Sources2D/neuron object) will be saved
%%%                 NOTE: if you don't want to save your data like this,
%%%                 you will have to make changes to how this code reads in
%%%                 the data

clear all;
close all;
addpath(genpath('./CNMF_E_FB'));
addpath(genpath('CelReg-master'));
addpath(genpath('NoRMCorre-master'));
addpath('functions');

% Adjust paths as needed. In this demo, paths are referenced relative to
% the main demo folder, so this script needs to be run from this folder.
dataPath = './test_data/';
outputPath = './processed_data';
mainFolder = 'm44'; % subject folder name
subFolder = 'day1'; % session folder name
fileName = 'ACC'; %  .tif file (without extension)
do_mc = false; % do motion corection
Fs = 34*3; % imaging sampling rate (sample data was collected at 34 Hz and downsampled by a factor of 3)
start_frame = 1; % first frame of imaging to take (when behavior starts)
subFolderPath = fullfile(dataPath,mainFolder,subFolder);

%%% This is where you would load in your behavior file in whaever format
%%% you use. I use a n_frames x n_variables table data structure for mine
behavload = load(fullfile(dataPath,mainFolder,subFolder,'behavior.mat'));
behavior = behavload.b; 
n_frames = 3000;%size(behavior,1);

stop_frame = n_frames;
%% Loading and motion correction
% This wll convert your .tif file to a .mat file. It's important that
% you've already figured out how to align your behavior file to the imaging
% (using triggers from the camera). So there should be a 1-to-1 alignment
% between the behavior file and images. start_frame will be the the first
% frame where your behavior starts and stop_frame will be the last.
% Here I picked 1001 as the start frame just for demonstration purposes, so
% the behavior and imaging aren't technically aligned in this demo. But this 
% will omit the first 1000 frames just as an example.
% This step also runs an edited version of NoRMCorre for motion correction
% if do_mc is set to true. 

preProcessImaging(subFolderPath,fileName,start_frame,stop_frame,do_mc);

%%
%%%%%%%%%% CNMFE options. these were selected empirically for the current
%%%%%%%%%% dataset. You may have to play around with some of
%%%%%%%%%% these to get the most out of your own datasets
CNMFE_opts = struct;

CNMFE_opts.Fs = 34; % temporal framerate
CNMFE_opts.do_residual = false; % decide whether to initialize neurons from residual. 
                     % Sometimes not necessary and takes a lot of time.
                     % worth making true if not detecting many neurons to
                     % start;
CNMFE_opts.tsub = 5; % temporal downsampling
CNMFE_opts.ssub = 2; % spatial downsampling
CNMFE_opts.bg_ssub = 2; % background downsampling                  
CNMFE_opts.gSig = 5;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
CNMFE_opts.gSiz = 15;      % pixel, neuron diameter
CNMFE_opts.bg_model = 'ring';  % model of the background {'ring', 'svd'(default), 'nmf'}
CNMFE_opts.nb = 1;             % number of background sources for each patch (only be used in SVD and NMF model)
CNMFE_opts.ring_radius = 16;  % when the ring model used, it is the radius of the ring used in the background model.
CNMFE_opts.merge_thr = 0.65;
CNMFE_opts.dmin = 5;       % minimum distances between two neurons. it is used together with merge_thr
CNMFE_opts.dmin_only = 5;  % merge neurons if their distances are smaller than dmin_only.
CNMFE_opts.min_corr =0.75; % minimum local correlation for a seeding pixel
CNMFE_opts.min_pnr = 6; %  minimum peak-to-noise ratio for a seeding pixel
CNMFE_opts.min_corr_res = 0.8;
CNMFE_opts.min_pnr_res = 0.7;
CNMFE_opts.with_dendrites = false;
CNMFE_opts.neuron_pixel_thresh = 5; % total number of pixels needed for a cell after using imerode.
                                    % This gets rid of spatially noisy sources that clearly aren't cells



%% Main source extraction scripts using CNMF-E
% first script will save all variables in a mat file called
% "fileName_results.mat'. This can be removed at the end to free up
% storage space if needed.

% for i = 1:nfiles
run_CNMFE_preprocess(subFolderPath,fileName,CNMFE_opts);
toc
%end

CNMFE_file = run_CNMFE_editsources(subFolderPath,fileName,CNMFE_opts);

%% Postprocess and save data in a new folder
[rois,Cn,C,dF_F,dF_F0,T] = postProcessImaging(CNMFE_file,Fs);
% rois = spatial footprints of neurons
% Cn = spatial correlation image (used for alignment)
% C = raw fluorescence output from CNMFE
% dF_F = baseline-normalized fluorescence
% dF_F0
% = baseline-normalized fluorescence with sub-threshold signal zeroed
% T = transient times 

savepath = fullfile(outputPath,mainFolder);
if ~isfolder(savepath)
    mkdir(savepath)
end
roifolder = fullfile(outputPath,mainFolder,'ROIs_unaligned');
if ~isfolder(roifolder)
    mkdir(roifolder)
end
roialignedfolder = fullfile(outputPath,mainFolder,'ROIs_aligned');
if ~isfolder(roialignedfolder)
    mkdir(roialignedfolder)
end
save(fullfile(savepath,[fileName,'_processed_',subFolder,'.mat']),'rois','Cn','C','dF_F','dF_F0','T');
save(fullfile(roifolder,[fileName,'_rois_',subFolder,'.mat']),'Cn','rois');
