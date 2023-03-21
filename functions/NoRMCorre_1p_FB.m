% demo file for applying the NoRMCorre motion correction algorithm on 
% 1-photon widefield imaging data using low memory (good for long datasets)
% Example file is provided from the miniscope project page
% www.miniscope.org

%clear;
function NoRMCorre_1p_lowram(name)
gcp;
%% read data and convert to double
%name = 'C:\Users\Andrew T\Desktop\MemoryConsolidation\singlecell\data\rawdata\b790rf\R15\ACC.mat';

[filepath,file_name,ext] = fileparts(name);
tmp_folder = fullfile(filepath,[file_name,'_tmp']);
if ~exist(tmp_folder,'dir')
    mkdir(tmp_folder);
end
movefile(name,tmp_folder);
m = matfile(fullfile(tmp_folder,file_name));

addpath(genpath('C:\Users\Andrew T\Documents\MATLAB\NoRMCorre-master'));

[d1,d2,~] = size(m.Y);
nFrames = size(m,'Y',3);
%% perform some sort of deblurring/high pass filtering
% The function does not load the whole file in memory. Instead it loads 
% chunks of the file and then saves the high pass filtered version in a 
% h5 file.

gSig = 5; 
gSiz = 3*gSig; 
psf = fspecial('gaussian', round(2*gSiz), gSig);
 ind_nonzero = (psf(:)>=max(psf(:,1)));
 psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;   % only use pixels within the center disk


h5_name = fullfile(filepath,[file_name,'_filtered_data.h5']);
chunksize = 2500;    % read 2500 frames at a time
cnt = 1;
if(isfile(h5_name))
    disp('h5 already exists, deleting');
    delete(h5_name);
end
w = size(m.Y,2);
h = size(m.Y,1);
[cols rows] = meshgrid(1:h, 1:w);
in_fov=(rows - h/2).^2 + (cols - w/2).^2 <= (w/2).^2;
                 [I,J] = ind2sub([h,w],find(in_fov));
                 [nI,nJ] = ind2sub([h,w],find(~in_fov));

f = @(x) imfilter(x,psf,'symmetric');
H = fspecial('Gaussian',20,10);
f2 = @(x) imfilter(x,H,'symmetric');
tic
while (1)  % read filter and save file in chunks
    crange = (cnt-1)*chunksize +1 : cnt*chunksize;
    if crange(end)>nFrames
        if crange(1) > size(m.Y,3)
            break
        else
            crange = crange(1):nFrames;
            Yf = m.Y(:,:,crange);% 
         Y  = Yf + max(Yf(:)).*uint16(~in_fov);
             for i = 1:size(Yf,3)
                 Y(:,:,i) = medfilt2(Y(:,:,i),[5,5]);
                 bgY = roifilt2(Y(:,:,i),in_fov,f2);
                 Y(:,:,i) = Y(:,:,i) - bgY;
             end

                      for i = 1:size(Yf,3)
             Y(:,:,i) = imfilter(Y(:,:,i),psf,'symmetric');
                      end
         Y = single(Y);
          Y = Y-min(Y(:)); Y = Y/max(Y(:));
          
            saveash5(Y,h5_name);
            break
        end
    else
         Yf = m.Y(:,:,crange);
         Y  = Yf + max(Yf(:)).*uint16(~in_fov);

             for i = 1:size(Yf,3)
                 Y(:,:,i) = medfilt2(Y(:,:,i),[5,5]);
                 bgY = roifilt2(Y(:,:,i),in_fov,f2);
                 Y(:,:,i) = Y(:,:,i) - bgY;
             end

                 for i = 1:size(Yf,3)
             Y(:,:,i) = imfilter(Y(:,:,i),psf,'symmetric');
                 end
         Y = single(Y);
         Y = Y-min(Y(:)); Y = Y/max(Y(:));
        saveash5(Y,h5_name);
    end
    cnt = cnt+1;
    disp(cnt)
end
disp(' preprocess time: ')
disp(toc)
%% first try out rigid motion correction
    % exclude boundaries due to high pass filtering effects

%% register using the high pass filtered data and apply shifts to original data

gss = [50];
for i = 1:numel(gss)
    gs = gss(i);
    ols = round([20]);
    for j = 1:numel(ols)
        clear options_r options_nr

        ol = ols(j);
        r_name1 = fullfile(filepath,sprintf('%s_reg1_rig.h5',file_name));
        r_name2 = fullfile(filepath,sprintf('%s.mat',file_name));

options_r = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',100,'init_batch',100,...
    'max_shift',5,'iter',1,'correct_bidir',false,'output_type','h5','h5_filename',r_name1);

        tic; [~,shifts1,template1] = normcorre_batch(h5_name,options_r); toc % register filtered data
%         %     % exclude boundaries due to high pass filtering effects
%         %     
%         % % if you save the file directly in memory make sure you save it with a 
%         % % name that does not exist. Change options_r.tiff_filename 
%         % % or options_r.h5_filename accordingly.
        options_r.output_type = 'memmap';
        options_r.mem_filename = r_name2;
         tic; apply_shifts(m,shifts1,options_r); toc % apply shifts to full dataset
%         
        delete(r_name1);
        m1 = matfile(r_name2);
        m1.Properties.Writable = true;
        m1.Yr = [];
        m1.Ysiz = size(m1,'Y');
    end
end
delete(h5_name);
rmdir(tmp_folder,'s');