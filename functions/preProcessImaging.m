function preProcessImaging(data_path,filename,start_frame,stop_frame,do_mc)

%
tif_name = [filename,'.tif'];
file_name = filename;
% Converts .tif to .mat
    m1 = tif2mat(fullfile(data_path,tif_name));
    m = matfile(m1);
    m.Properties.Writable = true;
    m.Y = m.Y(:,:,start_frame:stop_frame);
    m.Ysiz = size(m,'Y');
        m.Properties.Writable = false;

   m = matfile(fullfile(data_path,[filename,'.mat']));
% motion correction
if do_mc
    NoRMCorre_1p_FB(m.Properties.Source);
end

