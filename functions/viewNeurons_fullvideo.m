function [obj_out,ind_del] = viewNeurons_fullvideo(obj, ind, C2, folder_nm,path_to_file)
%% view all components and delete components manually. it shows spatial
%   components in the full-frame and zoomed-in view. It also shows temporal
%   components
%% input:
%   ind: vector, indices of components to be displayed, no bigger than the maximum
%       number of neurons
%   C2:  K*T matrix, another temporal component to be displayed together
%       with the esitmated C. usually it is C without deconvolution.

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016
global t;
global lplot;
global cmax;
[~,~,ext] = fileparts(path_to_file);
disp('Loading full vid');
f_start = 200;
f_end = size(obj.C,2);
switch ext
    case '.tif'
        m = TIFFStack(path_to_file);
        vload = m(:,:,f_start:round(obj.Fs/6):f_end);

    case '.mat'
        m = matfile(path_to_file);
        vload = m.Y(:,:,f_start:round(obj.Fs/6):f_end);

end
disp('done');

if ~exist('ind', 'var') || isempty(ind)
    % display all neurons if ind is not specified.
    ind = 1:size(obj.A, 2);
elseif ind==-1
    ind = size(obj.A,2):-1:1;
end
if ~exist('C2', 'var'); C2=[]; end

if exist('folder_nm', 'var')&&(~isempty(folder_nm))
    % create a folder to save images
    save_img = true;
    cur_cd = cd();
    if ~exist(folder_nm, 'dir'); mkdir(folder_nm);
    else
        fprintf('The folder has been created and old results will be overwritten. \n');
    end
    cd(folder_nm);
else
    save_img = false;
end

% obj.delete(sum(obj.A>0, 1)<max(obj.options.min_pixel, 1));

Amask = (obj.A~=0);
ind_trim = false(size(ind));    % indicator of trimming neurons
ind_del = false(size(ind));     % indicator of deleting neurons
ctr = obj.estCenter();      %neuron's center
gSiz = obj.options.gSiz;        % maximum size of a neuron
if isempty(obj.Coor) || (size(obj.A, 2)~=length(obj.Coor))   % contours of the neuron has not been calculated
    figure;
    obj.Coor = obj.get_contours(0.7);
end
obj.Coor = obj.get_contours(0.7);

Coor = obj.Coor; 
Cn = obj.Cn;   
% time
T = size(obj.C, 2);
t = f_start:T;
if ~isnan(obj.Fs)
    t = t/obj.Fs;
    str_xlabel = 'Time (Sec.)';
else
    str_xlabel = 'Frame';
end

%% keep the log information
if ~save_img
    try
        log_file =  obj.P.log_file;
        flog = fopen(log_file, 'a');
        log_data = matfile(obj.P.log_data, 'Writable', true); %#ok<NASGU>
        manual_intervention.before = obj.obj2struct();
        
        fprintf(flog, '[%s]\b', get_minute());
        fprintf(flog, 'Start manual interventions:\n');
    end
end

%% start viewing neurons
f1= figure('position', [0, 0, 1024, 800]);
m=1;
while and(m>=1, m<=length(ind))
    %% full-frame view
    figure(f1)
    subplot(231); cla;
    obj.image(Cn, [0,1]); hold on; colormap winter;
    axis equal off tight;
    for k=1:m
        % plot contour
        tmp_con = Coor{ind(k)};
        cont_del = (sum(tmp_con<=1, 1)>0);
        tmp_con(:, cont_del) = [];
        if isempty(tmp_con)
            plot(ctr(m, 2), ctr(m, 2));
            disp('empty tmp_con');
        else           
            if and(k<m, ~ind_del(k))
                plot(tmp_con(1, 1:end), tmp_con(2, 1:end), 'color','k', 'linewidth', 1);          
            elseif k==m
                plot(tmp_con(1, 1:end), tmp_con(2, 1:end), 'r', 'linewidth', 2);
            end
        end
        
        title(sprintf('Neuron %d', ind(m)));
    end
    axis equal; axis off;
    if ind_del(m)
        title(sprintf('Neuron %d', ind(m)), 'color', 'r');
    else
        title(sprintf('Neuron %d', ind(m)));
    end
    %% zoomed-in view
    figure(f1)
    subplot(232); cla;
    obj.image(obj.A(:, ind(m)).*Amask(:, ind(m))); %
    %     imagesc(reshape(obj.A(:, ind(m)).*Amask(:,ind(m))), obj.options.d1, obj.options.d2));
    axis equal; axis off;
    colormap(gca, jet);
    x0 = ctr(ind(m), 2);
    y0 = ctr(ind(m), 1);
    if ~isnan(x0)
        xlim(x0+[-gSiz, gSiz]*2);
        ylim(y0+[-gSiz, gSiz]*2);
    end
    
    %% temporal components
    figure(f1)
    subplot(2,3,4:6);cla;
    if ~isempty(C2)
        plot(t, C2(ind(m), f_start:end)*max(obj.A(:, ind(m))), 'linewidth', 2); hold on;
        plot(t, obj.C(ind(m), f_start:end)*max(obj.A(:, ind(m))), 'r');
    else
        
        plot(t, obj.C(ind(m), f_start:end)*max(obj.A(:, ind(m))));
    end
    
    cmax = max( obj.C(ind(m), f_start:end)*max(obj.A(:, ind(m))) );
    xlim([t(1), t(end)]);
    xlabel(str_xlabel);
    
    tmax = find(obj.C(ind(m), f_start:end) == max(obj.C(ind(m), f_start:end)));
    %% save images
    if save_img
        drawnow();
        saveas(gcf, sprintf('neuron_%d.png', ind(m)));
        m = m+1;
    else
        fprintf('Neuron %d, keep(k, default)/delete(d,l)/split(s)/trim(t)\n\t/trim cancel(tc)/delete all(da)/backward(b)/end(e):    ', ind(m));
%         f2 = figure('Position',[200,400,500,500]);
%         figure(f2);
        subplot(2,3,3);cla;
        ts1 = 0.5 / (size(obj.C, 2)/obj.Fs); % half a second
        ts2 = 5 / (size(obj.C, 2)/obj.Fs); % 5 seconds
        SliderHandle = uicontrol('Style', 'slider','SliderStep',[ts1,ts2],...
            'Position',[120,20,810,20],'Callback',@slider_callback);
        global currA;
        global ph;
        global pw;
        global vshow;
        global maskshow;
        ph = size(obj.Cn,1);
        pw = size(obj.Cn,2);
        currA = obj.A(:, ind(m));
        
        val = (tmax/numel(f_start:size(obj.C,2)));
        if val >1, val = 1; end
        set(SliderHandle,'Value',val);
        inds = find(currA> 0);
        [rs,cs] = ind2sub([ph,pw],inds);
        rmed = round(median(rs));
        cmed = round(median(cs));
        vidh = rmed - 45: rmed + 45; 
        vidh(vidh > ph) = []; vidh(vidh<1) = [];
        vidw = cmed - 45: cmed + 45; 
        vidw(vidw > ph) = []; vidw(vidw<1) = [];
        currmask = false(ph,pw);
        currmask(inds) = true;
        currmask = boundarymask(currmask);
        maskshow = currmask(vidh,vidw);
        vshow = double(vload(vidh,vidw,:));
        vshow = vshow-mean(vshow,3);
        slider_callback(SliderHandle)
        pause(0.1);
        
        val = ts1+(tmax/numel(f_start:size(obj.C,2)));
        if val>1, val = 1; end
        set(SliderHandle,'Value',val);
        inds = find(currA> 0);
        [rs,cs] = ind2sub([ph,pw],inds);
        rmed = round(median(rs));
        cmed = round(median(cs));
        vidh = rmed - 45: rmed + 45; 
        vidh(vidh > ph) = []; vidh(vidh<1) = [];
        vidw = cmed - 45: cmed + 45; 
        vidw(vidw > ph) = []; vidw(vidw<1) = [];
        currmask = false(ph,pw);
        currmask(inds) = true;
        currmask = boundarymask(currmask);
        maskshow = currmask(vidh,vidw);
        vshow = double(vload(vidh,vidw,:));
        vshow = vshow-mean(vshow,3);
        slider_callback(SliderHandle)      

         
        temp = input('', 's');
        if temp=='d'
            ind_del(m) = true;
            m = m+1;
        elseif strcmpi(temp,'l')
            ind_del(m) = true;
            m = m+1;
        elseif strcmpi(temp, 'b')
            m = m-1;
        elseif strcmpi(temp, 'da')
            ind_del(m:end) = true;
            break;
        elseif strcmpi(temp, 'k')
            ind_del(m) = false;
            m= m+1;
        elseif strcmpi(temp, 's')
            try
                subplot(232);
                temp = imfreehand();
                tmp_ind = temp.createMask();
                tmpA = obj.A(:, ind(m));
                obj.A(:, end+1) = tmpA.*tmp_ind(:);
                obj.C(end+1, :) = obj.C(ind(m), :);
                obj.A(:, ind(m)) = tmpA.*(1-tmp_ind(:));
                obj.S(end+1, :) = obj.S(ind(m), :);
                obj.C_raw(end+1, :) = obj.C_raw(ind(m), :);
                obj.P.kernel_pars(end+1, :) = obj.P.kernel_pars(ind(m), :);
                k_ids = obj.P.k_ids;
                obj.ids(end+1) = k_ids+1;   % assign an neuron id
                obj.tags(end+1) = obj.tags(ind(m));
                obj.P.k_ids = k_ids+1;
                fprintf(flog, '\tSplit %d --> %d + %d\n', obj.ids(ind(m)),obj.ids(ind(m)), k_ids);
            catch
                fprintf('the neuron was not split\n');
            end
        elseif strcmpi(temp, 't')
             try
                subplot(232);
                temp = imfreehand();
                tmp_ind = temp.createMask();
                Amask(:, ind(m)) = tmp_ind(:);
                Coor_new = get_contour_single(obj,0.9,obj.A(:, ind(m)).*Amask(:, ind(m)));
                Coor{ind(m)} = Coor_new{1};
                ind_trim(m) = true;
            catch
                fprintf('the neuron was not trimmed\n');
            end
        elseif strcmpi(temp, 'tc')
            Amask(:, ind(m)) = (obj.A(:, ind(m)) > 0);
            ind_trim(m) = false;
        elseif strcmpi(temp, 'e')
            break;
        elseif ~isnan(str2double(temp))
            m = m + floor(str2double(temp));
            m = max(m, 1);
            m = min(m, length(ind));
            fprintf('jump to neuron %d / %d\n', m, length(ind));
        else
            m = m+1;
        end
    end
end
if save_img
    cd(cur_cd);
else
    if ~isempty(ind(ind_trim))
        obj.A(:, ind(ind_trim)) = obj.A(:,ind(ind_trim)).*Amask(:, ind(ind_trim));
        try
            fprintf(flog, '\n\tFollowing neurons were trimmed:\n');
            ids_trimmed = ind(ind_trim);
            for m=1:length(ids_trimmed)
                fprintf(flog, '%2d, ', ids_trimmed(m));
            end
            fprintf(flog, '\n');
        end
    end
    
    if ~isempty(ind(ind_del))
        try
            fprintf(flog, '\tDeleting manually selected neurons:\n');
        end
        obj.delete(ind(ind_del));
    end
         obj.Coor = obj.get_contours(0.9);
    
    obj_out = obj;

    return;
end
try
    fprintf(flog, '[%s]\b', get_minute());
    fprintf(flog, 'Finished the manual intervention.\n');
    fprintf(flog, '[%s]\b', get_minute());
    if obj.options.save_intermediate
        manual_intervention.after = obj.obj2struct(); %#ok<STRNU>
        tmp_str = get_date();
        tmp_str=strrep(tmp_str, '-', '_');
        eval(sprintf('log_data.manual_%s = manual_intervention;', tmp_str));
        
        fprintf(flog, '\tThe results were saved as intermediate_results.manual%s\n\n', tmp_str);
    end
    fclose(flog);
    obj_out = obj;

end
close(f1);
obj_out = obj;
end



function slider_callback(ObjH,EventData);
global currA;
global ph;
global pw;
global vshow;
global maskshow;
global t;
global lplot;
global cmax;
val = get(ObjH,'Value');
it = round(val*size(vshow,3));
l = round(val*numel(t));
if(it==0), it = 1; end;
if(l==0),l = 1; end;

try
subplot(2,3,4:6); hold on;
delete(lplot);
lplot = plot([t(l),t(l)],[-15,cmax],'k','linewidth',1);

subplot(2,3,3);
cla;
           axis equal; axis off;

if max(vshow(:)) >= 2*cmax
    crange = [min(vshow(:)),3*cmax];% [prctile(vshow(:),5),prctile(vshow(:),95)];
else
    crange = [min(vshow(:)),3*max(vshow(:))];
end


imagesc(vshow(:,:,it)); colormap(gca,'gray'); hold on;
           axis equal; axis off;

try
% caxis(crange);
catch;
%     crange = [-50,500];
%     caxis(crange);
end
imo = imagesc(maskshow);
set(imo,'Alphadata',maskshow); 
           axis equal; axis off;
drawnow;
hold off;
catch
   a = 0;
   obj_out = obj;

end
end