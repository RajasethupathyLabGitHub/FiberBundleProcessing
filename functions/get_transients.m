function [T,T_C,transient_idxs] = get_transients(dF,C,thresh,mindur,w)
% dF = dF/F
% C = C matrix
% thresh = # of standard deviations above mean
% mindur = minimum duration of a transient (could be 0.5 or 1 x Fs)
% w  = running mean width

% T = transient locations with height = height of dF/F
% T_C = transient locations with height = height of C

% finds noise cutoff (use thresh = 3)
tol = 0.02;
cutoff = std(dF);
sig = std(dF(dF<cutoff));
while abs(cutoff - thresh*sig)> tol
    if cutoff < thresh*sig
                    cutoff=1.1*cutoff;
                    sig = std(dF(dF<cutoff));
    elseif cutoff > thresh*sig
                    cutoff=0.9*cutoff;
                    sig = std(dF(dF<cutoff));
    end
end
dF_mean = movmean(dF,w);
baseline = mean(dF(dF<cutoff));
dF_cutoff = dF_mean>thresh*sig +baseline;
%T_cutoff = thresh*sig +baseline;

% finds frames where df/f is above the cutoff for at least "mindur" amount of time
ddFi = find([0,diff(dF_cutoff)]>0);
transient_idxs = {};
tcount = 1;
for it = 1:numel(ddFi)
    j = 0;
    while dF_cutoff(ddFi(it)+j)
        j=j+1;
        if ddFi(it) + j > numel(dF_mean)
           break; 
        end
    end
    if j >= mindur
        transient_idxs{tcount} = ddFi(it):ddFi(it)+j-1;
        tcount = tcount+1;
    end
end
% 


T = zeros(size(dF));
T_C = zeros(size(dF));
% places transients in appropriate place
if(numel(transient_idxs)>0)
    for it = 1:numel(transient_idxs)
        idxs = transient_idxs{it};
        [pks,locs] = findpeaks(dF_mean(idxs),'MinPeakProminence',thresh*sig/2,...
            'MinPeakHeight',thresh*sig+baseline,'MinPeakDistance',mindur/2,'MinPeakWidth',mindur/2);
        if numel(pks) > 1
           % If there are multiple peaks here (deals with multiple
           % transient peaks in the same window)
           for i_p = 1:numel(pks)
              if i_p == 1
                  start_loc = idxs(1);
                  T(start_loc) = dF(idxs(locs(1)));
                  T_C(start_loc) = C(idxs(locs(1)));
              else
                  curr_loc = locs(i_p);
                  search_range = idxs(locs(i_p-1)):idxs(locs(i_p));
                  [~,troughs] = findpeaks(-dF(search_range),'MinPeakProminence',sig,'MinPeakWidth',mindur/2);
                  if isempty(troughs)
                      continue;
                  end
                  start_loc = search_range(troughs(end));
                  
                  T(start_loc) = dF(idxs(locs(i_p)));
                  T_C(start_loc) = C(idxs(locs(i_p)));
              end
           end
        else
            % start locatioin will be the nearest time point where the df/f
            % reaches thre threshold 
             start_loc = find(1:numel(dF) < idxs(1) & dF<=(thresh*sig/2),1,'last');
             if isempty(start_loc), start_loc = idxs(1); end
             if it > 1
                 if start_loc < transient_idxs{it-1}(end)
                     start_loc = idxs(1);
                 end
             end
             %start_loc = idxs(1);
           T(start_loc) = max(dF(idxs(1):idxs(end)));
           T_C(start_loc) = max(C(idxs(1):idxs(end)));
        end
        last_loc = start_loc;
    end
end


%%
% finds any transients that still reach the threshold but may have been
% missed initially
[pks,locs] = findpeaks(dF_mean,'MinPeakProminence',3*thresh*sig/4,...
    'MinPeakHeight',thresh*sig+baseline,'MinPeakDistance',mindur/2,'MinPeakWidth',mindur);
currpks = find(T>0);
for  i = 1:numel(locs)
    search_range = locs(i) - 6*mindur:locs(i);
    search_range(search_range<1) = [];
    
    if ~isempty(intersect(currpks,search_range))
        % if this method found a transient already labeled, then skip
        continue
    end
    start_loc = find(1:numel(dF) < locs(i) & dF<=(thresh*sig/2),1,'last');
    if isempty(start_loc)
        [~,troughs] = findpeaks(-dF_mean(search_range),'MinPeakProminence',0.1,'MinPeakWidth',mindur/2);
        if isempty(troughs)
            continue;
        end
        start_loc = search_range(troughs(end));
                              
        T(start_loc) = dF(locs(i));
        T_C(start_loc) = C(locs(i));
    else
        T(start_loc) = dF(locs(i));
        T_C(start_loc) = C(locs(i));
    end

    
    
end
todel = find(T<thresh*sig);
T(todel) = 0;
T_C(todel) = 0;

