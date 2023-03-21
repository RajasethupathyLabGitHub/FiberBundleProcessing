function dF_F = DF_F_runningW(C_in,w)
dF_F = zeros(size(C_in));
nC = size(C_in,1);
T = size(C_in,2);

if mod(w,2)==0
    w = w+1;
end
hw = (w-1)/2;

step = 10;
for iC = 1:nC
    curr_C = C_in(iC,:);
    bsd = std(curr_C(curr_C<median(curr_C))); % noise level
    baseline = zeros(size(curr_C));
    padval = prctile(curr_C,8);
    padC = [padval*ones(1,hw),curr_C,padval*ones(1,hw)];
    parfor it = 1:T
        baseline(it) = prctile(padC(hw + (it-hw):step:(it+hw) + hw),8);
    end
    dF_F(iC,:) = (curr_C-baseline)./bsd;
end

end

