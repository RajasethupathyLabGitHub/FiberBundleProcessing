function neuron_out = remove_rough_sources(neuron,thresh)

A = neuron.A;

nC = size(A,2);
d = sqrt(size(A,1));
se = strel('disk',2);
to_keep = true(nC,1);
for i = 1:nC 
    Ai = imerode(full(reshape(A(:,i),[d,d])),se);
    if sum(Ai(:)>0) <thresh
        to_keep(i) = false;
    end
end

neuron_out = update_source_ids(neuron,to_keep);

end