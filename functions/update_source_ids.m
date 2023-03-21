function neuron_out = update_source_ids(neuron,to_keep)

neuron_out = neuron;
neuron_out.A_prev = neuron.A;
neuron_out.C_prev = neuron.C_raw;

neuron_out.A = neuron.A(:,to_keep);
neuron_out.C = neuron.C(to_keep,:);
neuron_out.C_raw = neuron.C_raw(to_keep,:);
neuron_out.S = neuron.S(to_keep,:);
neuron_out.ids = neuron.ids(to_keep);
neuron_out.tags = neuron.tags(to_keep);
neuron_out.Coor = neuron.Coor(to_keep);
neuron.P.kernel_pars = neuron.P.kernel_pars(to_keep);

