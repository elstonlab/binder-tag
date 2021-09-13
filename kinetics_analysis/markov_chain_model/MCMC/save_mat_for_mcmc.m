% save all data in .mat for MCMC runs
% from markov_distributions

%experimental
all_exp_data.tag_codiff.X = 260:60:1400;
all_exp_data.tag_codiff.Y = tagsrc_codiff_data_h;
all_exp_data.tag_codiff.Y = transpose(tagsrc_codiff_data_h);
all_exp_data.bin_codiff.Y = transpose(binder_codiff_data_h);
all_exp_data.bin_codiff.X = 260:60:1400;
all_exp_data.tag_all.X = 260:60:1400;
all_exp_data.tag_all.Y = transpose(tagsrc_all_data_h);
all_exp_data.cat.Y = exp_cat;
all_exp_data.cat.X = 1:4;
all_exp_data.proc_1.X = (10:20:1490);
all_exp_data.proc_2.X = (10:20:1490);
all_exp_data.proc_3.X = (10:20:1490);
all_exp_data.proc_4.X = (10:20:1490);
all_exp_data.proc_1.Y = h_i;
all_exp_data.proc_2.Y = h_ii;
all_exp_data.proc_3.Y = h_iii;
all_exp_data.proc_4.Y = h_iv;
save('BT_all_exp_data.mat', 'all_exp_data')

%simulated
all_sim_data.tag_codiff.X = 260:60:1400;
all_sim_data.tag_codiff.Y = tag_codiff_dist;
all_sim_data.bin_codiff.Y = binder_codiff_dist;
all_sim_data.bin_codiff.X = 260:60:1400;
all_sim_data.tag_all.X = 260:60:1400;
all_sim_data.tag_all.Y = dist_all_tag;
all_sim_data.cat.Y = sim_cat;
all_sim_data.cat.X = 1:4;
all_sim_data.proc_1.X = (10:20:1490);
all_sim_data.proc_2.X = (10:20:1490);
all_sim_data.proc_3.X = (10:20:1490);
all_sim_data.proc_4.X = (10:20:1490);
all_sim_data.proc_1.Y = sim_process_1_dist;
all_sim_data.proc_2.Y = sim_process_2_dist;
all_sim_data.proc_3.Y = sim_process_3_dist;
all_sim_data.proc_4.Y = sim_process_4_dist;
save('BT_all_sim_data.mat', 'all_sim_data')