clear; close all;
addpath('.\util');

%% Load infos of seqs and trackers
seqs = config_occ_seqs();
num_seqs = length(seqs);
name_seqs = cell(num_seqs,1); % name of all seqs
for idxSeq=1:num_seqs
    t = seqs{idxSeq};
    name_seqs{idxSeq}=t.name;
end

trackers = config_trackers();
num_trackers = length(trackers);
name_trackers = cell(num_trackers,1); % name of all trackers
for idxTrk=1:num_trackers
    t = trackers{idxTrk};
    name_trackers{idxTrk}=t.namePaper;
end

success_mat = load(".\perfMat\SuccessMat_8algs_OPE.mat");
success = success_mat.avg_success(:,:,11); % (track_idx, seq_idx, threshold_idx)