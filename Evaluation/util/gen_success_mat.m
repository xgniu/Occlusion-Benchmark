function gen_success_mat(seqs, trackers, eval_type, name_trackers, save_path) %#ok<*INUSL>
    
    results_path = '.\results\'; % file path of tracking results
    num_trackers = length(trackers);
    num_seqs = length(seqs);
    
    overlap_thresholds = 0:0.05:1;
    avg_success = zeros(num_trackers,num_seqs,length(overlap_thresholds));
    
    for idxSeq = 1:num_seqs
        seq = seqs{idxSeq};
        gt = dlmread(['D:\Visual Tracking\Benchmark\OCC\' seq.name '\groundtruth_rect.txt']);
        for idxTrk = 1:num_trackers
            tracker = trackers{idxTrk};
            results = load([results_path seq.name '_' tracker.name '.mat']);% load seq_tracker_results
            res = results.bboxes;
            disp([seq.name ' ' tracker.name]);
                                          
            success_scores = zeros(1,length(overlap_thresholds));
            err_overlap = calc_overlap_err(res, gt);
            for tIdx=1:length(overlap_thresholds)
                success_scores(1,tIdx) = sum( err_overlap > overlap_thresholds(tIdx) );
            end
                
            avg_success(idxTrk, idxSeq, :) = success_scores/(size(gt,1)+eps);
            
        end
    end
    
    dataName=[save_path 'SuccessMat_' num2str(num_trackers) 'algs_' eval_type '.mat'];
    save(dataName,'avg_success','name_trackers');
    
end