function gen_precision_mat(seqs, trackers, eval_type, name_trackers, save_path) %#ok<*INUSL>
    
    results_path = '.\results\'; % file path of tracking results
    num_trackers = length(trackers);
    num_seqs = length(seqs);

    center_error_thresholds = 0:0.05:1;    
    avg_precision = zeros(num_trackers,num_seqs,length(center_error_thresholds));
    
    for idxSeq=1:num_seqs
        seq = seqs{idxSeq};
        gt = dlmread(['D:\Visual Tracking\Benchmark\OCC\' seq.name '\groundtruth_rect.txt']);
        			
        for idxTrk=1:num_trackers
            tracker = trackers{idxTrk};
            results = load([results_path seq.name '_' tracker.name '.mat']);% load seq_tracker_results
            res = results.bboxes;
            disp([seq.name ' ' tracker.name]);
            
            precision_scores = zeros(1,length(center_error_thresholds));
			cle = calc_center_err(res, gt);                
            for tIdx=1:length(center_error_thresholds)
                precision_scores(1,tIdx) = sum( cle <= center_error_thresholds(tIdx) );
            end
            
            avg_precision(idxTrk, idxSeq, :) = precision_scores/(size(gt,1)+eps);
        end
    end

    dataName=[save_path 'PrecisionMat_' num2str(num_trackers) 'algs_' eval_type '.mat'];
    save(dataName,'avg_precision','name_trackers');
end