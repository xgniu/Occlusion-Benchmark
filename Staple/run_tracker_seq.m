function run_tracker_seq(seqName)
    
    base_path = 'D:\Visual Tracking\Benchmark\OTB13\';
    if nargin==0,
        seqName = choose_video(base_path);
    end
    
    tracker = Staple(base_path, seqName);
    num_frames = tracker.seq.num_frames;
    bboxes = zeros( num_frames, 4 );  %tracking results

    time = 0;
    for frame = 2:num_frames,
        im = imread(tracker.seq.img_files{frame});
        tic();
        
        track_output = tracker.track(im);
        bboxes( frame,:) = track_output;
        
        tracker.update(im);
        
        time = time + toc();
    end
    
    save(['./OTB13_results/' seqName '_bboxes'],'bboxes');
end

