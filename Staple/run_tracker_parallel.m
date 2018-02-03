function run_tracker_parallel()
    base_path ='D:\Visual Tracking\Benchmark\OTB13\';
    
    dirs = dir(base_path);
    videos = {dirs.name};
    videos(strcmp('.', videos) | strcmp('..', videos) | ...
        strcmp('anno', videos) | ~[dirs.isdir]) = [];
    
    parfor k = 1:numel(videos),
        disp(videos{k});
        if exist(['./OTB13_results/' videos{k} '_bboxes.mat'],'file')~=2,
            run_tracker_seq(videos{k});            
        end
    end
        
end
