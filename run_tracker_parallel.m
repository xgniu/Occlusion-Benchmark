function run_tracker_parallel()
    base_path ='D:\Visual Tracking\Benchmark\OCC\';
    
    dirs = dir(base_path);
    videos = {dirs.name};
    videos(strcmp('.', videos) | strcmp('..', videos) | ...
        strcmp('anno', videos) | ~[dirs.isdir]) = [];

%     parfor k = 1:numel(videos)
%         disp(videos{k});
%         if exist(['.\Results1\' videos{k} '_KCF.mat'],'file')~=2
%             run_tracker_seq_KCF(videos{k});
%         end
%     end
    
    parfor k = 1:numel(videos)
        disp(videos{k});
        if exist(['.\Results1\' videos{k} '_DSST_OD.mat'],'file')~=2
            run_tracker_seq_DSST_OD(videos{k});
        end
    end
end
