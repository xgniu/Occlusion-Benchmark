function plot_draw_save(numTrk,plotDrawStyle,avg_scores,idxSeqSet,rankNum,rankingType,rankIdx,...
        nameTrkAll,thresholdSet,titleName,xLabelName,yLabelName,figName)
    
    for idxTrk=1:numTrk
        tmp=avg_scores(idxTrk, idxSeqSet,:);
        aa=reshape(tmp,[length(idxSeqSet),size(avg_scores,3)]);
        aa=aa(sum(aa,2)>eps,:);
        bb=mean(aa);
        switch rankingType
            case 'AUC'
                perf(idxTrk) = mean(bb);
            case 'threshold'
                perf(idxTrk) = bb(rankIdx);
        end
    end
    [~,indexSort]=sort(perf,'descend');
        
    figure1 = figure;   
    axes1 = axes('Parent',figure1,'FontSize',14);
    i=1;
    for idxTrk=indexSort(1:rankNum)
        tmp=avg_scores(idxTrk,idxSeqSet,:);
        aa=reshape(tmp,[length(idxSeqSet),size(avg_scores,3)]);
        aa=aa(sum(aa,2)>eps,:);
        bb=mean(aa);
        
        switch rankingType
            case 'AUC'
                score = mean(bb);
                tmp=sprintf('%.3f', score);
            case 'threshold'
                score = bb(rankIdx);
                tmp=sprintf('%.3f', score);
        end
        
        tmpName{i} = [nameTrkAll{idxTrk} ' [' tmp ']'];
        plot(thresholdSet,bb,'color',plotDrawStyle{i}.color, 'lineStyle', plotDrawStyle{i}.lineStyle,'lineWidth', 4,'Parent',axes1);
        hold on
        i=i+1;
    end
    
    legend(tmpName,'Interpreter', 'none','fontsize',10);
    title(titleName,'fontsize',16);
    xlabel(xLabelName,'fontsize',16);
    ylabel(yLabelName,'fontsize',16);
    
    hold off
    
    saveas(gcf,figName,'png');
    
end