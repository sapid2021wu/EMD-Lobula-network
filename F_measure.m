% F-measure evaluation.

function FMeasure = F_measure(imBinary, imGT)
    
    TP = sum(sum(imGT==255&imBinary==255));		% True Positive 
    TN = sum(sum(imGT<=50&imBinary==0));		% True Negative
    FP = sum(sum((imGT<=50)&imBinary==255));	% False Positive
    FN = sum(sum(imGT==255&imBinary==0));		% False Negative
     
    recall = TP / (TP + FN);
    precision = TP / (TP + FP);
    
    FMeasure = 2.0 * (recall * precision) / (recall + precision);        
    
end
