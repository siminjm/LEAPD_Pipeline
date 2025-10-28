function metrics = evaluate_classification(scores, labels)
% EVALUATE_CLASSIFICATION - threshold at 0.5 and compute metrics
labels = labels(:);
scores = scores(:);
pred = double(scores>=0.5);
ACC = mean(pred==labels)*100;
TP = sum(pred==1 & labels==1);
TN = sum(pred==0 & labels==0);
FP = sum(pred==1 & labels==0);
FN = sum(pred==0 & labels==1);
SEN = TP / max(TP+FN,1);
SPC = TN / max(TN+FP,1);
PPV = TP / max(TP+FP,1);
NPV = TN / max(TN+FN,1);
OR  = (TP*TN) / max(FP*FN,1);
LRp = SEN / max(1-SPC, eps);

% Requires Statistics and Machine Learning Toolbox
try
    [~,~,~,AUC] = perfcurve(labels, scores, 1);
catch
    AUC = NaN;
end

metrics = struct('ACC',ACC,'AUC',AUC,'SEN',SEN,'SPC',SPC,'PPV',PPV,'NPV',NPV,'OR',OR,'LRp',LRp,...
    'TP',TP,'TN',TN,'FP',FP,'FN',FN);
end
