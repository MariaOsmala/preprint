
%http://www0.cs.ucl.ac.uk/staff/ucacbbl/roc/

%In practice to use a classifier one normally has to chose an operating point, a threshold. This fixes a point on the ROC. 
%In some cases the area may be misleading. That is, when comparing classifiers, the one with the larger area may not be the 
%one with the better performance at the chosen threshold (or limited range).

%The ROC curve can be used to choose the best operating point. The best operating point might be chosen so that the 
%classifier gives the best trade off between the costs of failing to detect positives against the costs of raising false alarms. 
%These costs need not be equal, however this is a common assumption.

%If the costs associated with classification are a simple sum of the cost of misclassifying positive and negative cases
%then all points on a straight line (whose gradient is given by the importance of the positive and negative examples) 
%have the same cost. If the cost of misclassifying positive and negative cases are the same, and positive and negative cases 
%occur equally often then the line has a slope of 1. That is ,it is at 45 degrees. The best 
%place to operate the classifier is the point on its ROC which lies on a 45 degree line closest to the north-west corner (0,1) of the ROC plot.
clear all;
path_to_dir='/m/cs/scratch/csb/projects/enhancer_prediction/experiments/RProjects/preprint/'
comparisons_path=strcat(path_to_dir, '/results/model_promoters_and_random_combined/K562/');
results=cell(2*2*3,3);

N=1000;
window=2000;
bin_size=100;


results_ROC=cell(2,2,3);

i=1;
%%

for distance_measure={'ML', 'Bayes_estimated_priors'} %'correlation'
    
    
    path=strcat(comparisons_path,'/',distance_measure);
    format_rep=45;
    sample_number=4000;
    
    clear true_labels
    clear predictions
    for k=1:5
        
        filename=strcat(path,'/5-fold_CV_',num2str(k),'/NSamples_',num2str(N),'_window_',...
            num2str(window),'_bin_',num2str(bin_size),'_5fold_cv_',num2str(k),'_test_data.txt.predict');
        
        
        predicted=importdata(filename{1}, ' ', 1);
        predictions( ((k-1)*(sample_number/5)+1):(k*(sample_number/5)) )=predicted.data(:,2);
        
                filename=strcat(path,'/5-fold_CV_',num2str(k),'/NSamples_',num2str(N),'_window_',...
                    num2str(window),'_bin_',num2str(bin_size),'_5fold_cv_',num2str(k),'_test_data.txt');

        fid=fopen(filename{1});
        format=strcat('%d8', strjoin(repmat({'%d8:%f'},format_rep,1)'  ,''));
        testdata=textscan(fid, format);
        fclose(fid);
        true_labels( ((k-1)*(sample_number/5)+1):(k*(sample_number/5)) )=testdata{1,1};
        
    end
    
    [X,Y,T,AUC,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(true_labels,predictions,1);
    
    results_ROC{i,1}=X;
    results_ROC{i,2}=Y;
    %X FPR
    %Y TPR
    
    
    threshold_0_01=T(find(X==0.01,1));
    TPR_0_01=Y(find(X==0.01,1));
    
    X_coord=find(X==OPTROCPT(1));
    Y_coord=find(Y==OPTROCPT(2));
    
   
    [c,m]=hist([X_coord; Y_coord], (min([X_coord; Y_coord]):1:max([X_coord; Y_coord])));
    thres_ind=m(find(c==2));
    
    tbop=T(thres_ind);
    FPR=OPTROCPT(1);
    TPR=OPTROCPT(2);
    %AUC
    
    results{i}=strjoin([ distance_measure, num2str(tbop), num2str(FPR), num2str(TPR), num2str(AUC),   num2str(threshold_0_01), num2str(TPR_0_01)]);
    i=i+1;
end

%% Plot ROC


figure(1)
wx = 10; wy = 10; set(gcf,'PaperSize',[wx wy], 'PaperPosition',[0 0 wx wy]); %centimeters
plot(results_ROC{1,1}, results_ROC{1,2}, 'color','red','lineWidth',2)
hold on
plot(results_ROC{2,1}, results_ROC{2,2}, 'color','blue','lineWidth',2 )
legend('ML','Bayes')
title(strcat('K562, 5-fold CV in the training data, combined'))
print(strcat('K562_trainingCV_ROC_combined.png'),'-dpng', '-r300');

%%

clear all;
path_to_dir='/m/cs/scratch/csb/projects/enhancer_prediction/experiments/RProjects/preprint'
comparisons_path=strcat(path_to_dir,'/results/model_promoters_and_random_combined/GM12878');


results=cell(2*2*2,3); %2*2*3
i=1;

N=1000;
window=2000;
bin_size=100;

results_ROC=cell(2,2,3);
%%


for distance_measure={'ML', 'Bayes_estimated_priors'}
    
    
    path=strcat(comparisons_path,'/',distance_measure);
    
    format_rep=45;
    sample_number=4000;
    
    
    clear true_labels
    clear predictions
    
    
    filename=strcat(path,'/NSamples_1000_window_2000_bin_100_test_data.txt.predict');
    predicted=importdata(filename{1}, ' ', 1);
  
    predictions=predicted.data(:,2);
    
    filename=strcat(path,'/NSamples_1000_window_2000_bin_100_test_data.txt');
    
    
    fid=fopen(filename{1});
    format=strcat('%d8', strjoin(repmat({'%d8:%f'},format_rep,1)'  ,''));
    testdata=textscan(fid, format);
    fclose(fid);
    true_labels=testdata{1,1};
    
    
    
    
    [X,Y,T,AUC,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(true_labels,predictions,1);
    
    results_ROC{i,1}=X;
    results_ROC{i,2}=Y;
    %X FPR
    %Y TPR
    
    threshold_0_01=T(find(X==0.01,1)); %0.0095
    
    
    
    TPR_0_01=Y(find(X==0.01,1));
    
   
    
    
  
    
    
    
    
    
    X_coord=find(X==OPTROCPT(1));
    Y_coord=find(Y==OPTROCPT(2));
    
    
    
    
    [c,m]=hist([X_coord; Y_coord], (min([X_coord; Y_coord]):1:max([X_coord; Y_coord])));
    thres_ind=m(find(c==2));
    
    tbop=T(thres_ind);
    FPR=OPTROCPT(1);
    TPR=OPTROCPT(2);
    %AUC
    
    results{i}=strjoin([ distance_measure, num2str(tbop), num2str(FPR), num2str(TPR), num2str(AUC), num2str(threshold_0_01), num2str(TPR_0_01)]);
    i=i+1;

end
%% Load GM12878 predictions and true labels for the test data


filename=strcat(path_to_dir,'/results/RFECS_combined/GM12878/predictions.txt');
 
predicted=importdata(filename,' ',0);
predictions=predicted;
    
filename=strcat(path_to_dir,'/results/RFECS_combined/GM12878/true_labels.txt');
true_labels=importdata(filename,' ',0);
[X_RFECS,Y_RFECS,T_RFECS,AUC_RFECS,OPTROCPT_RFECS,SUBY_RFECS,SUBYNAMES_RFECS] = perfcurve(true_labels,predictions,1);

% Plot ROC
% The returned values X and Y
%are coordinates for the performance curve and can be visualized  with
%PLOT(X,Y). By default, X is false positive rate, FPR and Y is true positive rate, TPR
%results_ROC{i,1}=X;
%results_ROC{i,2}=Y;
figure
wx = 10; wy = 10; set(gcf,'PaperSize',[wx wy], 'PaperPosition',[0 0 wx wy]); %centimeters
plot(results_ROC{1,1,1}, results_ROC{1,2,1}, 'color','red','lineWidth',2)
hold on
plot(results_ROC{2,1,1}, results_ROC{2,2,1}, 'color','blue','lineWidth',2 )
hold on
plot(X_RFECS, Y_RFECS, 'color','green','lineWidth',2 )

legend('ML','Bayes', 'RFECS')

title(strcat('GM12878, test data'))


print(strcat('GM12878_test_ROC_combined.png'),'-dpng', '-r300');



%%
