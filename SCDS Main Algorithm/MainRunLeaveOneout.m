%% Initialization
clear;
clc;
close all;
set(groot,'defaultfigureposition',[100 100 700 900])
ColorCodes = ['#0072BD';'#D95319';'#7E2F8E';'#77AC30';'#4DBEEE';'#A2142F'];

PlotFlag = 0;
plotDifferencesFlag = 0;
PlotHypFlag = 1;
PlotCorrelationFlag = 0;

Flag_RemoveStageWise = 0;
Flag_RemoveStageWise_Center = 0;

Flag_GradedCorrelation = 0; GradedCodingVal = 0.5;
Flag_AllSubjects = 0;


load EEGFeatEDF_MNNormalized.mat
load FeatLabels.mat

Cnum = 5;
FeatData = Data.EEG;
Ch_num = 2;
Fnum = 62;
FnumTot = Fnum*Ch_num;
% for i=1:39
%     FeatData{i} = FeatData{i}(:,Fnum+1:124);
% end

MDF_Fnum = [15,13,11,11,50];
ACCs_CV_MDF = [];
ACCs_CV_LDF = [];
ACCs_CV_CFF = [];
ACC_KNN = [];
ACC_RF = [];
ACC_KNN_MDF = [];
ACC_RF_MDF = [];
FeatIDX_CV = cell(1,5);
FeatIDX_CV_MDF = [];

IDsForTrain = [0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,14,14,15,15,16,16,17,17,18,18,19,19];
OriginalID = 1:39;
for TestID=0:19 %unique(IDsForTrain)
% for itt=0:20 %unique(IDsForTrain)
    
    
% IDstotest = randperm(20,15)-1;
IDstotest = TestID;
TrainIDs = OriginalID(~ismember(IDsForTrain,IDstotest));
TestIDs = OriginalID(ismember(IDsForTrain,IDstotest));
SubsetIds = TrainIDs;
Optimum_order = [5,3,2,4,1];
% Optimum_order = [5,4,3,1,2]; % alternative Order

ColorCodes = ['#0072BD';'#D95319';'#7E2F8E';'#77AC30';'#4DBEEE';'#A2142F'];
ColorCodeOrange = '#EDB120';
Ch_Names = ["Fpz-Cz","Pz-Oz"];

SleepStages = {'NREM1',...
               'NREM2',...
               'NREM3',...
               'REM',...
               'Wake'};

CH_Lab_total = [];
for i=1:Ch_num
    CH_Lab_total = cat(2,CH_Lab_total,CH_Lab);
end
indx = 1;
for i=1:Ch_num
    str = "_"+Ch_Names(i);
    for j=1:length(CH_Lab)
        CH_Lab_total{indx} = CH_Lab_total{indx}+str;
        indx = indx+1;
    end
end

%% Training

%% Discriminative Subspaces
DiscriminativeSubspaces_BothChannels;
for cidx=1:Cnum
    FeatIDX_CV{cidx}=[FeatIDX_CV{cidx};IDX{cidx}'];
end
FeatIDX_CV_MDF = cat(1,FeatIDX_CV_MDF,FeatIDXMDF);
%% Supervised Methods
TrainSupervised;
%% Testing
SubsetIds = TestIDs;
TestSupervised;
ACC_KNN = [ACC_KNN;Accuracy_KNN,AccuracyStages_KNN];
ACC_KNN_MDF = [ACC_KNN_MDF;Accuracy_KNN_MDF,AccuracyStages_KNN_MDF];
ACC_RF = [ACC_RF;Accuracy_RF,AccuracyStages_RF];
ACC_RF_MDF = [ACC_RF_MDF;Accuracy_RF_MDF,AccuracyStages_RF_MDF];

%% Sub Clustering With Feature Selection
% Labels are Wake:1 / NREM1:2 / NREM2:3 / NREM3:4 / REM:5 


disp('===========================================================')
disp('============================MDF============================')
SCDS_5stages_BothChannels
ACCs_CV_MDF = [ACCs_CV_MDF;mean(ACCsall(2:end,:),2)'];

% disp('===========================================================')
% disp('============================LDF============================')
% SCDS_5stages_BothChannels_LDF
% ACCs_CV_LDF = [ACCs_CV_LDF;mean(ACCsall(2:end,:),2)'];
% 
% disp('===========================================================')
% disp('============================CFF============================')
% SCDS_5stages_BothChannels_CFF
% ACCs_CV_CFF = [ACCs_CV_CFF;mean(ACCsall(2:end,:),2)'];
% close all
end
ACC = cat(1,ACC_KNN, ACC_RF, ACC_KNN_MDF, ACC_RF_MDF,...
            ACCs_CV_MDF, ACCs_CV_LDF, ACCs_CV_CFF);
L = length(ACC_KNN(:,1));
Method = cat(1,repmat("KNN",L,1),repmat("RF",L,1),...
               repmat("KNN_MDF",L,1),repmat("RF_MDF",L,1),...
               repmat("SCDS_MDF",L,1),repmat("SCDS_LDF",L,1),repmat("SCDS_CFF",L,1));
ACC_tbl = table(Method,ACC(:,1),ACC(:,2),ACC(:,3),ACC(:,4),ACC(:,5),ACC(:,6),...
            'VariableNames',{'Method','Overal','NREM1','NREM2','NREM3','REM','WAKE'});
write(ACC_tbl,'Accuracy_results_oneIn2.csv')