%% Initialization
clear;
clc;
close all;
set(groot,'defaultfigureposition',[100 100 700 900])
ColorCodes = ['#0072BD';'#D95319';'#7E2F8E';'#77AC30';'#4DBEEE';'#A2142F'];

PlotFlag = 0;
plotDifferencesFlag = 1;
PlotHypFlag = 1;
PlotCorrelationFlag = 0;

Flag_RemoveStageWise = 0;
Flag_RemoveStageWise_Center = 0;

Flag_GradedCorrelation = 0; GradedCodingVal = 0.5;
Flag_AllSubjects = 1;


load EEGFeatEDF_MNNormalized.mat
load FeatLabels.mat

Cnum = 5;
FeatData = Data.EEG;
Ch_num = 2;
Fnum = 62;
FnumTot = Fnum*Ch_num;
MDF_Fnum = [3,16,11,12,11];%[3,16,11,12,11]
% [1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1];
IDsForTrain = [1,2,1,1,2,2,1,2,1,1,1,1,1,2,1,2,2,1,2,2,1,1,1,2,2,1,2,1,2,2,2,1,1,2,1,2,2,1,2];%repmat(2,1,39);%
OriginalID = 1:39;
% For itterative exploration of centers
% Subsetnum = 39;
SubsetIds = 1:39;%OriginalID(IDsForTrain==2);
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

%%
HypDat = [];
FeatDat = [];
for i=SubsetIds
    HypDat = cat(1,HypDat,Data.Hyp{i});
    FeatDat = cat(1,FeatDat,FeatData{i});
end
HypDat(HypDat==4)=3;
HypDat(HypDat==5)=4;
HypDat(HypDat==6)=5;
HypDat(HypDat==0)=6;
%% Discriminative Subspaces
DiscriminativeSubspaces_BothChannels;
%% Sub Clustering With Feature Selection
% Labels are Wake:1 / NREM1:2 / NREM2:3 / NREM3:4 / REM:5 
SubsetIds = OriginalID(IDsForTrain==2);
disp('===========================================================')
disp('============================MDF============================')
SCDS_5stages_BothChannels
%% Sub Clustering With Feature Selection
% Labels are Wake:1 / NREM1:2 / NREM2:3 / NREM3:4 / REM:5 
disp('===========================================================')
disp('============================LDF============================')
SCDS_5stages_BothChannels_LDF
%% Sub Clustering With Feature Selection
% Labels are Wake:1 / NREM1:2 / NREM2:3 / NREM3:4 / REM:5 
disp('===========================================================')
disp('============================CFF============================')
SCDS_5stages_BothChannels_CFF