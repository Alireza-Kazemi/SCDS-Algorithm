%% Sub Clustering With Feature Selection
% Labels are Wake:1 / NREM1:2 / NREM2:3 / NREM3:4 / REM:5 
clear;
clc;
% close all;

%% Initialization 

precision = @(confusionMat) diag(confusionMat)./sum(confusionMat,2);
recall = @(confusionMat) diag(confusionMat)./sum(confusionMat,1)';
f1Scores = @(confusionMat) 2*(precision(confusionMat).*recall(confusionMat))./(precision(confusionMat)+recall(confusionMat)+(1-sign(recall(confusionMat))));
meanF1 = @(confusionMat) mean(f1Scores(confusionMat),'omitnan');

Fnum = 62;
Subsetnum = 19;
Classnum=5;
EpochTime = 30;
Nfilt = 0;



load EEGFeatAll_150lag_MNNormalized.mat

Indsorig = 1:Subsetnum; randperm(39,Subsetnum);%31; Works awesome for 31st participant
ACCFPZ = zeros(7,Subsetnum);
ACCPZ = zeros(7,Subsetnum);
Kappa = zeros(2,Subsetnum);
F1score = zeros(2,Subsetnum);
ACCBoth = zeros(7,Subsetnum);
Mylabels = cell(2,Subsetnum);
for itt = 1:Subsetnum
%     load FeatIDXnew_Jun2020.mat
    load FeatIDXnew.mat
    Channel = 1;
    FeatIDX = FeatIDX08(Channel,:);%{1:11,1:11,1:11,1:11,1:11};%FeatIDXleast(Channel,:);%
    Centers = Centers(Channel,:);
    for i=1:5
        Centers{i}=Centers{i}(FeatIDX{i},:);
    end
    clc
    Inds = Indsorig(itt);
    FeatData = Data.EEG;
    Hyp = [];
    Feats = [];
    disp(Inds)
    for i=Inds
        Hyp = cat(1,Hyp,Data.Hyp{i}(1:end));
        Feats = cat(1,Feats,FeatData{i});
    end
    Hyp(Hyp==4)=3;
    Hyp(Hyp==5)=4;
    Hyp(Hyp==6)=5;
    Hyp(Hyp==0)=6;
    Featsall = Feats;
    Feats = Featsall(:,1:Fnum);
    LabelString ={ 'NREM1', 'NREM2', 'NREM3', 'REM','Wake'};

    %% SubClustering
    disp('=========================================FPz-Cz')
    Xnum = length(Feats(:,1));
    Feats = (Feats - repmat(min(Feats),Xnum,1))./(repmat(max(Feats),Xnum,1)-repmat(min(Feats),Xnum,1));

    idx1 = zeros(Classnum,Xnum);
    DataInd1 = 1:Xnum;
    for Classind = [5,4,3,2,1]%[1,2,5,4,3]%[1,4,3,5,2]%Optimum Order [5,4,3,2,1]
        Cents = Centers{Classind};
        Features = Feats(DataInd1,:);
        Features = Features(:,FeatIDX{Classind});
        opts = statset('MaxIter',1000);
        if (size(Features,1)<2)
                break;
        end
        tmp1 = kmeans(Features,2,'Start',Cents','MaxIter',1000);%
        tmp1(tmp1==2) = 0;
        tmp1(tmp1==1) = Classind;
        idx1(Classind,DataInd1) = tmp1;
        DataInd1(tmp1~=0)=[];
    end

    Labtemp = sum(idx1);
    HypextFPz = Labtemp;
    Mylabels{1,itt} = HypextFPz;
    L1 = zeros(Classnum,Xnum*EpochTime);
    for i=1:Classnum
        idx1tmp = repmat(idx1(i,:)',1,EpochTime);
        if (Nfilt>0)
            L1(i,:) = round(medfilt1(reshape(idx1tmp',1,EpochTime*length(idx1tmp(:,1))),Nfilt));      
        else
            L1(i,:) = reshape(idx1tmp',1,EpochTime*length(idx1tmp(:,1)));
        end
    end

    idx1 = sum(L1);

   %% Results
   Labels = {'NREM1',...
              'NREM2',...
              'NREM3',...
              'REM',...
              'Wake'};

    Hyp(Hyp==6)=0;
    Inds = and((Labtemp~=0),(Hyp~=0)');
    Unrecognized = (1-length(find(Inds))/length(Inds))*100;
    disp("Unrecognized(%):")
    disp(Unrecognized)
    cm = confusionmat(Hyp(Inds),Labtemp(Inds));
    disp(cm)
    Accuracy = sum(diag(cm))/sum(sum(cm))*100;
    disp("Accuracy(%):")
    disp(Accuracy)
    AccuracyStages = diag(cm)'./sum(cm)*100;
    disp(AccuracyStages)
    ACCFPZ(3:end,itt) = AccuracyStages;
    ACCFPZ(2,itt) = Accuracy;
    ACCFPZ(1,itt) = Unrecognized;
    Kappa(Channel,itt) = kappa(cm);
    Kappa(Channel+2,itt) = kappaindex(Labtemp(Inds),Hyp(Inds),5);
    F1score(Channel,itt) = meanF1(cm);
    %% Channel Pz-Oz
    %#########################################################################
    %#########################################################################
    %#########################################################################
    %#########################################################################
    %#########################################################################
    %#########################################################################
    %#########################################################################
    %#########################################################################
    %#########################################################################
    %% Initialization 
%     load FeatIDXnew_Jun2020.mat
    load FeatIDXnew.mat
    Hyp(Hyp==0)=6;
    Channel = 2;
    FeatIDX = FeatIDX08(Channel,:);%{1:11,1:11,1:11,1:11,1:11};%FeatIDXleast(Channel,:);%
    Centers = Centers(Channel,:);
    for i=1:5
        Centers{i}=Centers{i}(FeatIDX{i},:);
    end

    Feats = Featsall(:,Fnum+1:end);
    
    %% SubClustering
    disp('=========================================Pz-Oz')
    Xnum = length(Feats(:,1));
    Feats = (Feats - repmat(min(Feats),Xnum,1))./(repmat(max(Feats),Xnum,1)-repmat(min(Feats),Xnum,1));

    idx1 = zeros(Classnum,Xnum);
    DataInd1 = 1:Xnum;
    for Classind = [5,4,3,2,1]%[1,2,5,4,3]%[1,4,3,5,2]%Optimum Order [5,4,3,2,1]
        Cents = Centers{Classind};
        Features = Feats(DataInd1,:);
        Features = Features(:,FeatIDX{Classind});
        opts = statset('MaxIter',1000);
        if (size(Features,1)<2)
                break;
        end
        tmp1 = kmeans(Features,2,'Start',Cents','MaxIter',1000);%
        tmp1(tmp1==2) = 0;
        tmp1(tmp1==1) = Classind;

        idx1(Classind,DataInd1) = tmp1;
        DataInd1(tmp1~=0)=[];

    end

    Labtemp = sum(idx1);
    HypextPz = Labtemp;
    Mylabels{2,itt} = HypextPz;
    L1 = zeros(Classnum,Xnum*EpochTime);
    for i=1:Classnum
        idx1tmp = repmat(idx1(i,:)',1,EpochTime);
        if (Nfilt>0)
            L1(i,:) = round(medfilt1(reshape(idx1tmp',1,EpochTime*length(idx1tmp(:,1))),Nfilt));      
        else
            L1(i,:) = reshape(idx1tmp',1,EpochTime*length(idx1tmp(:,1)));
        end
    end

    idx1 = sum(L1);

    %% Results
    
    Hyp(Hyp==6)=0;
    Inds = and((Labtemp~=0),(Hyp~=0)');
    Unrecognized = (1-length(find(Inds))/length(Inds))*100;
    disp("Unrecognized(%):")
    disp(Unrecognized)
    cm = confusionmat(Hyp(Inds),Labtemp(Inds));
    disp(cm)
    Accuracy = sum(diag(cm))/sum(sum(cm))*100;
    disp("Accuracy(%):")
    disp(Accuracy)
    AccuracyStages = diag(cm)'./sum(cm)*100;
    disp(AccuracyStages)
    ACCPZ(3:end,itt) = AccuracyStages;
    ACCPZ(2,itt) = Accuracy;
    ACCPZ(1,itt) = Unrecognized;
    %% Compare two channels
    Diffs = HypextPz - HypextFPz;
    inds = find(Diffs);
    % Hyp(inds) = [];
    HypextPz(inds) = 0;
    HypextFPz(inds) = 0;
    Labtemp = HypextFPz;
    disp("=========================================Both")
    Inds = and((Labtemp~=0),(Hyp~=0)');
    Unrecognized = (1-length(find(Inds))/length(Inds))*100;
    disp("Unrecognized(%):")
    disp(Unrecognized)
    cm = confusionmat(Hyp(Inds),HypextPz(Inds));
    disp(cm)
    Accuracy = sum(diag(cm))/sum(sum(cm))*100;
    disp("Accuracy(%):")
    disp(Accuracy)
    AccuracyStages = diag(cm)'./sum(cm)*100;
    AccuracyStages(isnan(AccuracyStages))=0;
    if(length(AccuracyStages)<5)
        acctemp = zeros(1,5);
        acctemp(unique(Hyp(Inds))) = AccuracyStages;
        AccuracyStages = acctemp;
    end
    disp(AccuracyStages)
    ACCBoth(3:end,itt) = AccuracyStages;
    ACCBoth(2,itt) = Accuracy;
    ACCBoth(1,itt) = Unrecognized;
    Kappa(Channel,itt) = kappa(cm);
    Kappa(Channel+2,itt) = kappaindex(Labtemp(Inds),Hyp(Inds),5);
    F1score(Channel,itt) = meanF1(cm);
end


%% plots

close all


figure
Unrec = cat(1,ACCFPZ(1,:),ACCPZ(1,:),ACCBoth(1,:))';
boxplot(Unrec)
ylim([0,100])
LabelString ={'FPz-Cz', 'Pz-Oz', 'Both'};
ylabel("Average percent(%)")
xticklabels(LabelString)
title("Percent of Unrecognized Epochs")


figure
FPZ = cat(1,ACCFPZ(2,:)/100,Kappa(1,:),F1score(1,:))';
PZ  = cat(1,ACCPZ(2,:)/100,Kappa(2,:),F1score(2,:))';
x = (1:3)-.2;
y = mean(FPZ,'omitnan');
err = std(FPZ,'omitnan');
errorbar(x,y,err,'x','Linewidth',1.2)
hold on
x = (1:3)+.2;
y = mean(PZ,'omitnan');
err = std(PZ,'omitnan');
errorbar(x,y,err,'x','Linewidth',1.2)
ylim([0,1])
xticks(1:3)
LabelString ={'FPz-Cz', 'Pz-Oz'};
xticklabels({'Accuracy','Cohen''s Kappa','F1-Score'})
ylabel("Performance")
legend(LabelString,'Location','NorthWest')



figure
LabelString ={'Overal', 'NREM1', 'NREM2', 'NREM3', 'REM','Wake'};
x = (1:6)-.2;
y = mean(ACCFPZ(2:end,:)','omitnan');
err = .5*std(ACCFPZ(2:end,:)','omitnan');
errorbar(x,y,err,'x','Linewidth',1.2)
hold on
x = 1:6;
y = mean(ACCPZ(2:end,:)','omitnan');
err = .5*std(ACCPZ(2:end,:)','omitnan');
errorbar(x,y,err,'x','Linewidth',1.2)

x = (1:6)+.2;
y = mean(ACCBoth(2:end,:)','omitnan');
err = .5*std(ACCBoth(2:end,:)','omitnan');
errorbar(x,y,err,'x','Linewidth',1.2)
ylim([0,100])
xticks(1:6)
xticklabels(LabelString)
ylabel("Average Accuracy (percent correct)")
legend({'FPz-Cz','Pz-Oz','Both'},'Location','NorthWest')

figure
LabelString ={'Overal', 'NREM1', 'NREM2', 'NREM3', 'REM','Wake'};
x = (1:6)-.2;
y = mean(ACCFPZ(2:end,:)','omitnan');
err = 1/sqrt(19)*std(ACCFPZ(2:end,:)','omitnan');
errorbar(x,y,err,'x','Linewidth',1.2)
hold on
x = 1:6+0.2;
y = mean(ACCPZ(2:end,:)','omitnan');
err = 1/sqrt(19)*std(ACCPZ(2:end,:)','omitnan');
errorbar(x,y,err,'x','Linewidth',1.2)

ylim([0,100])
xlim([0,7])
xticks(1:6)
xticklabels(LabelString)
ylabel("Average Accuracy (percent correct)")
legend({'FPz-Cz','Pz-Oz'},'Location','NorthWest')