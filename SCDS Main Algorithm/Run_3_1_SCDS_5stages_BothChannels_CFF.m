%	This file is part of SCDS Algorithm.
%
%    SCDS Algorithm is free: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    any later version.
%
%    SCDS Algorithm is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
% Designed and developed by Alireza Kazemi 2020
% Address your comments and questions to alireza.kzmi@gmail.com

%% Sub Clustering With Feature Selection
% Labels are Wake:1 / NREM1:2 / NREM2:3 / NREM3:4 / REM:5 
clear;
clc;
close all;
set(groot,'defaultfigureposition',[100 100 700 900])

%% Inline functions:
precision = @(confusionMat) diag(confusionMat)./sum(confusionMat,2);
recall = @(confusionMat) diag(confusionMat)./sum(confusionMat,1)';
f1Scores = @(confusionMat) 2*(precision(confusionMat).*recall(confusionMat))./(precision(confusionMat)+recall(confusionMat)+(1-sign(recall(confusionMat))));
meanF1 = @(confusionMat) mean(f1Scores(confusionMat),'omitnan');

%% Initialization 

ChNum = 2;
Fnum = 62*ChNum;
Subsetnum = 39;
Classnum=5;
EpochTime = 30;
Nfilt = 0;

ColorCodes = ['#0072BD';'#D95319';'#7E2F8E';'#77AC30';'#4DBEEE';'#A2142F'];
 
load EEGFeatEDF_MNNormalized.mat
load FeatLabels.mat
load FeatIDXBoth.mat


% FeatIDX = FeatIDX2;
FeatIDX = ;

for i=1:Classnum
    MDF_Centers{i}=Centers{i}(FeatIDX{i},:);
end
% FeatIDX = FeatIDX1;
    
Indsorig =1:39; randperm(39,Subsetnum);
ACCsall = zeros(7,Subsetnum);
ACCs3 = zeros(4,Subsetnum);
Kappa = zeros(1,Subsetnum);
Kappa2 = zeros(1,Subsetnum);
F1score = zeros(1,Subsetnum);
Mylabels = cell(1,Subsetnum);
%% Itteration
for itt = 1:Subsetnum
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
    LabelString ={ 'NREM1', 'NREM2', 'NREM3', 'REM','Wake'};

    %% SubClustering
    disp('=========================================Both Channels')
    Xnum = length(Feats(:,1));

    idx1 = zeros(Classnum,Xnum);
    DataInd1 = 1:Xnum;
    for Classind = Optimum_order
        Cents = MDF_Centers{Classind};
        Features = Feats(DataInd1,:);
        Features = Features(:,FeatIDX{Classind});
        opts = statset('MaxIter',1000);
        if (size(Features,1)<2)
            break;
        end
        tmp1 = kmeans(Features,2,'Start',Cents','MaxIter',1000);
        
        tmp1(tmp1==2) = 0;
        tmp1(tmp1==1) = Classind;

        idx1(Classind,DataInd1) = tmp1;
        DataInd1(tmp1~=0)=[];
    end
    Mylabels{itt} = sum(idx1);

    %% Plots
    Labelstr = {'NREM1',...
                  'NREM2',...
                  'NREM3',...
                  'REM',...
                  'Wake'};

    
    figure;
    subplot(2,1,1)
    [H,C]= hist(Mylabels{itt},1:Classnum+1);
    bar(C,H/30,.7);
    xlabel('Cluster Index')
    ylabel('Number of Members')    
    set(gca,'XTickLabel',Labelstr)

    subplot(2,1,2)
    [H,C]= hist(Hyp',1:Classnum);
    bar(C,H/30,.7);
    xlabel('Cluster Index')
    ylabel('Number of Members')
    set(gca,'XTickLabel',Labelstr(1:5))

    figure
    subplot(2,1,1)
    tmp = Mylabels{itt};
    Inds = 1:length(tmp);
    title('Kmeans')    
    tmp(tmp==0) = 2.5;
    Flag = ~isempty(find(tmp==2.5,1));
    H = plot(Inds(tmp==2.5),tmp(tmp==2.5),'*k');
    hold on
    Inds(tmp==2.5)=[];
    tmp(tmp==2.5)=[];
    DiscreteFlatPlot(Inds,tmp,Labelstr)
    title('Clustered Hypmogram')
    
    subplot(2,1,2)
    Inds = 1:length(Hyp);
    DiscreteFlatPlot(Inds,Hyp',Labelstr)
    xlabel('Time (s)');
    title('Real Hypmogram')
    
    %% Results
    Hyp(Hyp==6)=0;
    Labtemp = Mylabels{itt};
    Inds = and((Labtemp~=0),(Hyp~=0)');
    Unrecognized = (1-length(find(Inds))/length(Inds))*100;
    disp("Unrecognized(%):")
    disp(Unrecognized)
    cm = confusionmat(Hyp(Inds),Labtemp(Inds));
    disp(cm)
    Accuracy = sum(diag(cm))/sum(sum(cm))*100;
    disp("Accuracy all(%):")
    disp(Accuracy)
    AccuracyStages = diag(cm)'./sum(cm)*100;
    disp(AccuracyStages)
    ACCsall(3:end,itt) = AccuracyStages;
    ACCsall(2,itt) = Accuracy;
    ACCsall(1,itt) = Unrecognized;
    Kappa(itt) = kappa(cm);
    Kappa2(itt) = kappaindex(Labtemp(Inds),Hyp(Inds),5);
    F1score(itt) = meanF1(cm);
    
    
    Hyp(Hyp==1)=3;
    Hyp(Hyp==2)=3;
    Hyp(Hyp==3)=3;
    Labtemp(Labtemp==1)=3;
    Labtemp(Labtemp==2)=3;
    Labtemp(Labtemp==3)=3;
    cm = confusionmat(Hyp(Inds),Labtemp(Inds));
    disp(cm)
    Accuracy = sum(diag(cm))/sum(sum(cm))*100;
    disp("Accuracy 3(%):")
    disp(Accuracy)
    AccuracyStages = diag(cm)'./sum(cm)*100;
    ACCs3(2:end,itt) = AccuracyStages;
    ACCs3(1,itt) = Accuracy;
end
%% Plot final results
clc
close all;
LabelString ={'Overal', 'NREM1', 'NREM2', 'NREM3', 'REM','Wake'};

figure
x = (1:6);
y = mean(ACCsall(2:end,:)','omitnan');
err = 1/sqrt(size(ACCsall,2))*std(ACCsall(2:end,:)','omitnan');
disp(LabelString)
disp(x)
disp(y)
errorbar(x,y,err,'o','Linewidth',1.2)
ylim([0,100])
xticks(1:6)
xticklabels(LabelString)
xlim([0,7])
ylabel("Average Accuracy (percent correct)")

figure
for i=2:7
    subplot(6,1,i-1)
    plot(Indsorig,ACCsall(i,:),'-.','Color',ColorCodes(i-1,:))
    hold on
    plot(Indsorig,ACCsall(i,:),'*','Color',ColorCodes(i-1,:))
    ylabel(LabelString{i-1})
    ylim([0,100])
end
xlabel('Participants'' ID')

%---------- only 3 stages
LabelString ={'Overal', 'NREM', 'REM', 'Wake'};

figure
x = (1:4);
y = mean(ACCs3','omitnan');
err = 1/sqrt(size(ACCs3,2))*std(ACCs3','omitnan');
errorbar(x,y,err,'o','Linewidth',1.2)
ylim([0,100])
xticks(1:4)
xticklabels(LabelString)
xlim([0,5])
ylabel("Average Accuracy (percent correct)")
