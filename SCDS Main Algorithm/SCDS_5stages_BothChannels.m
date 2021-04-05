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


%% Inline functions:
precision = @(confusionMat) diag(confusionMat)./sum(confusionMat,2);
recall = @(confusionMat) diag(confusionMat)./sum(confusionMat,1)';
f1Scores = @(confusionMat) 2*(precision(confusionMat).*recall(confusionMat))./(precision(confusionMat)+recall(confusionMat)+(1-sign(recall(confusionMat))));
meanF1 = @(confusionMat) mean(f1Scores(confusionMat),'omitnan');

%% Initialization 

% FeatIDX = FeatIDX2;
FeatIDX = FeatIDXMDF;

for i=1:Cnum
    MDF_Centers{i}=Centers{i}(FeatIDX{i},:);
end
    
% SubsetIds = OriginalID(IDsForTrain==2);

Subsetnum = length(SubsetIds);
ACCsall = zeros(7,Subsetnum);
ACCs3 = zeros(4,Subsetnum);
Kappa = zeros(1,Subsetnum);
Kappa2 = zeros(1,Subsetnum);
F1score = zeros(1,Subsetnum);
Mylabels = cell(1,Subsetnum);
%% Itteration
for itt = 1:Subsetnum
    Inds = SubsetIds(itt);
    Hyp = [];
    Feats = [];
%     disp(Inds)
    for i=Inds
        Hyp = cat(1,Hyp,Data.Hyp{i});
        Feats = cat(1,Feats,FeatData{i});
    end
    Hyp(Hyp==4)=3;
    Hyp(Hyp==5)=4;
    Hyp(Hyp==6)=5;

    %% SubClustering
%     disp('=========================================Both Channels')
    Xnum = length(Feats(:,1));

    idx1 = zeros(Cnum,Xnum);
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
    if(PlotHypFlag == 1)
        Labelstr = {'NREM1',...
                      'NREM2',...
                      'NREM3',...
                      'REM',...
                      'Wake'};

% 
%         figure;
%         subplot(2,1,1)
%         [H,C]= hist(Mylabels{itt},1:Cnum+1);
%         bar(C,H/30,.7);
%         xlabel('Cluster Index')
%         ylabel('Number of Members')    
%         set(gca,'XTickLabel',Labelstr)
% 
%         subplot(2,1,2)
%         [H,C]= hist(Hyp',1:Cnum);
%         bar(C,H/30,.7);
%         xlabel('Cluster Index')
%         ylabel('Number of Members')
%         set(gca,'XTickLabel',Labelstr(1:5))

        tmp = Mylabels{itt};
        tmp = repmat(tmp,30,1);
        tmp = reshape(tmp,1,numel(tmp));
        tmpH = repmat(Hyp',30,1);
        tmpH = reshape(tmpH,1,numel(tmpH));
        t = (0:(length(tmpH)-1))/3600;
        IndU = (tmp==0);
        tmp(IndU) = tmpH(IndU);
        figure
        subplot(2,1,1)
        plot(t,tmp);
        hold on
        plot(t(IndU),tmp(IndU),'*k');
        ylim([0,6])
        yticks([0,1,2,3,4,5,6])
        yticklabels({' ','N1','N2','N3','REM','W',' '})
        xlim([0,t(end)])
        title('Clustered Hypmogram')
        
        subplot(2,1,2)
        plot(t,tmpH)
        ylim([0,6])
        yticks([0,1,2,3,4,5,6])
        yticklabels({' ','N1','N2','N3','REM','W',' '})
        xlim([0,t(end)])
        xlabel('Time (Hour)');
        title('Real Hypmogram')
        a=0;
%         figure
%         subplot(2,1,1)
%         tmp = Mylabels{itt};
%         Inds = 1:length(tmp);
%         title('Kmeans')    
%         tmp(tmp==0) = 2.5;
%         Flag = ~isempty(find(tmp==2.5,1));
%         H = plot(Inds(tmp==2.5),tmp(tmp==2.5),'*k');
%         hold on
%         Inds(tmp==2.5)=[];
%         tmp(tmp==2.5)=[];
%         DiscreteFlatPlot(Inds,tmp,Labelstr)
%         title('Clustered Hypmogram')
% 
%         subplot(2,1,2)
%         Inds = 1:length(Hyp);
%         DiscreteFlatPlot(Inds,Hyp',Labelstr)
%         xlabel('Time (s)');
%         title('Real Hypmogram')
    end
    %% Results
    Hyp(Hyp==6)=0;
    Labtemp = Mylabels{itt};
    Inds = and((Labtemp~=0),(Hyp~=0)');
    Unrecognized = (1-length(find(Inds))/length(Inds))*100;
%     disp("Unrecognized(%):")
%     disp(Unrecognized)
    cm = confusionmat(Hyp(Inds),Labtemp(Inds));
%     disp(cm)
    Accuracy = sum(diag(cm))/sum(sum(cm))*100;
%     disp("Accuracy all(%):")
%     disp(Accuracy)
    AccuracyStages = (diag(cm)./sum(cm,2)*100)';
%     disp(AccuracyStages)
    ACCsall(3:end,itt) = AccuracyStages;
    ACCsall(2,itt) = Accuracy;
    ACCsall(1,itt) = Unrecognized;
%     Kappa(itt) = kappa(cm);
    Kappa2(itt) = kappaindex(Labtemp(Inds),Hyp(Inds),5);
    F1score(itt) = meanF1(cm);
    
    
    Hyp(Hyp==1)=3;
    Hyp(Hyp==2)=3;
    Hyp(Hyp==3)=3;
    Labtemp(Labtemp==1)=3;
    Labtemp(Labtemp==2)=3;
    Labtemp(Labtemp==3)=3;
    cm = confusionmat(Hyp(Inds),Labtemp(Inds));
%     disp(cm)
    Accuracy = sum(diag(cm))/sum(sum(cm))*100;
%     disp("Accuracy 3(%):")
%     disp(Accuracy)
    AccuracyStages = (diag(cm)./sum(cm,2)*100)';
    ACCs3(2:end,itt) = AccuracyStages;
    ACCs3(1,itt) = Accuracy;
end
%% Plot final results
ACCs3_MDF = ACCs3';
ACCall_MDF = ACCsall';
save('ACCMDF.mat','ACCs3_MDF','ACCall_MDF');

LabelString ={'Overal', 'NREM1', 'NREM2', 'NREM3', 'REM','Wake'};

figure
x = (1:6);
y = mean(ACCsall(2:end,:),2,'omitnan')';
err = 1/sqrt(size(ACCsall,2))*std(ACCsall(2:end,:),0,2,'omitnan')';
errorbar(x,y,err,'o','Linewidth',1.2)
ylim([0,100])
xticks(1:6)
xticklabels(LabelString)
xlim([0,7])
ylabel("Average Accuracy (percent correct)")
disp(LabelString)
disp(x)
disp(y)

figure
for i=2:7
    subplot(6,1,i-1)
    plot(SubsetIds,ACCsall(i,:),'-.','Color',ColorCodes(i-1,:))
    hold on
    plot(SubsetIds,ACCsall(i,:),'*','Color',ColorCodes(i-1,:))
    ylabel(LabelString{i-1})
    ylim([0,100])
end
xlabel('Participants'' ID')

%---------- only 3 stages
LabelString ={'Overal', 'NREM', 'REM', 'Wake'};

figure
x = (1:4);
y = mean(ACCs3,2,'omitnan')';
err = 1/sqrt(size(ACCs3,2))*std(ACCs3,0,2,'omitnan')';
errorbar(x,y,err,'o','Linewidth',1.2)
ylim([0,100])
xticks(1:4)
xticklabels(LabelString)
xlim([0,5])
ylabel("Average Accuracy (percent correct)")
disp(LabelString)
disp(x)
disp(y)
