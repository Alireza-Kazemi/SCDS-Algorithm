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
%% Pull Pre-knowledge Data
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
%% Find optimum order
% % Select All Participants
% close all
% Hyp = [];
% Feats = [];
% for i=1:length(FeatData)
%     Hyp = cat(1,Hyp,Data.Hyp{i});
%     Feats = cat(1,Feats,FeatData{i});
% end
% Hyp(Hyp==4)=3;
% Hyp(Hyp==5)=4;
% Hyp(Hyp==6)=5;
% Hyp(Hyp==0)=6;
% 
% excludes = [5,4];
% indexes = ~ismember(Hyp,excludes);
% Hyp = Hyp(indexes);
% Feats = Feats(indexes,:);
% R = zeros(Cnum,Fnum*Ch_num);
% for Classind=[1,2,3]
%     disp(Classind);
%     Labels = Hyp;
%     pre = circshift(1:5,1);
%     post = circshift(1:5,-1);
%     seqlabs = [pre(Classind),Classind,post(Classind)];
%     Labels(~ismember(Labels,seqlabs)) = 0;
%     Labels(Labels == post(Classind)) = 0.4;
%     Labels(Labels == pre(Classind)) = 0.4;
%     Labels(Labels == Classind) = 1;
% %     Labels(Labels ~= Classind) = 0;
% %     Labels(Labels == Classind) = 1;
%     for F_idx=1:Fnum*Ch_num
%         F_values = Feats(:,F_idx);
%         temp = corrcoef(F_values,Labels);
%         R(Classind,F_idx) = abs(temp(2,1));
%     end
%     [~,bestF] = sort(R(Classind,:),'descend');
%     figure;
%     Inds = randperm(length(Labels),1000);
%     Inds = sort(Inds);
%     for i=1:5
%         subplot (5,1,i)
%         plot(Labels(Inds),'o')
%         hold on
%         plot(Feats(Inds,bestF(i)))
%         ylabel(CH_Lab_total{bestF(i)},'Interpreter', 'none')
%         title([SleepStages{Classind},'  R_Both = ',num2str(R(Classind,bestF(i)))],'Interpreter', 'none')
%     end
% end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Flag_AllSubjects == 1)
%% Find correaltion of each feature with different sleep stages Both Channels
% All participants
Feats = FeatDat;
Hyp = HypDat;

RBothChannels = zeros(Cnum,FnumTot);
for Classind=Optimum_order
    Labels = Hyp;
    if(Flag_GradedCorrelation == 1)
        pre = circshift(1:5,1);
        post = circshift(1:5,-1);
        seqlabs = [pre(Classind),Classind,post(Classind)];
        Labels(~ismember(Labels,seqlabs)) = 0;
        Labels(Labels == post(Classind)) = GradedCodingVal;
        Labels(Labels == pre(Classind)) = GradedCodingVal;
        Labels(Labels == Classind) = 1;
    else
        Labels(Labels ~= Classind) = 0;
        Labels(Labels == Classind) = 1;
    end

    for F_idx=1:FnumTot
        F_values = Feats(:,F_idx);
        temp = corrcoef(F_values,Labels);
        RBothChannels(Classind,F_idx) = abs(temp(2,1));
    end
    [~,bestF] = sort(RBothChannels(Classind,:),'descend');
    if(PlotCorrelationFlag == 1)
        figure;
        Inds = randperm(length(Labels),5000);
        Inds = sort(Inds);
        for i=1:5
            subplot (5,1,i)
            plot(Labels(Inds),'o')
            hold on
            plot(Feats(Inds,bestF(i)))
            ylabel(CH_Lab_total{bestF(i)},'Interpreter', 'none')
            title([SleepStages{Classind},'  R_Both = ',num2str(RBothChannels(Classind,bestF(i)))],'Interpreter', 'none')
        end
    end
    if(Flag_RemoveStageWise == 1)
        exInds = Hyp==Classind;
        Hyp(exInds) = [];
        Feats(exInds,:) = [];
    end
end

%% Find Centers
% All participants
Feats = FeatDat;
Hyp = HypDat;

Centers = cell(1,Cnum);
for Classind=Optimum_order
    Cents = zeros(FnumTot,2);
    Labels = Hyp;
    X1 = Feats;
    X = X1(Labels==Classind,:);
    Cents(:,1) = mean(X,1);
    X = X1(Labels~=Classind,:);
    Cents(:,2) = mean(X,1);
    Centers{Classind} = Cents;
    if(Flag_RemoveStageWise_Center == 1)
        exInds = Hyp==Classind;
        Hyp(exInds) = [];
        Feats(exInds,:) = [];
    end
end

%% Sort Features based on their correlation
% All participants
RhoValues = zeros(FnumTot,Cnum);
IDX = cell(1,Cnum);
for i=1:Cnum
    Rs = RBothChannels(i,:);
    RhoValues(:,i) = Rs;
    [tmp,IDX{i}] = sort(Rs,'descend');
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
%% Find correaltion of each feature with different sleep stages Both Channels
% individual Subjects
% Feats = FeatDat;
% Hyp = HypDat;
RBothChannels = zeros(Cnum,FnumTot,length(SubsetIds));
for SubjID = 1:length(SubsetIds)
    
    Feats = FeatData{SubsetIds(SubjID)};
    Hyp = Data.Hyp{SubsetIds(SubjID)};
    Hyp(Hyp==4)=3;
    Hyp(Hyp==5)=4;
    Hyp(Hyp==6)=5;
    Hyp(Hyp==0)=6;
    
    for Classind=1:Cnum
        Labels = Hyp;
        if(Flag_GradedCorrelation == 1)
            pre = circshift(1:5,1);
            post = circshift(1:5,-1);
            seqlabs = [pre(Classind),Classind,post(Classind)];
            Labels(~ismember(Labels,seqlabs)) = 0;
            Labels(Labels == post(Classind)) = GradedCodingVal;
            Labels(Labels == pre(Classind)) = GradedCodingVal;
            Labels(Labels == Classind) = 1;
        else
            Labels(Labels ~= Classind) = 0;
            Labels(Labels == Classind) = 1;
        end
        for F_idx=1:FnumTot
            F_values = Feats(:,F_idx);
            temp = corrcoef(F_values,Labels);
            RBothChannels(Classind,F_idx,SubjID) = abs(temp(2,1));
        end
        if(Flag_RemoveStageWise == 1)
            exInds = Hyp==Classind;
            Hyp(exInds) = [];
            Feats(exInds,:) = [];
        end
    end
end

%% Find Centers
% individual Subjects
Centers = cell(1,Cnum);
Centers_SE = cell(1,Cnum);
Centers_Ind = zeros(FnumTot,2,Cnum,length(SubsetIds));

for SubjID = 1:length(SubsetIds)
    
    Feats = FeatData{SubsetIds(SubjID)};
    Hyp = Data.Hyp{SubsetIds(SubjID)};
    Hyp(Hyp==4)=3;
    Hyp(Hyp==5)=4;
    Hyp(Hyp==6)=5;
    Hyp(Hyp==0)=6;
    
    for Classind=Optimum_order
        Cents = zeros(FnumTot,2);
        Labels = Hyp;
        X1 = Feats;
        X = X1(Labels==Classind,:);
        Cents(:,1) = mean(X,1);
        X = X1(Labels~=Classind,:);
        Cents(:,2) = mean(X,1);
        Centers_Ind(:,:,Classind,SubjID) = Cents;
        if(Flag_RemoveStageWise_Center == 1)
            exInds = Hyp==Classind;
            Hyp(exInds) = [];
            Feats(exInds,:) = [];
        end
    end
    
end

for i=1:Cnum
    tmp = squeeze(Centers_Ind(:,:,i,:));
    Centers{i} = mean(tmp,3,'omitnan');
    Centers_SE{i} = std(tmp,0,3,'omitnan')/sqrt(size(tmp,3));
end
%% Sort Features based on their correlation
% individual Subjects
RhoValues = zeros(FnumTot,Cnum);
RhoValues_SE = zeros(FnumTot,Cnum);
IDX = cell(1,Cnum);
for i=1:Cnum
    tmp = squeeze(RBothChannels(i,:,:));
    if(size(tmp,1)==1)
        tmp=tmp';
    end
    Rs = mean(tmp,2,'omitnan');
    RhoValues(:,i) = Rs;
    RhoValues_SE(:,i) =std(tmp,0,2,'omitnan')/sqrt(size(tmp,2));
    [tmp,IDX{i}] = sort(Rs,'descend');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%% Find Indexes for Both Channels

FeatIDXMDF = cell(1,Cnum);
FeatIDXLDF  = cell(1,Cnum);
Inds = 1:FnumTot;
CFFInds = 1:11;
for i=2:Ch_num
    CFFInds = [CFFInds,CFFInds+Fnum];
end
for i=1:Cnum
    Rs = RhoValues(:,i);
    FeatIDXMDF{i} = IDX{i}(1:MDF_Fnum(i))';%Inds(Rs>(max(Rs)-1*std(Rs)) | Rs>=0.45);%
    FeatIDXLDF{i} = IDX{i}(end-9:end)';
    FeatIDXCFF{i} = CFFInds;
end

%% Find the Accuracy of each discriminative subspaces
if(PlotFlag ==1)
    Hyp = [];
    Feats = [];
    PartID = 1:39;
    Inds = randperm(length(PartID),2);
    Inds = PartID(Inds);
    disp(Inds)
    for i=Inds
        Hyp = cat(1,Hyp,Data.Hyp{i}(1:end));
        Feats = cat(1,Feats,FeatData{i}(:,1:FnumTot));
    end

    Acc = zeros(FnumTot,Cnum);

    for i=1:Cnum
        Labels = Hyp;
        Labels(Labels ~= i) = -1;
        Labels(Labels == i) = 1;
        F_Inds = IDX{i};
        X = Feats;
        CNames = unique(Labels);
        F = false(1,FnumTot);
        disp(i)
        for j=1:length(F_Inds)
            F(F_Inds(1:j))=true;
            [~,Acc(j,i)] = train1channel124([X,Labels],F,CNames);
        end
    end

    % Figures

    Labels = { 'NREM1',...
               'NREM2',...
               'NREM3',...
               'REM',...
               'Wake'};


    figure
    for i=1:Cnum    
        subplot(Cnum,1,i)
        plot(RhoValues(IDX{i},i),'LineWidth',2);
        hold on
        plot(RhoValues(IDX{i},i),'*');
        xticks(1:FnumTot)
        xticklabels(string(IDX{i}(1:FnumTot)))
        title([Labels{i},' VS other Stages in BothChannels'])
        ylabel('Correlations')
        xtickangle(90)
        grid on
    end
    xlabel('Feature Index')


    FBest = zeros(1,Cnum);
    figure
    for i=1:Cnum
        subplot(2,3,i)
        plot(Acc(:,i),'LineWidth',2);
        hold on
        plot(Acc(:,i),'*');
        [~,FBest(i)] = max(Acc(:,i));

        title([Labels{i},' VS other Stages in BothChannels'])
        ylabel('Accuracy')
        xlabel('The most high ranked Features')
        grid on
    end
    
    for i=1:5
        disp(Labels{i});
        disp((CH_Lab_total(FeatIDXMDF{i})'));
    end
end





%% Save

save('FeatIDXBoth.mat','IDX','FeatIDXMDF','FeatIDXLDF','FeatIDXCFF','Centers','Optimum_order')

%% plot Differences
if(Flag_AllSubjects ==1)
% All participants
if (plotDifferencesFlag ==1)
    Labels = { 'NREM1',...
               'NREM2',...
               'NREM3',...
               'REM',...
               'Wake'};
           
    figure;
    RhoValues = RBothChannels';
    subplot 211
    hold on
    FnumTot = 62;
    Shapes = {'^','s','d','p','h'};
    for i=1:Cnum
        plot((1:FnumTot)',RhoValues(1:FnumTot,i),['-',Shapes{i}]);
    end
    % grid on;
    xlabel('Feature index')
    ylabel('Correlation coefficients on FPz-Cz')
    xlim([0,63])
    xticks(1:FnumTot)
    xtickangle(90)
    legend(Labels,'Location','west')
    box on
    
    subplot 212
    hold on
    Shapes = {'^','s','d','p','h'};
    for i=1:Cnum
        plot((1:FnumTot)',RhoValues((1+FnumTot):end,i),['-',Shapes{i}]);
    end
    % grid on;
    xlabel('Feature index')
    ylabel('Correlation coefficients on Pz-Oz')
    xlim([0,63])
    xticks(1:FnumTot)
    xtickangle(90)
    legend(Labels,'Location','west')
    box on
end
%% plot Differences
else
% individual Subjects
if (plotDifferencesFlag ==1)
    Labels = { 'NREM1',...
               'NREM2',...
               'NREM3',...
               'REM',...
               'Wake'};

    figure;
    subplot 211
    hold on
    Inds = (1:Fnum);
    for i=1:Cnum
        x = (1:Fnum);
        y = RhoValues(Inds,i);
        err = RhoValues_SE(Inds,i);
        errorbar(x,y,err,'o','Linewidth',1.2,'Color',ColorCodes(i+1,:))
    end
    xlabel('Feature index')
    ylabel('Correlation coefficients on FPz-Cz')
    xlim([0,Fnum+1])
    xticks(1:FnumTot)
    xtickangle(90)
    legend(Labels,'Location','west')
    box on

    subplot 212
    hold on
    Inds = Fnum+1:FnumTot;
    for i=1:Cnum
        x = (1:Fnum);
        y = RhoValues(Inds,i);
        err = RhoValues_SE(Inds,i);
        errorbar(x,y,err,'o','Linewidth',1.2,'Color',ColorCodes(i+1,:))
    end
    xlabel('Feature index')
    ylabel('Correlation coefficients on Pz-Oz')
    xlim([0,Fnum+1])
    xticks(1:Fnum)
    xtickangle(90)
    legend(Labels,'Location','west')
    box on
end
end