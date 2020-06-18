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

clear;
clc;
close all;

load EEGFeatAll_150lag_Normalized.mat
load FeatLabels.mat

Cnum = 5;
FeatData = Data.EEG;
Fnum = 62;


ColorCodes = ['#0072BD';'#D95319';'#7E2F8E';'#77AC30';'#4DBEEE';'#A2142F'];
ColorCodeOrange = '#EDB120';

%% Choose a subset of participants and evaluate the best discriminability For Fpz-Cz
Subsetnum = 39;
ittnum = 1;
Centers = cell(1,Cnum);
DetailedCenters = cell(1,Cnum);
for Classind=1:Cnum
    disp(Classind);
    % Find the Centers
    Cents = zeros(Fnum,2,ittnum);
    for itt=1:ittnum
        % Select Subset of Participants
        Hyp = [];
        Feats = [];
        Inds = randperm(39,Subsetnum);
        for i=Inds
            Hyp = cat(1,Hyp,Data.Hyp{i}(1:end));
            Feats = cat(1,Feats,FeatData{i}(:,1:Fnum));
        end
        Hyp(Hyp==4)=3;
        Hyp(Hyp==5)=4;
        Hyp(Hyp==6)=5;
        Hyp(Hyp==0)=6;
        
        Labels = Hyp;
        X1 = Feats;
        X = X1(Labels==Classind,:);
        Cents(:,1,itt) = mean(X,1);
        X = X1(Labels~=Classind,:);
        Cents(:,2,itt) = mean(X,1);
    end
    Centers{Classind} = mean(Cents,3);
    DetailedCenters{Classind} = Cents;
end
%% plot samples

% MM = mean(Cents,3);
% figure;hold on;i=20;scatter(1:ittnum,Cents(i,1,:));
% scatter(1:ittnum,Cents(i,2,:));scatter(ittnum/2,MM(i,1),'filled');scatter(ittnum/2,MM(i,2),'filled')
% figure;hold on;i=10;scatter(1:ittnum,Cents(i,1,:));
% scatter(1:ittnum,Cents(i,2,:));scatter(ittnum/2,MM(i,1),'filled');scatter(ittnum/2,MM(i,2),'filled')

Labels = { 'NREM1',...
           'NREM2',...
           'NREM3',...
           'REM',...
           'Wake'};
figure
for i=1:Cnum
    subplot(3,2,i)
    MM = Centers{i};    
    scatter(1:length(MM(:,1)),MM(:,1),'filled','MarkerFaceColor',ColorCodes(i,:));
    hold on;
    scatter(1:length(MM(:,2)),MM(:,2),'MarkerEdgeColor',ColorCodeOrange);
    grid on;
%     grid minor;
    xlabel('Feature Index')
    ylabel('Estimated Centers')
    xlim([0,63])
%     ylim([0,1])
    title([Labels{i},' vs. other stages in Fpz-Cz'])
end
subplot(3,2,6)
hold on
for i=1:5
    scatter(0,1,'filled','MarkerFaceColor',ColorCodes(i,:))
end
xlim([10,11])
scatter(0,1,'MarkerEdgeColor',ColorCodeOrange)
legend(Labels{1},Labels{2},Labels{3},Labels{4},Labels{5},'Other Stages','Location','west')
axis off


%% Find significantly different centers by t-test
% T = zeros(Fnum,Cnum);
% for i=1:Cnum
%     Cents = DetailedCenters{i};
%     for Find = 1:Fnum
%         T(Find,i) = ttest(Cents(Find,1,:),Cents(Find,2,:),'alpha',.001);
%     end
% end
%% Sort Centers based on a Euclidean Distance
Diffs = zeros(Fnum,Cnum);
IDX = cell(1,Cnum);
for i=1:Cnum
    Cents = Centers{i};
    Diffs(:,i) = abs(Cents(:,1)-Cents(:,2));
    [~,IDX{i}] = sort(Diffs(:,i),'descend');
end

%% Find the Accuracy of each discriminative subspaces
% Fnum = 62;
% Hyp = [];
% Feats = [];
% PartID = 1:39;
% Inds = randperm(length(PartID),2);
% Inds = PartID(Inds);
% disp(Inds)
% for i=Inds
%     Hyp = cat(1,Hyp,Data.Hyp{i}(1:end));
%     Feats = cat(1,Feats,FeatData{i}(:,1:Fnum));
% end
% 
% Hyp(Hyp==4)=3;
% Hyp(Hyp==5)=4;
% Hyp(Hyp==6)=5;
% Hyp(Hyp==0)=6;
% Acc = zeros(Fnum,Cnum);
% 
% for i=1:Cnum
%     Labels = Hyp;
%     Labels(Labels ~= i) = -1;
%     Labels(Labels == i) = 1;
%     F_Inds = IDX{i};
%     X = Feats;
%     CNames = unique(Labels);
%     F = false(1,Fnum);
%     disp(i)
%     for j=1:length(F_Inds)
%         F(F_Inds(1:j))=true;
%         [~,Acc(j,i)] = train1channel([X,Labels],F,CNames);
%     end
% end
% 
% %% Figures
% Labels = { 'NREM1',...
%            'NREM2',...
%            'NREM3',...
%            'REM',...
%            'Wake'};
%    
%    
% figure
% Fnum = 62;
% for i=1:Cnum    
%     subplot(Cnum,1,i)
%     plot(Diffs(IDX{i},i),'LineWidth',2);
%     hold on
%     plot(Diffs(IDX{i},i),'*');
%     xticks(1:Fnum)
%     xticklabels(num2str(IDX{i}(1:Fnum)))
%     title([Labels{i},' VS other Stages in Fpz-Cz'])
%     ylabel('Difference')
%     xtickangle(90)
%     grid on
% end
% xlabel('Feature Index')
%        
%        
% FBest = zeros(1,Cnum);
% figure
% for i=1:Cnum
%     subplot(2,3,i)
%     plot(Acc(:,i),'LineWidth',2);
%     hold on
%     plot(Acc(:,i),'*');
%     [~,FBest(i)] = max(Acc(:,i));
%     
%     title([Labels{i},' VS other Stages in Fpz-Cz'])
%     ylabel('Accuracy')
%     xlabel('The most high ranked Features')
%     grid on
% end

%% Find Indexes for Fpz-Cz

FeatIDX05 = cell(1,Cnum);
FeatIDX03 = cell(1,Cnum);
FeatIDXleast  = cell(1,Cnum);
Inds = 1:Fnum;
for i=1:Cnum
    DD = Diffs(:,i);
    FeatIDX05{i}= Inds(DD>=.5);
    FeatIDX03{i}= Inds(DD>=.3);
    FeatIDXleast{i} = Inds(DD<.006);
end



for i=1:5
    disp(Labels{i});
    disp((CH_Lab(FeatIDX03{i})'));
end
for i=1:5
    disp(Labels{i});
    disp((CH_Lab(FeatIDX05{i})'));
end


Centers_FPz = Centers;
FeatIDX05_FPz = FeatIDX05;
FeatIDX03_FPz = FeatIDX03;
FeatIDXleast_FPz = FeatIDXleast;

%%                          Pz-Oz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose a subset of participants and evaluate the best discriminability For Pz-Oz
Subsetnum = 10;
ittnum = 300;
Centers = cell(1,Cnum);
DetailedCenters = cell(1,Cnum);
for Classind=1:Cnum
    disp(Classind);
    % Find the Centers
    Cents = zeros(Fnum,2,ittnum);
    for itt=1:ittnum
        % Select Subset of Participants
        Hyp = [];
        Feats = [];
        Inds = randperm(39,Subsetnum);
        for i=Inds
            Hyp = cat(1,Hyp,Data.Hyp{i}(1:end));
            Feats = cat(1,Feats,FeatData{i}(:,Fnum+1:end));
        end
        Hyp(Hyp==4)=3;
        Hyp(Hyp==5)=4;
        Hyp(Hyp==6)=5;
        Hyp(Hyp==0)=6;
        
        Labels = Hyp;
        X1 = Feats;
        X = X1(Labels==Classind,:);
        Cents(:,1,itt) = mean(X,1);
        X = X1(Labels~=Classind,:);
        Cents(:,2,itt) = mean(X,1);
    end
    Centers{Classind} = mean(Cents,3);
    DetailedCenters{Classind} = Cents;
end
%% plot samples

% MM = mean(Cents,3);
% figure;hold on;i=20;scatter(1:ittnum,Cents(i,1,:));
% scatter(1:ittnum,Cents(i,2,:));scatter(ittnum/2,MM(i,1),'filled');scatter(ittnum/2,MM(i,2),'filled')
% figure;hold on;i=10;scatter(1:ittnum,Cents(i,1,:));
% scatter(1:ittnum,Cents(i,2,:));scatter(ittnum/2,MM(i,1),'filled');scatter(ittnum/2,MM(i,2),'filled')

Labels = { 'NREM1',...
           'NREM2',...
           'NREM3',...
           'REM',...
           'Wake'};
figure
for i=1:Cnum
    subplot(3,2,i)
    MM = Centers{i};    
    scatter(1:length(MM(:,1)),MM(:,1),'filled','MarkerFaceColor',ColorCodes(i,:));
    hold on;
    scatter(1:length(MM(:,2)),MM(:,2),'MarkerEdgeColor',ColorCodeOrange);
    grid on;
%     grid minor;
    xlabel('Feature Index')
    ylabel('Estimated Centers')
    xlim([0,63])
%     ylim([0,1])
    title([Labels{i},' vs. other stages in Pz-Oz'])
end
subplot(3,2,6)
hold on
for i=1:5
    scatter(0,1,'filled','MarkerFaceColor',ColorCodes(i,:))
end
xlim([10,11])
scatter(0,1,'MarkerEdgeColor',ColorCodeOrange)
legend(Labels{1},Labels{2},Labels{3},Labels{4},Labels{5},'Other Stages','Location','west')
axis off


%% Find significantly different centers by t-test
% T = zeros(Fnum,Cnum);
% for i=1:Cnum
%     Cents = DetailedCenters{i};
%     for Find = 1:Fnum
%         T(Find,i) = ttest(Cents(Find,1,:),Cents(Find,2,:),'alpha',.001);
%     end
% end
%% Sort Centers based on a Euclidean Distance
Diffs = zeros(Fnum,Cnum);
IDX = cell(1,Cnum);
for i=1:Cnum
    Cents = Centers{i};
    Diffs(:,i) = abs(Cents(:,1)-Cents(:,2));
    [~,IDX{i}] = sort(Diffs(:,i),'descend');
end

%% Find the Accuracy of each discriminative subspaces
% Fnum = 62;
% Hyp = [];
% Feats = [];
% PartID = 1:39;
% Inds = randperm(length(PartID),2);
% Inds = PartID(Inds);
% disp(Inds)
% for i=Inds
%     Hyp = cat(1,Hyp,Data.Hyp{i}(1:end));
%     Feats = cat(1,Feats,FeatData{i}(:,1:Fnum));
% end
% 
% Hyp(Hyp==4)=3;
% Hyp(Hyp==5)=4;
% Hyp(Hyp==6)=5;
% Hyp(Hyp==0)=6;
% Acc = zeros(Fnum,Cnum);
% 
% for i=1:Cnum
%     Labels = Hyp;
%     Labels(Labels ~= i) = -1;
%     Labels(Labels == i) = 1;
%     F_Inds = IDX{i};
%     X = Feats;
%     CNames = unique(Labels);
%     F = false(1,Fnum);
%     disp(i)
%     for j=1:length(F_Inds)
%         F(F_Inds(1:j))=true;
%         [~,Acc(j,i)] = train1channel([X,Labels],F,CNames);
%     end
% end
% 
% %% Figures
% Labels = { 'NREM1',...
%            'NREM2',...
%            'NREM3',...
%            'REM',...
%            'Wake'};
%    
%    
% figure
% Fnum = 62;
% for i=1:Cnum    
%     subplot(Cnum,1,i)
%     plot(Diffs(IDX{i},i),'LineWidth',2);
%     hold on
%     plot(Diffs(IDX{i},i),'*');
%     xticks(1:Fnum)
%     xticklabels(num2str(IDX{i}(1:Fnum)))
%     title([Labels{i},' VS other Stages in Pz-Oz'])
%     ylabel('Difference')
%     xtickangle(90)
%     grid on
% end
% xlabel('Feature Index')
%        
%        
% FBest = zeros(1,Cnum);
% figure
% for i=1:Cnum
%     subplot(2,3,i)
%     plot(Acc(:,i),'LineWidth',2);
%     hold on
%     plot(Acc(:,i),'*');
%     [~,FBest(i)] = max(Acc(:,i));
%     
%     title([Labels{i},' VS other Stages in Pz-Oz'])
%     ylabel('Accuracy')
%     xlabel('The most high ranked Features')
%     grid on
% end

%% Find Indexes for Pz-Oz

FeatIDX05 = cell(1,Cnum);
FeatIDX03 = cell(1,Cnum);
FeatIDXleast  = cell(1,Cnum);
Inds = 1:Fnum;
for i=1:Cnum
    DD = Diffs(:,i);
    FeatIDX05{i}= Inds(DD>=.5);
    FeatIDX03{i}= Inds(DD>=.3);
    FeatIDXleast{i} = Inds(DD<.006);
end



for i=1:5
    disp(Labels{i});
    disp((CH_Lab(FeatIDX03{i})'));
end
for i=1:5
    disp(Labels{i});
    disp((CH_Lab(FeatIDX05{i})'));
end


Centers_Pz = Centers;
FeatIDX05_Pz = FeatIDX05;
FeatIDX03_Pz = FeatIDX03;
FeatIDXleast_Pz = FeatIDXleast;


%% Save
FeatIDX05 = cat(1,FeatIDX05_FPz,FeatIDX05_Pz);
FeatIDX03 = cat(1,FeatIDX03_FPz,FeatIDX03_Pz);
FeatIDXleast = cat(1,FeatIDXleast_FPz,FeatIDXleast_Pz);
Centers = cat(1,Centers_FPz,Centers_Pz);

x = 1; ACC(1,x) = Acc(length(FeatIDX03_Pz{x}),x);
x = 2; ACC(1,x) = Acc(length(FeatIDX03_Pz{x}),x);
x = 3; ACC(1,x) = Acc(length(FeatIDX03_Pz{x}),x);
x = 4; ACC(1,x) = Acc(length(FeatIDX03_Pz{x}),x);
x = 5; ACC(1,x) = Acc(length(FeatIDX03_Pz{x}),x);
x = 5; ACC(2,x) = Acc(length(FeatIDX03_FPz{x}),x);
x = 4; ACC(2,x) = Acc(length(FeatIDX03_FPz{x}),x);
x = 3; ACC(2,x) = Acc(length(FeatIDX03_FPz{x}),x);
x = 2; ACC(2,x) = Acc(length(FeatIDX03_FPz{x}),x);
x = 1; ACC(2,x) = Acc(length(FeatIDX03_FPz{x}),x);

save('FeatIDXnew.mat','FBest','IDX','FeatIDX05','FeatIDX03','Centers','FeatIDXleast')
