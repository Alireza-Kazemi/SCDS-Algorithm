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

%% Initialization 


Fnum = 62;
Subsetnum = 5;
Classnum=5;
EpochTime = 30;
Nfilt = 0;


 
load EEGFeatAll_150lag.mat
load FeatLabels.mat
Indsorig = [11,13,15,17,19]; randperm(39,Subsetnum);%31; Works awesome for 31st participant
Uinds = cell(2,Subsetnum);
UGroup = cell(2,Subsetnum);
Unumber = zeros(2,Subsetnum);
Mylabels = cell(2,Subsetnum);
CenterUs = cell(2,Subsetnum);
for itt = 1:Subsetnum
    load FeatIDXnew.mat
    Channel = 1;
    FeatIDX = FeatIDX08(Channel,:);
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
        tmp1 = kmeans(Features,2,'Start',Cents','MaxIter',1000);%
    %     [~,Ind] = min(mean(abs(Features - repmat(Cents(:,1)',size(Features,1),1)),2));
    %     Cents(:,1) = Features(Ind,:)';
    %     [~,Ind] = min(mean(abs(Features - repmat(Cents(:,2)',size(Features,1),1)),2));
    %     Cents(:,2) = Features(Ind,:)';
    %     tmp1 = kmedoids(Features,2,'Options',opts,'Replicates',20,'Algorithm','clara');

        tmp1(tmp1==2) = 0;
        tmp1(tmp1==1) = Classind;

        idx1(Classind,DataInd1) = tmp1;
        CenterUs{Channel,itt}(Classind,1:62) = mean(Feats(DataInd1(tmp1~=0),:),1);
        DataInd1(tmp1~=0)=[];
    end

    Thresh = round(.01*Xnum);
    Classind = 6;
    while(length(DataInd1)>Thresh)
        Features = Feats(DataInd1,:);
        opts = statset('MaxIter',1000);
        tmp1 = kmedoids(Features,2,'Options',opts,'Replicates',90,'Algorithm','pam');
        tmpind = tmp1;
        tmpind = hist(tmpind,1:2);
        [~,tmpind] = max(tmpind);
        tmp1(tmp1==tmpind) = 0;
        tmp1(tmp1~=0) = Classind;
        idx1(Classind,DataInd1) = tmp1;
        CenterUs{Channel,itt}(Classind,1:62) = mean(Feats(DataInd1(tmp1~=0),:),1);
        DataInd1(tmp1~=0)=[];
        Classind = Classind+1;
    end
    
    idx1(Classind,DataInd1) = Classind;
    CenterUs{Channel,itt}(Classind,1:62) = mean(Feats(DataInd1,:),1);
    
    Labtemp = sum(idx1);
    HypextFPz = Labtemp;
    Mylabels{Channel,itt} = HypextFPz;
    L1 = zeros(Classnum,Xnum*EpochTime);
    for i=1:Classnum
        idx1tmp = repmat(idx1(i,:)',1,EpochTime);
        if (Nfilt>0)
            L1(i,:) = round(medfilt1(reshape(idx1tmp',1,EpochTime*length(idx1tmp(:,1))),Nfilt));      
        else
            L1(i,:) = reshape(idx1tmp',1,EpochTime*length(idx1tmp(:,1)));
        end
    end

%     idx1 = sum(L1);


    %% Compare the performance
    Hypext = repmat(Hyp,1,EpochTime);
    Hypext = reshape(Hypext',1,EpochTime*length(Hyp));
    idx1 = repmat(Labtemp',1,EpochTime);
    idx1 = reshape(idx1',1,EpochTime*length(Labtemp));


    %% Results


    Labelstr = {'NREM1',...
                  'NREM2',...
                  'NREM3',...
                  'REM',...
                  'Wake'};
    j=1;
    for i=6:max(Labtemp)
        Labelstr{i} = char("U"+num2str(j));
        j=j+1;
    end

    figure
    subplot(2,1,1)
    Inds = 1:length(idx1);
    title('Kmeans')
    tmp = idx1;
    hold on
    DiscreteFlatPlot(Inds,tmp,Labelstr)
    title('Clustered Hypmogram based on FPz-Cz')
    subplot(2,1,2)
    DiscreteFlatPlot(Inds,Hypext,Labelstr)
    xlabel('Time (s)');
    title('Real Hypmogram')

    Uinds{Channel,itt} = find(Labtemp>5);
    UGroup{Channel,itt} = Hyp(Labtemp>5);
    Unumber(Channel,itt) = length(unique(UGroup{Channel,itt}));
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
    load FeatIDXnew.mat
    Channel = 2;
    FeatIDX = FeatIDX08(Channel,:);
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
    Featsall = Feats;
    Feats = Featsall(:,Fnum+1:end);
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
        tmp1 = kmeans(Features,2,'Start',Cents','MaxIter',1000);%
    %     [~,Ind] = min(mean(abs(Features - repmat(Cents(:,1)',size(Features,1),1)),2));
    %     Cents(:,1) = Features(Ind,:)';
    %     [~,Ind] = min(mean(abs(Features - repmat(Cents(:,2)',size(Features,1),1)),2));
    %     Cents(:,2) = Features(Ind,:)';
    %     tmp1 = kmedoids(Features,2,'Options',opts,'Replicates',20,'Algorithm','clara');

        tmp1(tmp1==2) = 0;
        tmp1(tmp1==1) = Classind;

        idx1(Classind,DataInd1) = tmp1;
        CenterUs{Channel,itt}(Classind,1:62) = mean(Feats(DataInd1(tmp1~=0),:),1);
        DataInd1(tmp1~=0)=[];

    end

    Thresh = round(.01*Xnum);
    Classind = 6;
    while(length(DataInd1)>Thresh)
        Features = Feats(DataInd1,:);
        opts = statset('MaxIter',1000);
        tmp1 = kmedoids(Features,2,'Options',opts,'Replicates',90,'Algorithm','pam');
        tmpind = tmp1;
        tmpind = hist(tmpind,1:2);
        [~,tmpind] = max(tmpind);
        tmp1(tmp1==tmpind) = 0;
        tmp1(tmp1~=0) = Classind;
        idx1(Classind,DataInd1) = tmp1;
        CenterUs{Channel,itt}(Classind,1:62) = mean(Feats(DataInd1(tmp1~=0),:),1);
        DataInd1(tmp1~=0)=[];
        Classind = Classind+1;
    end
    
    idx1(Classind,DataInd1) = Classind;
    CenterUs{Channel,itt}(Classind,1:62) = mean(Feats(DataInd1,:),1);
    
    Labtemp = sum(idx1);
    HypextFPz = Labtemp;
    Mylabels{Channel,itt} = HypextFPz;
    L1 = zeros(Classnum,Xnum*EpochTime);
    for i=1:Classnum
        idx1tmp = repmat(idx1(i,:)',1,EpochTime);
        if (Nfilt>0)
            L1(i,:) = round(medfilt1(reshape(idx1tmp',1,EpochTime*length(idx1tmp(:,1))),Nfilt));      
        else
            L1(i,:) = reshape(idx1tmp',1,EpochTime*length(idx1tmp(:,1)));
        end
    end

%     idx1 = sum(L1);


    %% Compare the performance
    Hypext = repmat(Hyp,1,EpochTime);
    Hypext = reshape(Hypext',1,EpochTime*length(Hyp));
    idx1 = repmat(Labtemp',1,EpochTime);
    idx1 = reshape(idx1',1,EpochTime*length(Labtemp));


    %% Results


    Labelstr = {'NREM1',...
                  'NREM2',...
                  'NREM3',...
                  'REM',...
                  'Wake'};
    j=1;
    for i=6:max(Labtemp)
        Labelstr{i} = char("U"+num2str(j));
        j=j+1;
    end

%     figure
%     subplot(2,1,1)
%     Inds = 1:length(idx1);
%     title('Kmeans')
%     tmp = idx1;
%     hold on
%     DiscreteFlatPlot(Inds,tmp,Labelstr)
%     title('Clustered Hypmogram based on FPz-Cz')
%     subplot(2,1,2)
%     DiscreteFlatPlot(Inds,Hypext,Labelstr)
%     xlabel('Time (s)');
%     title('Real Hypmogram')

    Uinds{Channel,itt} = find(Labtemp>5);
    UGroup{Channel,itt} = Hyp(Labtemp>5);
    Unumber(Channel,itt) = length(unique(UGroup{Channel,itt}));
end
close all



Freqs = zeros(3,6);
for ch = 1:2
    for i=1:size(UGroup,2)
        Hyp = Data.Hyp{i};
        Freqs(ch,:) = Freqs(ch,:)+hist(UGroup{ch,i},0:5);
        Freqs(3,1) = Freqs(3,1)+ length(Hyp(Hyp==0));
    end
end
Freqs(3,1) = Freqs(3,1)/2;
Freqs = Freqs/39;
Freqs(3,:)=[];
figure
bar(Freqs')
xticklabels({'Movement','NREM1',...
                  'NREM2',...
                  'NREM3',...
                  'REM',...
                  'Wake'})
xlabel('Expert''s labels')
ylabel('Average Number of Epochs')
legend({'FPz-Cz','Pz-Oz'})%,'Expert''s Label'})

%%
Freqs = zeros(2,5,size(UGroup,2));
for ch = 1:2
    for i=1:size(UGroup,2)
        Hyp = Data.Hyp{i};
        Freqs(ch,:,i) = hist(UGroup{ch,i},1:5);
    end
end

figure
x=1:5;
y = mean(Freqs,3,'omitnan');
err = std(Freqs,[],3,'omitnan')/5;
errorbar(x-.1,y(1,:),err(1,:),'o','Linewidth',1);
hold on 
errorbar(x+.1,y(2,:),err(2,:),'o','Linewidth',1);
xticks(1:5)
xticklabels({'NREM1',...
                  'NREM2',...
                  'NREM3',...
                  'REM',...
                  'Wake'})
xlabel('Expert''s labels')
ylabel('Average Number of Epochs')
ylim([0,52])
legend({'FPz-Cz','Pz-Oz'})%,'Expert''s Label'})

%% plots

figure
x = 1:2;
y = mean(Unumber','omitnan');
err = std(Unumber','omitnan');
errorbar(x,y,err,'o','Linewidth',3)
% boxplot(Unumber')
xlim([0,3])
xticks([1,2])
xticklabels({'FPz-Cz', 'Pz-Oz'})
ylabel("Number of new clusters")
% title("Number of Independent lusters")


figure;
Snum = 1;
chind = 1;
Fnum=1:11;
x=(Fnum);
y = mean(CenterUs{chind,Snum}(1:5,Fnum),'omitnan');
err = std(CenterUs{chind,Snum}(1:5,Fnum),'omitnan');
errorbar(x,y,err,'o','Linewidth',1)
hold on;
% x=(1:62)+.2;
% y = mean(CenterUs{2}(6:end,:),'omitnan');
% err = std(CenterUs{2}(6:end,:),'omitnan');
% errorbar(x,y,err,'o','Linewidth',1)
plot(CenterUs{chind,Snum}(6:end,Fnum)','.','MarkerSize',20)
Unum = size(CenterUs{chind,Snum}(6:end,Fnum),1);
Labelstr = cell(1,1+Unum);
Labelstr{1}='SleepStages';
for i=1:Unum
    Labelstr{i+1}=char("U"+num2str(i));
end
xticks(Fnum);
xlim([0,max(Fnum+1)])
xticklabels(CH_Lab(Fnum))
set(gca,'TicklabelInterpreter','None')
xtickangle(-45)
legend(Labelstr)

Msleep = zeros(19,62,2);
MUs = zeros(19,62,2);
for chind = 1:2
    for Snum=1:Subsetnum
        figure
        Fnum=1:62;
        x=(Fnum);
        y = mean(CenterUs{chind,Snum}(1:5,Fnum),'omitnan');
        Msleep(Snum,:,chind) = y;
        err = std(CenterUs{chind,Snum}(1:5,Fnum),'omitnan');
        errorbar(x,y,err,'o','Linewidth',1)
        hold on;
        plot(CenterUs{chind,Snum}(6:end,Fnum)','.','MarkerSize',20)
        MUs(Snum,:,chind) = mean(CenterUs{chind,Snum}(6:end,Fnum));
        Unum = size(CenterUs{chind,Snum}(6:end,Fnum),1);
        Labelstr = cell(1,1+Unum);
        Labelstr{1}='SleepStages';
        for i=1:Unum
            Labelstr{i+1}=char("U"+num2str(i));
        end
        xticks(Fnum);
        xlim([0,max(Fnum+1)])
        xticklabels(CH_Lab(Fnum))
        set(gca,'TicklabelInterpreter','None')
        xtickangle(-45)
        legend(Labelstr)
        grid on
        title("Participant "+num2str(Snum)+" Channel "+ num2str(chind) )
    end
end
clc
[H1,P,CI,STATS] = ttest(Msleep(:,:,1),MUs(:,:,1));
[H2,P,CI,STATS] = ttest(Msleep(:,:,2),MUs(:,:,2));
H = H1.*H2;
disp(CH_Lab(H==1)')
disp(CH_Lab(H1==1)')
disp(CH_Lab(H2==1)')
figure
plot(H1,'.','MarkerSize',20)
hold on
plot(H2,'.','MarkerSize',20)
xticks(Fnum);
xlim([0,max(Fnum+1)])
xticklabels(CH_Lab(Fnum))
set(gca,'TicklabelInterpreter','None')
xtickangle(-45)