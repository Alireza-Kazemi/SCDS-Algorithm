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

%% RemoveOutliers
clear;
clc;
load EEGFeatAll_150lag_June2020.mat
load FeatLabels.mat
%% Normalize Features and Remove Outliers (>3sd)
FeatData = Data.EEG;
HypData = Data.Hyp;
for SID = 1:length(FeatData)
    Dat = FeatData{SID};
    H_Dat = HypData{SID};
    H_Dat(H_Dat==4)=3;
    H_Dat(H_Dat==5)=4;
    H_Dat(H_Dat==6)=5;
    H_Dat(H_Dat==0)=6;
    DatNormalized = zscore(Dat);
    Index = DatNormalized;
    Index(abs(DatNormalized)<3)=0;
    Index(abs(DatNormalized)>3)=1;
    for F_idx= 1:size(Index,2)
        Idx = logical(Index(:,F_idx));
        SleepStage = unique(H_Dat(Idx));
        Hyp_idx = H_Dat;
        Hyp_idx(~Idx)=-1;
        for S_idx = SleepStage'
            if(~isempty(DatNormalized(H_Dat == S_idx & Hyp_idx~=S_idx, F_idx)))
                MeanVal = mean(DatNormalized(H_Dat == S_idx & Hyp_idx~=S_idx, F_idx));
                DatNormalized(Hyp_idx==S_idx, F_idx) = MeanVal;
            end
        end
    end
    if(~isempty(DatNormalized(isnan(DatNormalized))))
        disp(SID)
        Inds = 1:124;
        Inds = Inds(isnan(mean(DatNormalized)));
        Inds = Inds(Inds<63);
        disp(CH_Lab(Inds))
    end
    Feats = DatNormalized;
    Xnum = size(Feats,1);
    Feats = (Feats - repmat(min(Feats),Xnum,1))./(repmat(max(Feats),Xnum,1)-repmat(min(Feats),Xnum,1));
    FeatData{SID} = Feats;
end
Data.EEG = FeatData;
%% Normalize Features and Remove Outliers (>3sd)

FeatData = Data.EEG_BandNoised;
HypData = Data.Hyp;
for SID = 1:length(FeatData)
    Dat = FeatData{SID};
    H_Dat = HypData{SID};
    H_Dat(H_Dat==4)=3;
    H_Dat(H_Dat==5)=4;
    H_Dat(H_Dat==6)=5;
    H_Dat(H_Dat==0)=6;
    DatNormalized = zscore(Dat);
    Index = DatNormalized;
    Index(abs(DatNormalized)<3)=0;
    Index(abs(DatNormalized)>3)=1;
    for F_idx= 1:size(Index,2)
        Idx = logical(Index(:,F_idx));
        SleepStage = unique(H_Dat(Idx));
        Hyp_idx = H_Dat;
        Hyp_idx(~Idx)=-1;
        for S_idx = SleepStage'
            if(~isempty(DatNormalized(H_Dat == S_idx & Hyp_idx~=S_idx, F_idx)))
                MeanVal = mean(DatNormalized(H_Dat == S_idx & Hyp_idx~=S_idx, F_idx));
                DatNormalized(Hyp_idx==S_idx, F_idx) = MeanVal;
            end
        end
    end
    if(~isempty(DatNormalized(isnan(DatNormalized))))
        disp(SID)
        Inds = 1:124;
        Inds = Inds(isnan(mean(DatNormalized)));
        Inds = Inds(Inds<63);
        disp(CH_Lab(Inds))
    end
    Feats = DatNormalized;
    Xnum = size(Feats,1);
    Feats = (Feats - repmat(min(Feats),Xnum,1))./(repmat(max(Feats),Xnum,1)-repmat(min(Feats),Xnum,1));
    FeatData{SID} = Feats;
end
Data.EEG_BandNoised = FeatData;

%% Normalize Features and Remove Outliers (>3sd)

FeatData = Data.EEG_FullNoised;
HypData = Data.Hyp;
for SID = 1:length(FeatData)
    Dat = FeatData{SID};
    H_Dat = HypData{SID};
    H_Dat(H_Dat==4)=3;
    H_Dat(H_Dat==5)=4;
    H_Dat(H_Dat==6)=5;
    H_Dat(H_Dat==0)=6;
    DatNormalized = zscore(Dat);
    Index = DatNormalized;
    Index(abs(DatNormalized)<3)=0;
    Index(abs(DatNormalized)>3)=1;
    for F_idx= 1:size(Index,2)
        Idx = logical(Index(:,F_idx));
        SleepStage = unique(H_Dat(Idx));
        Hyp_idx = H_Dat;
        Hyp_idx(~Idx)=-1;
        for S_idx = SleepStage'
            if(~isempty(DatNormalized(H_Dat == S_idx & Hyp_idx~=S_idx, F_idx)))
                MeanVal = mean(DatNormalized(H_Dat == S_idx & Hyp_idx~=S_idx, F_idx));
                DatNormalized(Hyp_idx==S_idx, F_idx) = MeanVal;
            end
        end
    end
    if(~isempty(DatNormalized(isnan(DatNormalized))))
        disp(SID)
        Inds = 1:124;
        Inds = Inds(isnan(mean(DatNormalized)));
        Inds = Inds(Inds<63);
        disp(CH_Lab(Inds))
    end
    Feats = DatNormalized;
    Xnum = size(Feats,1);
    Feats = (Feats - repmat(min(Feats),Xnum,1))./(repmat(max(Feats),Xnum,1)-repmat(min(Feats),Xnum,1));
    FeatData{SID} = Feats;
end
Data.EEG_FullNoised = FeatData;

%% Save Data
save('EEGFeatAll_150lag_MNNormalized.mat','Data','-v7.3');