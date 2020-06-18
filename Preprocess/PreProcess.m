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

%% Preprocessing 
clear;
clc;

load('Data.mat')
Fs = 100;
ChanNum = 2;
EpochLength = 30;

%% Bandpass Filtering
for i=1:length(Data.EEG)
    [Dat,d] = bandpass(Data.EEG{i}(1:2,:)',[.5,48],Fs);
    level = 5;
    wname = 'sym6';
    npc = 'heur';
    [~,~, npc] = wmspca(Dat,level,wname,npc);
    npc(1:3) = zeros(1,3);
    [x_sim, qual, npc] = wmspca(Dat,level,wname,npc);
    Data.EEG_MSPCA{i} = x_sim';
    Data.EEG{i} = Dat';
end

%% Epoching

for i=1:length(Data.EEG)
    Hyp = Data.Hyp{i};
    Hyp(Hyp==9) = [];
    Hyp(Hyp~=6)=1;
    Hyp(Hyp==6)=0;
    
    Ones = find(Hyp);
    StartIndex = (min(Ones) - 120-1)*30*Fs+1;
    EndIndex = (max(Ones) + 120)*30*Fs;
    L = EndIndex - StartIndex + 1;
    L = rem(L,EpochLength);
    
    Hyp = Data.Hyp{i};
    Data.Hyp{i} = Hyp((min(Ones)-120):(max(Ones)+120))';
    Data.EEG{i} = Data.EEG{i}(1:2,StartIndex:(EndIndex-L));
    Data.EEG_MSPCA{i} = Data.EEG_MSPCA{i}(1:2,StartIndex:(EndIndex-L));
    
    SampleLength = (EndIndex-L) - StartIndex+1;
    SegEEG = zeros(SampleLength/(EpochLength*Fs),EpochLength*Fs,ChanNum);
    for j=1:ChanNum
        SegEEG(:,:,j) = reshape(Data.EEG{i}(j,:),EpochLength*Fs,SampleLength/(EpochLength*Fs))'; 
    end
    Data.EEG{i} = SegEEG;
    
    SampleLength = (EndIndex-L) - StartIndex+1;
    SegEEG = zeros(SampleLength/(EpochLength*Fs),EpochLength*Fs,ChanNum);
    for j=1:ChanNum
        SegEEG(:,:,j) = reshape(Data.EEG_MSPCA{i}(j,:),EpochLength*Fs,SampleLength/(EpochLength*Fs))'; 
    end
    Data.EEG_MSPCA{i} = SegEEG;
    
end



%% Artifact Remove
WaveletLvl = 5;
Wavename = 'sym6';

for i=1:length(Data.EEG)
    for j=1:ChanNum
        NumEpoch = length(Data.EEG{i}(:,1,j));
        for k = 1:NumEpoch
            Temp = Data.EEG{i}(k,:,j);
            [C,L] = wavedec(Temp,WaveletLvl,Wavename);
            temp = wrcoef('a',C,L,Wavename,WaveletLvl);
            Data.EEG_Filt{i}(k,:,j) = Temp-temp;
            
            Temp = Data.EEG_MSPCA{i}(k,:,j);
            [C,L] = wavedec(Temp,WaveletLvl,Wavename);
            temp = wrcoef('a',C,L,Wavename,WaveletLvl);
            Data.EEG_EEGMSPCA_Filt{i}(k,:,j) = Temp-temp;
        end
    end
end

save("ReadyData.mat","Data","-v7.3");