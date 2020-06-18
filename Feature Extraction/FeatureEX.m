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


%% Feature Extraction
clear;
clc;

%% Parameters
Fs = 100;
ChanNum = 2;
Phi = 1;
Rho = 0;
nLag = 150; 
Fnum=62;
%% Initialize
load EEG.mat % Load Data here; Data should be in the format of cell array on participants as {[Epochs, Samples, Channels]}

DataEEG = EEG;
FeatData = cell(1,length(DataEEG));

for Participant = 1:length(DataEEG)
    RawData = DataEEG{Participant};
    Xnum=length(RawData(:,1,1));
    Feats = [];
    disp(['Step:',num2str(Participant),' of ',num2str(length(DataEEG))]);
    for Ch=1:ChanNum        
        CH_Data = RawData(:,:,Ch);   
        CH_feats = zeros(Xnum,Fnum);
        n = length(CH_Data(1,:));  %number of samples in each epoch sample
        for ix=1:Xnum
            X = CH_Data(ix,:);
            %% Power related features: RSP SWI and HP  total number 28
            %================= RSP computation
            Ptot=bandpower(X);
            CH_feats(ix,1) = bandpower(X,Fs,[.5,2])/Ptot;
            CH_feats(ix,2) = bandpower(X,Fs,[2,4])/Ptot;
            CH_feats(ix,3) = bandpower(X,Fs,[4,6])/Ptot;
            CH_feats(ix,4) = bandpower(X,Fs,[6,8])/Ptot;
            CH_feats(ix,5) = bandpower(X,Fs,[8,10])/Ptot;
            CH_feats(ix,6) = bandpower(X,Fs,[10,12])/Ptot;
            CH_feats(ix,7) = bandpower(X,Fs,[12,14])/Ptot;
            CH_feats(ix,8) = bandpower(X,Fs,[14,16])/Ptot;
            CH_feats(ix,9) = bandpower(X,Fs,[16,24])/Ptot;
            CH_feats(ix,10)= bandpower(X,Fs,[24,32])/Ptot;
            CH_feats(ix,11)= bandpower(X,Fs,[32,45])/Ptot;

            %================= HP computation  Fc and S(Fc)
            [P,f] = pwelch(X,[],[],0:10^-1:Fs/2,Fs);
            
            %----> Fc    computation
            CH_feats(ix,12)=sum(P(6:40).*f(6:40))/sum(P(6:40));
            CH_feats(ix,13)=sum(P(41:80).*f(41:80))/sum(P(41:80));
            CH_feats(ix,14)=sum(P(81:120).*f(81:120))/sum(P(81:120));
            CH_feats(ix,15)=sum(P(121:160).*f(121:160))/sum(P(121:160));
            CH_feats(ix,16)=sum(P(161:351).*f(161:351))/sum(P(161:351));

            %----> Fsigma    computation
            CH_feats(ix,17)=sqrt(sum(P(6:40).*(f(6:40)-CH_feats(ix,11)).^2)/sum(P(6:40)));
            CH_feats(ix,18)=sqrt(sum(P(41:80).*(f(41:80)-CH_feats(ix,12)).^2)/sum(P(41:80)));
            CH_feats(ix,19)=sqrt(sum(P(81:120).*(f(81:120)-CH_feats(ix,13)).^2)/sum(P(81:120)));
            CH_feats(ix,20)=sqrt(sum(P(121:160).*(f(121:160)-CH_feats(ix,14)).^2)/sum(P(121:160)));
            CH_feats(ix,21)=sqrt(sum(P(161:351).*(f(161:351)-CH_feats(ix,13)).^2)/sum(P(161:351)));

            %----> S(fc) computation
            CH_feats(ix,22)=P(round(10*CH_feats(ix,11))+1);
            CH_feats(ix,23)=P(round(10*CH_feats(ix,12))+1);
            CH_feats(ix,24)=P(round(10*CH_feats(ix,13))+1);
            CH_feats(ix,25)=P(round(10*CH_feats(ix,14))+1);
            CH_feats(ix,26)=P(round(10*CH_feats(ix,15))+1);

            %================= SWI computation
            bspD = bandpower(X,Fs,[.6,4]);% subband spectrul power of delta
            bspD2 = bandpower(X,Fs,[2,4]);% subband spectrul power of delta
            bspT = bandpower(X,Fs,[4,8]);% subband spectrul power of theta
            bspA = bandpower(X,Fs,[8,11.5]);% subband spectrul power of aplha based on Jobert et. al. 1994
            %----> DSI    computation
            CH_feats(ix,27)=bspD/(bspT+bspA);
            %----> TSI    computation
            CH_feats(ix,28)=bspT/(bspD+bspA);
            %----> ASI    computation
            CH_feats(ix,29)=bspA/(bspD2+bspT);
            %% Hjorth features
            Xp=diff(X)*Fs;         %1st Differentiate
            Xpp=diff(X,2)*Fs^2;    %2nd Differentiate
            %----> Activity
            CH_feats(ix,30)=var(X);
            %----> Mobility
            CH_feats(ix,31)=sqrt(var(Xp)/var(X));
            %----> Complexity
            CH_feats(ix,32)=sqrt(var(Xpp)*var(X))/var(Xp);

            
            %% Skewness and Kortusis
            M2=0; %Moment 2
            M3=0; %Moment 3  
            M4=0; %Moment 4
            m=mean(X);
            for i=1:n
                M2=M2+(X(i)-m)^2;
                M3=M3+(X(i)-m)^3;
                M4=M4+(X(i)-m)^4;
            end
            M2=M2/n;
            M3=M3/n;
            M4=M4/n;
            %----> Skewness
            CH_feats(ix,33)=M3/sqrt(M2^3);
            %----> Kurtosis
            CH_feats(ix,34)=M4/M2^2;     

            %% Bispectrum Time Series
            [BTS,FFF] = btsestimate(X',Phi,Rho,nLag,Fs);
            CH_feats(ix,35:44) = abs(BTS([53,57,61,65,70,71,72,73,74,76]));
            CH_feats(ix,45:54) = angle(BTS([53,57,61,65,70,71,72,73,74,76]));
            
            %% Wavelet Coefficients
            WaveletLvl = 8;
            Wavename = 'sym6';
            type = 'd';
            [C,L] = wavedec(X,WaveletLvl,Wavename);
            PwaveletCoefs = zeros(1,WaveletLvl);
            for i=1:WaveletLvl
                PwaveletCoefs(i) = bandpower(wrcoef('d',C,L,Wavename,i))/Ptot;
            end
            CH_feats(ix,55:62)=PwaveletCoefs;
        end  
        Feats = cat(2,Feats,CH_feats);
    end
%     Feats = (Feats - repmat(min(Feats),Xnum,1))./(repmat(max(Feats),Xnum,1)-repmat(min(Feats),Xnum,1));
    FeatData{Participant} = Feats;
end

Data.EEG = FeatData;

%% Save
save('EEGFeatAll_50lag.mat','Data','-v7.3');
%%