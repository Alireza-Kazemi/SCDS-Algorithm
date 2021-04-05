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


function [Bspec,waxis] = bispectrum(x,nLag,fc,nFFT)
%% Initialize
Np = 2*nLag+1;
s=size(x);
if(s(1)==1)
    x=x';
end
x=double(x);
Cum = cumulant(x,[nLag,nLag]); % Cumulants
%% Lag Window 
Window = parzenwin(Np);
% 
% indx = (-nLag:nLag)';
% Window = sinc(indx/nLag);
% Window(nLag+1)=1;
%          
BWind = zeros(Np,Np);

for i=nLag:-1:0
    ind=i:Np-1;
    BWind(ind*Np+ind+1-i) = Window(nLag-i+1);
end
BWind = BWind + triu(BWind,1)';
Window = Window*Window';
BWind = BWind.*Window;
% figure;surf(abs(BWind));
% figure;surf(abs(fftshift(fft2(BWind))));
WCum = Cum.*BWind;
%% Bispectrum 
nfft = nFFT;
% Bspec = fft2(WCum, nfft, nfft); 
Bspec = fft2(WCum, nfft, nfft); 
Bspec = fftshift(Bspec); 

if (rem(nfft,2) == 0) 
    waxis = (-nfft/2:(nfft/2-1))*fc/nfft; 
else
    waxis = (-(nfft-1)/2:(nfft-1)/2)*fc/nfft; 
end
% figure
% contour(waxis,waxis,abs(Bspec),4), grid on 
% title('Bispectrum estimated via the indirect method')
% xlabel('f1'), ylabel('f2') 
% set(gcf,'Name','Hosa BISPECI')


end
