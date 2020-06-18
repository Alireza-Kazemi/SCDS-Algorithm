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

function [BTS,Waxis] = btsestimate(Sig,Phi,Rho,nLag,Fc)

L=length(Sig(1,:));
Np = 2*nLag+1;
Sig = double(Sig);
BTS = cell(L,1);


for i=1:L
    x = Sig(:,i);
    [B,A] = butter(5,[1,45]/(Fc*.5));
    x = filter(B,A,x);
    [Bispec,Waxis]= bispectrum(x,nLag,Fc);
    W1 = Waxis;
    W2 = Phi*Waxis+Rho;
    W1(abs(W2)>max(Waxis))=[];
    W2(abs(W2)>max(Waxis))=[];
    W1=round(W1*Np/Fc+nLag+1);
    W2=round(W2*Np/Fc+nLag+1);
    BTS{i} = Bispec(W1+(W2-1)*Np);
end
x = Sig(:,1);
[~,Waxis]= bispectrum(x,nLag,Fc);
BTS = cell2mat(BTS);

end