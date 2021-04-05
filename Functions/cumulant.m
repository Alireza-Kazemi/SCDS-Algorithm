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

function C_X = cumulant(X,T)

% Cumulant(X,T) is a function to calculate 3rd Cumulant which is
% C(t1,t2) = E[ X(k)*X(k+t1)*X(k+t2) ]
% X is vector contain the series and T contains the lags values

N=length(X);
X = X-mean(X);   
    
t1=T(1);
t2=T(2);
L = 2*t1+1;
Xr = repmat(X,1,t1+1);
Xdf = zeros(size(Xr));
% Xdb = zeros(size(Xr));
C_X = zeros(L,L);
for i=t1:-1:0
    Xdf(1:N-i,i+1)=X(i+1:N);  
    
    ind1 = 1;
    ind2 = i+1;
    for j=(t2+i+1):L
        C_X(t1+i+1,j)= (Xr(:,ind1).*Xdf(:,i+1))'*...
                       Xdf(:,ind2);
         ind1 = ind1+1;
         ind2 = ind2+1;
    end                  
                  
    m=i:t2;        
    n=ones(1,t2-i+1)*i;
    C_X((t2+1-m)+(t1+n-m)*L)= C_X(t1+i+1,t2+i+1:L);
    C_X((t2+1-n)+(t1+m-n)*L)= C_X(t1+i+1,t2+i+1:L);
end  
C_X = 1/N*C_X;

C_X = C_X + triu(C_X,1)';

end

