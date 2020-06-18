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

% Cumulant(X,T) is a function to calculate 2nd and 3rd Cumulant which is
% C(t) = E[ X(k)*X(k+t) ]  and C(t1,t2) = E[ X(k)*X(k+t1)*X(k+t2) ]
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
    C_X(t1+i+1,t2+i+1:L)= 1/N*sum(Xr(:,1:t2-i+1).*...
                          repmat(Xdf(:,i+1),1,t2+1-i).*...
                          Xdf(:,i+1:t2+1) );                               
    m=i:t2;        
    n=ones(1,t2-i+1)*i;
    C_X((t2+1-m)+(t1+n-m)*L)= C_X(t1+i+1,t2+i+1:L);
    C_X((t2+1-n)+(t1+m-n)*L)= C_X(t1+i+1,t2+i+1:L);
end  

% for i=-1:-1:-1*t1
%     Xdb(-i+1:N,-i+1)=X(1:N+i);
%     C_X(t1+i+1,L+i+1:L)= 1/N*sum(Xr(:,1:-i).*...
%                           repmat(Xdb(:,-i+1),1,-i).*...
%                           Xdf(:,t2+i+2:t2+1) );                             
% end 
C_X = C_X + triu(C_X,1)';

end

% 
%%
% function C_X = cumulant(X,T)
% 
% % Cumulant(X,T) is a function to calculate 2nd and 3rd Cumulant which is
% % C(t) = E[ X(k)*X(k+t) ]  and C(t1,t2) = E[ X(k)*X(k+t1)*X(k+t2) ]
% % X is vector contain the series and T contains the lags values
% 
% N=length(X);
% X = X-mean(X);
% if (length(T)==1)   
%     
%     C_X = zeros(1,2*T+1);
%     for i=-T:T
%         if (i>=0)
%             C_X(T+i+1)= 1/N*sum((X-mean(X)).*([(X(i+1:end)-mean(X));zeros(i,1)]));
%         elseif (i<0)
%             C_X(T+i+1)= 1/N*sum((X-mean(X)).*([zeros(abs(i),1);(X(1:end+i)-mean(X))]));        
%         end
%     end    
% elseif (length(T)==2)        
%     tic
%     t1=T(1);
%     t2=T(2);
%     C_X = zeros(2*t1+1,2*t2+1);
%     for i=0:t1
%         for j=i+1:t2
%             C_X(t1+i+1,t2+j+1)= 1/N*sum(X.*...
%                                 ([X(i+1:end);zeros(i,1)]).*...
%                                 ([X(j+1:end);zeros(j,1)]));
% %             C_X(t2+1+j,t1+1+i)= C_X(t1+i+1,t2+j+1);
%             C_X(t2+1-j,t1+1+i-j)= C_X(t1+i+1,t2+j+1);
% %             C_X(t2+1+i-j,t1+1-j)= C_X(t1+i+1,t2+j+1);
%             C_X(t2+1-i,t1+1+j-i)= C_X(t1+i+1,t2+j+1);
% %             C_X(t2+1+j-i,t1+1-i)= C_X(t1+i+1,t2+j+1);
%         end
%     end 
%     for i=0:t1
%         j=i;
%         C_X(t1+i+1,t2+j+1)= 1/N*sum(X.*...
%                             ([X(i+1:end);zeros(i,1)]).*...
%                             ([X(j+1:end);zeros(j,1)]));                
%         C_X(t2+1-i,t1+1+j-i)= C_X(t1+i+1,t2+j+1);
% %         C_X(t2+1+j-i,t1+1-i)= C_X(t1+i+1,t2+j+1);
%     end 
%     
%     for i=-1*t1:-1
%         for j=t1+1+i:t2
%             C_X(t1+i+1,t2+j+1)= 1/N*sum(X.*...
%                                     ([zeros(abs(i),1);X(1:end+i)]).*...
%                                     ([X(j+1:end);zeros(j,1)]));
% %             C_X(t2+1+j,t1+1+i)= C_X(t1+i+1,t2+j+1);                       
%         end
%     end 
%     C_X = C_X + triu(C_X,1)';
%     toc
% else
%     disp('Error! check Order or T');
%     disp('Order=2 -> T=t       &&      T<length(X)');
%     disp('Order=3 -> T=[t1,t2] && max(T)<length(X)');
% end
% 
% end


%% cumulant matrix
% function C_X = cumulant(X,T)
% 
% % Cumulant(X,T) is a function to calculate 2nd and 3rd Cumulant which is
% % C(t) = E[ X(k)*X(k+t) ]  and C(t1,t2) = E[ X(k)*X(k+t1)*X(k+t2) ]
% % X is vector contain the series and T contains the lags values
% 
% N=length(X);
% if (length(T)==1)   
%     
%     C_X = zeros(1,2*T+1);
%     for i=-T:T
%         if (i>=0)
%             C_X(T+i+1)= 1/N*sum((X-mean(X)).*([(X(i+1:end)-mean(X));zeros(i,1)]));
%         elseif (i<0)
%             C_X(T+i+1)= 1/N*sum((X-mean(X)).*([zeros(abs(i),1);(X(1:end+i)-mean(X))]));        
%         end
%     end    
% elseif (length(T)==2)        
%     
%     t1=T(1);
%     t2=T(2);
%     C_X = zeros(2*t1+1,2*t2+1);
%     for i=-t1:t1
%         for j=-t2:t2
%             
%             if (i>=0 && j>=0)
%                 
%                 C_X(t1+i+1,t2+j+1)= 1/N*sum((X-mean(X)).*...
%                                     ([(X(i+1:end)-mean(X));zeros(i,1)]).*...
%                                     ([(X(j+1:end)-mean(X));zeros(j,1)]));
%             elseif (i<0 && j>=0)
%                 
%                 C_X(t1+i+1,t2+j+1)= 1/N*sum((X-mean(X)).*...
%                                     ([zeros(abs(i),1);(X(1:end+i)-mean(X))]).*...
%                                     ([(X(j+1:end)-mean(X));zeros(j,1)]));
%             elseif (i>=0 && j<0)
%                 
%                 C_X(t1+i+1,t2+j+1)= 1/N*sum((X-mean(X)).*...
%                                     ([(X(i+1:end)-mean(X));zeros(i,1)]).*...
%                                     ([zeros(abs(j),1);(X(1:end+j)-mean(X))]));
%             elseif (i<0 && j<0)
%                 
%                 C_X(t1+i+1,t2+j+1)= 1/N*sum((X-mean(X)).*...
%                                     ([zeros(abs(i),1);(X(1:end+i)-mean(X))]).*...
%                                     ([zeros(abs(j),1);(X(1:end+j)-mean(X))]));
%             end
%         end
%     end 
% else    
%     disp('Error! check Order or T');
%     disp('Order=2 -> T=t       &&      T<length(X)');
%     disp('Order=3 -> T=[t1,t2] && max(T)<length(X)');
% end
% 
% end





%%
% function C_X = cumulant(X,T)
% 
% % Cumulant(X,T) is a function to calculate 2nd and 3rd Cumulant which is
% % C(t) = E[ X(k)*X(k+t) ]  and C(t1,t2) = E[ X(k)*X(k+t1)*X(k+t2) ]
% % X is vector contain the series and T contains the lags values
% 
% N=length(X);
% if (length(T)==1)   
% %     C_X = 1/((N-T)*var(X))*sum((X-mean(X)).*([(X(T+1:end)-mean(X));zeros(T,1)]));            
%     C_X = 1/N*sum((X-mean(X)).*([(X(T+1:end)-mean(X));zeros(T,1)])); 
% elseif (length(T)==2)    
%     t1=T(1);
%     t2=T(2);
% %     C_X = 1/((N-max(T))*moment(X,3))*sum((X-mean(X)).*([(X(t1+1:end)-mean(X));zeros(t1,1)]).*([(X(t2+1:end)-mean(X));zeros(t2,1)]));
%     C_X = 1/N*sum((X-mean(X)).*([(X(t1+1:end)-mean(X));zeros(t1,1)]).*([(X(t2+1:end)-mean(X));zeros(t2,1)]));
% else    
%     disp('Error! check Order or T');
%     disp('Order=2 -> T=t       &&      T<length(X)');
%     disp('Order=3 -> T=[t1,t2] && max(T)<length(X)');
% end
% 
% end
% 
