function [m] = action2m(r,m0,delta_T,delta_S)
%ACTION2M 此处显示有关此函数的摘要
%   此处显示详细说明
[T,N] = size(r);
N = N-1;
r = round(r*delta_T/delta_S);
m = zeros(T,N+1);
for i = 1:1:T
    if i==1
        for j=1:1:N+1
            new_j = r(i,j)+j;
            if new_j > N+1
                new_j = N+1;
            end
            m(i,new_j) = m(i,new_j)+m0(i,j);      
        end
%         m(i,:) = m(i,:)/sum(m(i,:));
    else
        for j=1:1:N+1
            new_j = r(i,j)+j;
            if new_j > N+1
                new_j = N+1;
            end
            m(i,new_j) = m(i,new_j)+m(i-1,j);           
        end
%         m(i,:) = m(i,:)/sum(m(i,:));
    end
end
end




% function [m] = action2m(r,m0,delta_T,delta_S,epsilon)
% %ACTION2M 此处显示有关此函数的摘要
% %   此处显示详细说明
% [T,N] = size(r);
% % r = round(r*delta_T/delta_S);
% 
% m = zeros(T,N);
% for i = 1:1:T
%     if i==1
%         for j=1:1:N
% %             new_j = r(i,j)+j;
% %             m(i,new_j) = m(i,new_j)+m0(i,j);
%             if j==1
%                 m(i,j) = m0(j) - delta_T/delta_S*(r(i,j+1)*m0(j+1)-r(i,j)*m0(j))+epsilon*(m0(j+1)-m0(j));
%             elseif j==N
%                 m(i,j) = m0(j) - delta_T/delta_S*(r(i,j)*m0(j)-r(i,j-1)*m0(j-1))+epsilon*(-m0(j)+m0(j-1));
%             else
%                 m(i,j) = m0(j) - delta_T/2/delta_S*(r(i,j+1)*m0(j+1)-r(i,j-1)*m0(j-1))+epsilon*(m0(j+1)-2*m0(j)+m0(j-1));
%             end           
%         end
%         m(i,:) = m(i,:)/sum(m(i,:));
%     else
%         for j=1:1:N
% %             new_j = r(i,j)+j;
% %             m(i,new_j) = m(i,new_j)+m(i-1,j);           
%             if j==1
%                 m(i,j) = m(i-1,j) - delta_T/delta_S*(r(i,j+1)*m(i-1,j+1)-r(i,j)*m(i-1,j))+epsilon*(m(i-1,j+1)-m(i-1,j));
%             elseif j==N
%                 m(i,j) = m(i-1,j) - delta_T/delta_S*(r(i,j)*m(i-1,j)-r(i,j-1)*m(i-1,j-1))+epsilon*(-m(i-1,j)+m(i-1,j-1));
%             else
%                 m(i,j) = m(i-1,j) - delta_T/2/delta_S*(r(i,j+1)*m(i-1,j+1)-r(i,j-1)*m(i-1,j-1))+epsilon*(m(i-1,j+1)-2*m(i-1,j)+m(i-1,j-1));
%             end
%         end
%         m(i,:) = m(i,:)/sum(m(i,:));
%     end
% end
% % m = m.*(m>0)-m.*(m<=0);
% % for i = 1:1:24
% %     m(i,:) = m(i,:)/sum(m(i,:));
% % end
% end
