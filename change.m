function  [R_cf,Point]=change(mm,which,ii,tau,Point,M,K,BETAA,Pp,Pu,U,Phii_cf,Monte,prelogFactor)
wu=zeros(tau,K);
for s=1:tau
    iii=1;
    for i=1:K
        if Point(i)==s
            wu(s,iii)=i;
            iii=iii+1;
        end
    end
end
none_zeros=(wu~=0);
for i=1:tau
    use_number(i)=sum(none_zeros(i,:));
end
[~,minindex]=min(use_number);
% mm(which(ii))=0;
% co=zeros(1,tau);
% [~,label]=size(wu);
% for i=1:tau
%     for j=1:label
%         if wu(i,j)~=0
%             co(i)=co(i)+mm(wu(i,j));
%         end
%     end
% end
%  co(Point(which(ii)))=1;
%  [col,colo]=min(co);
% % colo=randperm(S,1);
%  Point(which(ii))=colo;
%  Phii_cf(:,which(ii))=U(:,colo);
%  R_cf=CalSINR(Phii_cf,M,K,BETAA,tau,Pp,Pu);
 
 Point(which(ii))=minindex;
 Phii_cf(:,which(ii))=U(:,minindex);
 R_cf=CalSINR(Phii_cf,M,K,BETAA,tau,Pp,Pu);
% R_cf=functionZFcombiner(M,K,tau,Phii_cf,BETAA,Pp,Monte,prelogFactor,Pu);         
