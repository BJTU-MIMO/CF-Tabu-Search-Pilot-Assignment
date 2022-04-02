function [RR_cfta,Gammaa,Phii_cf]=functiontabu(K,M,N_iter,U,Phii_cf,BETAA,tau,Pp,Pu,LEN,Point,Monte,prelogFactor)
max_tabu=0;
    i_iter=1;
    TabuList=zeros(K,K+1);                      %(tabu list)
    TabuList(:,1)=inf;                      %(tabu list)
    Best=zeros(1,N_iter);
    while i_iter<N_iter
        %i_iter
        for i=1:K
            if TabuList(i,1)<=0
                TabuList(i,:)=[inf,zeros(1,K)];                      %(tabu list)
            end
        end
%         SE_ZF_0=functionZFcombiner(M,K,tau,Phii_cf,BETAA,Pp,Monte,prelogFactor,Pu);
        R_cfta=CalSINR(Phii_cf,M,K,BETAA,tau,Pp,Pu);
        [SR,which]=sort(R_cfta);
            %[M,which]=min(R);
        RRr=zeros(1,LEN);
        PPiont=zeros(LEN,K);
            for i=1:K
                mm(i)=sum(BETAA(i,:));
            end
            for ii=1:LEN
                [cfta,colo]=change(mm,which,ii,tau,Point,M,K,BETAA,Pp,Pu,U,Phii_cf,Monte,prelogFactor);
%                 [cfta,colo]=change_random(mm,which,ii,tau,Point,M,K,BETAA,Pp,Pu,U,Phii_cf);
                piont=colo;
                rr(ii,:)=cfta;
                RRr(ii)=sum(cfta);
                PPiont(ii,:)=piont;
            end
            [mr,index]=sort(RRr,'descend');
            %if和else之间表示特赦准则的情况
            %TabuList(i,:)=[which,chang_point];
              if mr(1)>max_tabu
                  max_tabu=mr(1);
                  Best(i_iter)=mr(1);
                  RR_cfta=rr(index(1),:);
                  Point=PPiont(index(1),:);
                  max_point=Point;
                  Phii_cf(:,which(index(1)))=U(:,PPiont(index(1),which(index(1))));
                  TabuList(:,1)=TabuList(:,1)-1;
%                   TabuList(max(exit),:)=[5,Point];
                  TS_list=TabuList(:,[2:K+1]);
                  exit=ismember(TS_list,Point,'rows');
                  if max(exit)==0
                     TabuList(min(find(TabuList(:,1)==inf)),:)=[5,Point];
                  else
                     TabuList(max(exit),:)=[5,Point];
                  end
%                   flag=999
                  tabu_index = i_iter;
              else
                 
                 for i=1:LEN
                    TS_list=TabuList(:,[2:K+1]);
                    eexit=ismember(TS_list,PPiont(index(i),:),'rows');
                    if max(eexit)==0
                        i
                        TabuList(:,1)=TabuList(:,1)-1;
                        TabuList(min(find(TabuList(:,1)==inf)),:)=[5,PPiont(index(i),:)];%加入这个分配
                        Best(i_iter)=mr(i);
                        RR_cfta=rr(index(i),:);
                        Point=PPiont(index(i),:);
                        Phii_cf(:,which(index(i)))=U(:,PPiont(index(i),which(index(i))));
                        break;
                    end
                end
%                     flag=2
              end
              i_iter=i_iter+1;
%               max_tabu=max(Best);
    end
for i=1:K
    Phii_cf(:,i) =U(:,max_point(i));
end
    Gammaa = zeros(M,K);
mau=zeros(M,K);
for m=1:M
    for k=1:K
        mau(m,k)=norm( (BETAA(m,:).^(1/2)).*(Phii_cf(:,k)'*Phii_cf))^2;
    end
end

for m=1:M
    for k=1:K
        Gammaa(m,k)=tau*Pp*BETAA(m,k)^2/(tau*Pp*mau(m,k) + 1);
    end
end
    
%     RC_CF_iter(n,:)=R_cf;