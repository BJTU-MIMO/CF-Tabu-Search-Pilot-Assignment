function [g,g_hat,C,y_p]=new_functionChannelEstimation(M,K,tau,Phii_cf,BETAA,Pp,Monte)

    h = sqrt(1/2)*(randn(Monte,M,K) +randn(Monte,M,K)*1i);
%     g = zeros(Monte,M,K);
    BETAA_mon = reshape(BETAA,[1,M,K]);
%     for monte=1:Monte
%         g(monte,:,:)=sqrt(BETAA_mon).*h(monte,:,:);  
%     end 
    g=sqrt(BETAA_mon).*h;  
    
    noise = sqrt(0.5*1)*(randn(Monte,M,K) + 1i*randn(Monte,M,K));
    pre_inter = zeros(Monte,K);
    inter = zeros(Monte,M,K);
    y_p = zeros(Monte,M,K);
    for m=1:M
        for k=1:K
            for kk=1:K
                pre_inter(:,kk)=g(:,m,kk)*(Phii_cf(:,k)'*Phii_cf(:,kk));
            end
            inter(:,m,k)=sum(pre_inter,2);
        end
    end
    
%     for k=1:K
%             for kk=1:K
%                 pre_inter(:,kk)=g(:,kk)*(Phii_cf(:,k)'*Phii_cf(:,kk));
%             end
%     end

%     inter=sum(pre_inter,2);
    for m=1:M            
        for k=1:K
%             inter(m,k)=sum(pre_inter);
            y_p(:,m,k)=sqrt(tau*Pp)*inter(:,m,k)+noise(:,m,k);
        end
    end
    
%     sqrt(tau*Pp).*g(:,:,kk)*(Phii_cf(:,k)'*Phii_cf(:,kk))+noise(:,:,)*Phii_cf(:,k);
    
    C = zeros(1,M,K);
    mau=zeros(M,K);
    for m=1:M
        for k=1:K
            mau(m,k)=norm((BETAA(m,:).^(1/2)).*(Phii_cf(:,k)'*Phii_cf))^2;
        end
    end
    for m=1:M
        for k=1:K
            C(1,m,k)=sqrt(tau*Pp)*BETAA(m,k)/(tau*Pp*mau(m,k) + 1);
        end
    end
    g_hat=C.*y_p;
%     g_hat=zeros(Monte,M,K);
%     for mon=1:Monte
%         g_hat(mon,:,:)=C.*y_p(mon,:,:);
%     end