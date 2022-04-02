function [g_hat,C,y_p]=functionChannelEstimation(M,K,tau,g,Phii_cf,BETAA,Pp,Monte)
% noise = zeros(M,tau);
%     for m=1:M
%         for i=1:tau
%             noise(m,i) =sqrt(1/2)*(randn(1,1) +randn(1,1)*1i);
%         end
%     end

    
    noise = sqrt(0.5)*(randn(tau,M) + 1i*randn(tau,M));
    pre_inter = zeros(1,K);
    inter = zeros(M,K);
    y_p = zeros(M,K);
    for m=1:M
        for k=1:K
            for kk=1:K
                pre_inter(1,kk)=g(m,kk)*(Phii_cf(:,k)'*Phii_cf(:,kk));
            end
            inter(m,k)=sum(pre_inter);
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
            y_p(m,k)=sqrt(tau*Pp)*inter(m,k)+Phii_cf(:,k)'*noise(:,m);
        end
    end
    
%     sqrt(tau*Pp).*g(:,:,kk)*(Phii_cf(:,k)'*Phii_cf(:,kk))+noise(:,:,)*Phii_cf(:,k);
    
    C = zeros(M,K);
    mau=zeros(M,K);
    for m=1:M
        for k=1:K
            mau(m,k)=norm((BETAA(m,:).^(1/2)).*(Phii_cf(:,k)'*Phii_cf))^2;
        end
    end
    for m=1:M
        for k=1:K
            C(m,k)=sqrt(tau*Pp)*BETAA(m,k)/(tau*Pp*mau(m,k) + 1);
        end
    end
    g_hat=C.*y_p;