function [R_cf,Gammaa]=CalSINR(Phii_cf,M,K,BETAA,tau,Pp,Pu)
%Phii_cf = Phii; % pilot set of cell-free systems
% %% Create Gamma matrix
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


%%% 1) Each Terminal transmits full eta_k=1

%%% Compute Rate
SINR=zeros(1,K);
R_cf=zeros(1,K);

%Pilot contamination
PC = zeros(K,K);
for ii=1:K
    for k=1:K
        PC(ii,k) = sum( (Gammaa(:,k)./BETAA(:,k)).*BETAA(:,ii)  )*Phii_cf(:,k)'*Phii_cf(:,ii);
    end
end
PC1=(abs(PC)).^2;

for k=1:K
    deno1=0;
    for m=1:M
        deno1=deno1 + Gammaa(m,k)*sum(BETAA(m,:));
    end
    SINR(k) = Pu*(sum(Gammaa(:,k)))^2/(sum(Gammaa(:,k)) + Pu*deno1 + Pu*sum(PC1(:,k)) - Pu*PC1(k,k));
    %Rate:
    R_cf(k) = log2(1+ SINR(k));
end