function R_cf_opt_user=cvx_opti(R_cf_min,K,BETAA,Pu,Gammaa,Phii_cf)
R_cf_opt_user=zeros(1,K);
%% 2) Max-Min power allocations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmin=2^R_cf_min-1;
% if R_cf_min(n)>0.5
%     tmax=2^(4*R_cf_min(n))-1;
% elseif ((R_cf_min(n)<=0.5) && (R_cf_min(n)>=0.3))
%     tmax=2^(R_cf_min(n)+1.2)-1;
% else
%     %tmax=tmin+0.6;
%     tmax=2^(R_cf_min(n)+0.8)-1;
% end
%epsi=max(0.1,(tmax-tmin)/8);

tmax=2^(2*R_cf_min+1.2)-1;
epsi=max(tmin/5,0.01);

BETAAn=BETAA*Pu;
Gammaan=Gammaa*Pu;
PhiPhi=zeros(K,K);
Te1 =zeros(K,K);
Te2 =zeros(K,K);

x_cf_opt=ones(K,1);

for ii=1:K
    for k=1:K
        PhiPhi(ii,k)=norm(Phii_cf(:,ii)'*Phii_cf(:,k));
    end
end

for ii=1:K
    for k=1:K
        Te1(ii,k)=sum(BETAAn(:,ii).*Gammaan(:,k));
        Te2(ii,k)=(sum((Gammaan(:,k)./BETAA(:,k)).*BETAA(:,ii)) )^2*PhiPhi(k,ii)^2;   
    end
end

%cvx_solver sedumi
cvx_quiet true
            while( tmax - tmin > epsi)

            tnext = (tmax+tmin)/2; 
           cvx_begin %sdp
              variables x(K,1) 
              minimize(0)
              subject to
                for k=1:K
                    Te1(:,k)'*x + [Te2(1:(k-1),k); Te2((k+1):K,k) ]'*[x(1:(k-1)); x((k+1):K)] + sum(Gammaan(:,k)) <= (1/tnext)*(sum(Gammaan(:,k)))^2*x(k) ;
                end              
                for k=1:K
                    x(k)>=0;
                    x(k)<=1;
                end
                            
            cvx_end


            % bisection
            if strfind(cvx_status,'Solved') % feasible
            fprintf(1,'Problem is feasible ',tnext);
            tmin = tnext;
            x_cf_opt=x;
            else % not feasible
            fprintf(1,'Problem not feasible ',tnext);
            tmax = tnext;   
            end

            end

            for k=1:K
                mu1k=(sum(Gammaan(:,k)))^2*x_cf_opt(k)/( Te1(:,k)'*x_cf_opt + [Te2(1:(k-1),k); Te2((k+1):K,k) ]'*[x_cf_opt(1:(k-1)); x_cf_opt((k+1):K)] + sum(Gammaan(:,k)) );
                R_cf_opt_user(1,k) = log2(1+mu1k);
            
            end              
            
%R_cf_opt_min(n) = log2(1+tmin);
