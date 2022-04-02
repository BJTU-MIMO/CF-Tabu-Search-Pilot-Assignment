clc
clear all

%Uplink number of terminals
M=100; %number of access points
B=20; %Mhz
D=1; %length of square area
%Consider a rectangular area with DxD km^2
%M distributed APs serves K terminals, they all randomly located in the area
K=40; % number of user equipments
tau_p=30;% number of pilots
tau_c=200;
% prelogFactor =B*(1-tau_p/tau_c)*(1/2);
prelogFactor =1;

[U,S,V]=svd(randn(tau_p,tau_p)+randn(tau_p,tau_p)*1i);%U includes tau_p orthogonal sequences 
[UU,SS,VV]=svd(randn(K,K)+randn(K,K)*1i);%U includes tau_p orthogonal sequences 

power_f=0.1; %downlink power: 200 mW
noise_p = 10^((-203.975+10*log10(20*10^6)+9)/10); %noise power
Pu = power_f/noise_p;%nomalized receive SNR
Pp=Pu;%pilot power
Pd=2*Pu;

N=200;
Monte=1000;  %number of monte carlo simulation
minn=zeros(1,N);
nn=zeros(1,N);
DS = zeros(N,K); %desired signal
BU = zeros(N,K); %beamforming uncerntenty
UI = zeros(N,K); %user interference
Extra_noise = zeros(N,K); %noise

R_MR_gra = zeros(N,K); %user interference
R_MR_ran = zeros(N,K); %user interference
R_MR_ortho = zeros(N,K); %user interference
R_MR_ex = zeros(N,K); %user interference
R_MR_greedy = zeros(N,K); %user interference

DS_mon = zeros(Monte,K); %desired signal
first_mon = zeros(Monte,K); %beamforming uncerntenty
inter_sum=zeros(Monte,K);
UI_mon = zeros(Monte,K); %user interference
Extra_noise_mon = zeros(Monte,K); %noise
SE_ZF_gra = zeros(N,K); %
SE_ZF_ran = zeros(N,K); %
SE_ZF_ortho = zeros(N,K); %
SE_ZF_ex = zeros(N,K); %
SE_ZF_greedy = zeros(N,K); %
SE_ZF_tabu = zeros(N,K); %
SE_MR_tabu = zeros(N,K); %
n=1;
Point=zeros(1,K);
piont=zeros(1,K);
N_iter=20;
LEN=K;
max_tabu=0;
while n < N+1
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(N)]);
    
    %Generate channel matrix
    [BETAA]=functionChannelScene(D,M,K,Monte);
    %Random pilot assignment
%     prelogFactor =B*(1-tau_p/tau_c)*(1/2);
%     Phii=zeros(tau_p,K);
%     for k=1:K
%         Point=randi([1,tau_p]);
%         Phii(:,k)=U(:,Point);
%     end
%     Phii_cf = Phii;
    
   
    Phii=zeros(tau_p,K);
    for k=1:K
        Point(k)=randi([1,tau_p]);
        Phii(:,k)=U(:,Point(k));
    end

    Phii_cf = Phii;
    
    
    [R_cf_ran_2,Gammaa]=CalSINR(Phii_cf,M,K,BETAA,tau_p,Pp,Pu);
%     R_MR_ran(n,:)=prelogFactor*R_cf_ran_2;
    R_cf_min_n=min(R_cf_ran_2);
    R_MR_ran(n,:)=cvx_opti(R_cf_min_n,K,BETAA,Pu,Gammaa,Phii_cf);
    R_MR_ran(n,:)=prelogFactor*R_MR_ran(n,:);
    
    rrr=zeros(1,K);
%% greedy pilot assignment
stepp=100;
Ratestep=zeros(stepp,K);
Ratestep(1,:)=R_cf_ran_2;
%%% Find the pilot sequences using greedy method
for st=2:stepp
    [minvalue minindex]=min(Ratestep(st-1,:));
    if (minn(n)==minindex)
        break;
    end
    minn(n)=minindex;
    nn(n)=nn(n)+1;
    %[maxvalue maxindex]=max(Ratestep(st-1,:));
    
        Mat=zeros(tau_p,tau_p)-Pu*sum(BETAA(:,minindex))*Phii_cf(:,minindex)*Phii_cf(:,minindex)';
        %Mat2=zeros(tau_p,tau_p)-Pu*sum(BETAA(:,maxindex))*Phii(:,maxindex)*Phii(:,maxindex)';
        for kk=1:K
            Mat = Mat + Pu*sum(BETAA(:,kk))*Phii_cf(:,kk)*Phii_cf(:,kk)';
            %Mat2 = Mat2 + Pu*sum(BETAA(:,kk))*Phii(:,kk)*Phii(:,kk)';
        end
        [U1,S1,V1] = svd(Mat);
        %[U2,S2,V2] = svd(Mat2);
        Phii_cf(:,minindex) = U1(:,tau_p);
        %Phii(:,maxindex) = U2(:,1);
    
 
  %% Create Gamma matrix
Gammaa = zeros(M,K);
mau=zeros(M,K);
for m=1:M
    for k=1:K
        mau(m,k)=norm( (BETAA(m,:).^(1/2)).*(Phii_cf(:,k)'*Phii_cf))^2;
    end
end

for m=1:M
    for k=1:K
        Gammaa(m,k)=tau_p*Pp*BETAA(m,k)^2/(tau_p*Pp*mau(m,k) + 1);
    end
end

%% Compute etaa(m): (each AP transmits equal power to K terminals)
etaa=zeros(M,1);
for m=1:M
    etaa(m)=1/(sum(Gammaa(m,:)));
end
%% Compute downlink Rate
SINR=zeros(1,K);

%Pilot contamination
PC = zeros(K,K);
for ii=1:K
    for k=1:K
%         for m=1:M
%         PC(ii,k) = PC(ii,k) + (etaa(m)^(1/2))*Gammaa(m,ii)*(BETAA(m,k)/BETAA(m,ii))*Phii(:,ii)'*Phii(:,k);
%         end
        PC(ii,k)=sum((etaa.^(1/2)).*((Gammaa(:,ii)./BETAA(:,ii)).*BETAA(:,k)))*Phii_cf(:,ii)'*Phii_cf(:,k);
    end
end
PC1=(abs(PC)).^2;

for k=1:K
%     num=0;
%     for m=1:M
%         num=num + (etaa(m)^(1/2))*Gammaa(m,k);
%     end
    num = (etaa.^(1/2))'*Gammaa(:,k);
    SINR(k) = Pd*num^2/(1 + Pd*sum(BETAA(:,k)) + Pd*sum(PC1(:,k)) - Pd*PC1(k,k));
    %Rate:
    Ratestep(st,k) = log2(1+ SINR(k));
end

       
end
    [R_cf_greedy,Gammaa]=CalSINR(Phii_cf,M,K,BETAA,tau_p,Pp,Pu);
%     R_MR_greedy(n,:)=prelogFactor*R_cf_greedy;
    R_cf_min_n=min(R_cf_greedy);
    R_MR_greedy(n,:)=cvx_opti(R_cf_min_n,K,BETAA,Pu,Gammaa,Phii_cf);
    R_MR_greedy(n,:)=prelogFactor*R_MR_greedy(n,:);
%         rrr=zeros(1,K);

    %% orthogonal
    Phii_ortho=zeros(K,K);
    for k=1:K
        Phii_ortho(:,k)=UU(:,k);
    end
%     prelogFactor=B*(1-K/tau_c)*(1/2);
%     SE_ZF_4=functionZFcombiner(M,K,tau_p,Phii_ortho,BETAA,Pp,Monte,prelogFactor,Pu);
%     SE_ZF_ortho(n,:)=SE_ZF_4;
    
    [R_cf_ran_4,Gammaa]=CalSINR(Phii_ortho,M,K,BETAA,K,Pp,Pu);
        R_cf_min_n=min(R_cf_ran_4);
    R_MR_ortho(n,:)=cvx_opti(R_cf_min_n,K,BETAA,Pu,Gammaa,Phii_ortho);
    R_MR_ortho(n,:)=prelogFactor*R_MR_ortho(n,:);
% %     
% %% exhausted
%     max_ex=0;
%     max_A=zeros(1,K);
%     maxR_ex=zeros(1,K);
%     maxPhii_ex=zeros(tau_p,K);
%     
%     max_ex_ZF=0;
%     max_A_ZF=zeros(1,K);
%     maxR_ex_ZF=zeros(1,K);
%     maxPhii_ex_ZF=zeros(tau_p,K);
%     for num_ex=0:(tau_p^K-1)
%         A=dec2base(num_ex,tau_p);
%         AAA=zeros(1,K);
%         numberr=0;
%         for ii=K:-1:2
%             if length(A)>numberr
%                 AAA(ii)=str2num(A(end-(K-ii)));
%                 numberr=numberr+1;
%             end
%         end
%         AAA=AAA+1;
%         Phii_ex=zeros(tau_p,K);
%         for k=1:K
%             Phii_ex(:,k)=U(:,AAA(k));
%         end
%         [R_cf_5,Gammaa]=CalSINR(Phii_ex,M,K,BETAA,tau_p,Pp,Pu);
%         if sum(R_cf_5)>max_ex
%             max_ex=sum(R_cf_5);
%             max_A=AAA;
%             maxR_ex=R_cf_5;
%             maxPhii_ex = Phii_ex;
%             maxGammaa = Gammaa;
%         end
%         
% %         [R_cf_ZF,Gammaa]=CalSINR(Phii_ex,M,K,BETAA,tau_p,Pp,Pu);
% %         if sum(R_cf_ZF)>max_ex_ZF
% %             max_ex_ZF=sum(R_cf_ZF);
% %             max_A_ZF=AAA;
% %             maxR_ex_ZF=R_cf_ZF;
% %             maxPhii_ex_ZF = Phii_ex;
% %             maxGammaa_ZF = Gammaa;
% %         end
%     end
% %     SE_ZF_ex(n,:)=maxR_ex_ZF;
% %     R_MR_ex(n,:)=maxR_ex;
%     R_cf_min_n=min(maxR_ex);
%     R_MR_ex(n,:)=cvx_opti(R_cf_min_n,K,BETAA,Pu,maxGammaa,maxPhii_ex);
%     R_MR_ex(n,:)=prelogFactor*R_MR_ex(n,:);


    n=n+1;
end
prelogFactor =B*(1-tau_p/tau_c)*(1/2);
prelogFactor1 =B*(1-K/tau_c)*(1/2);
% prelogFactor=1;
% prelogFactor1 =1;
Y1=linspace(0,1,N*K);
SE_CF_ZF_tabu=reshape(SE_ZF_tabu,N*K,1);
SE_CF_ZF_gra=reshape(SE_ZF_gra,N*K,1);
SE_CF_ZF_greedy=reshape(SE_ZF_greedy,N*K,1);
SE_CF_ZF_ran = reshape(SE_ZF_ran,N*K,1);
SE_CF_ZF_ortho =reshape(SE_ZF_ortho,N*K,1);
SE_CF_ZF_ex=reshape(SE_ZF_ex,N*K,1);
SE_CF_MR_ran=reshape(R_MR_ran,N*K,1);
SE_CF_MR_greedy =reshape(R_MR_greedy,N*K,1);
SE_CF_MR_graph =reshape(R_MR_gra,N*K,1);
SE_CF_MR_ex=reshape(R_MR_ex,N*K,1);
SE_CF_MR_ortho =reshape(R_MR_ortho,N*K,1);
SE_CF_MR_tabu =reshape(SE_MR_tabu,N*K,1);

ZF_tabu = sort(SE_CF_ZF_tabu);
ZF_gra = sort(SE_CF_ZF_gra);
ZF_greedy = sort(SE_CF_ZF_greedy);
ZF_ran = sort(SE_CF_ZF_ran);
ZF_ortho = sort(SE_CF_ZF_ortho);
ZF_ex = sort(SE_CF_ZF_ex);

MR_ran = sort(SE_CF_MR_ran);
MR_greedy= sort(SE_CF_MR_greedy);
MR_graph= sort(SE_CF_MR_graph);
MR_ex= sort(SE_CF_MR_ex);
MR_ortho= sort(SE_CF_MR_ortho);
MR_tabu= sort(SE_CF_MR_tabu);

figure

line2=plot(prelogFactor*MR_ran,Y1(:),'k--','LineWidth',1.5);
hold on

line4=plot(prelogFactor*MR_greedy,Y1(:),'b-.','LineWidth',1.5);
hold on

line6=plot(prelogFactor*MR_tabu,Y1(:),'r--','LineWidth',1.5);
hold on
line10=plot(prelogFactor*MR_ex,Y1(:),'k-','LineWidth',1.5);
hold on

line8=plot(prelogFactor*MR_graph,Y1(:),'g--','LineWidth',1.5);
hold on
line9=plot(prelogFactor1*MR_ortho,Y1(:),'r:','LineWidth',1.5);
line=[line2,line4,line6,line8,line9,line10];
xlabel(['Per-User Uplink SE (bits/s/Hz), number of tau_p£∫', num2str(tau_p), '  , M/K:',num2str(M/K), '  , K/tau_p:',num2str(K/tau_p)])
ylabel('Cumulative Distribution')
legend(line,{'Random [1]','Greedy [1]','TS','Graph','Orthogonal','Exhaustive'},'location','northwest');%ÃÌº”±Í ∂
grid on
