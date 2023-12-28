function Cov=RC(cs_sig)

%Input data
%cs_sig: the formal error of any institution (SHCs*months).
%The order of cs_sig: S22,S21,C20,C21,C22,S33...

%Output data
%Reconstructed covariance matrices (SHCs*SHCs*months)


% load data
load Common_Eigenvector.mat; % Common eigenvectors: V_gfz06,V_tj2018.
load Fitted_Cofficients.mat; % Fitted coefficients from formal errors to eigenvalues: x_gfz,x_tj.
load Scale_Factor.mat; % Correction factors of fitted coffients.
load Fusion_Weight.mat; % Fusion weights for the reconstructed results from GFZ RL06 and Tongji-Grace2018 covariance

QQ=cs_sig.^2;

load time.mat;
t0=2002;
tt=(time-t0)';
p=size(cs_sig,1);
n=size(cs_sig,2);

for i=1:p
ss1(:,1)=log(QQ(i,:))';
A6(:,1)=ones(157,1);
A6(:,2)=ss1;
A6(:,3)=tt;
A6(:,4)=tt.*ss1;
A6(:,5)=tt.^2.*ss1;
A6(:,6)=tt.^2;

log_dd_tj2itsg(:,i)=A6*x_tj(:,i);
log_dd_gfz2itsg(:,i)=A6*x_gfz(:,i);
end

for i=1:2
Cov_tj2itsg=V_tj2018*diag(k_tj2itsg.*exp(log_dd_tj2itsg(i,:))')*V_tj2018';
Cov_gfz2itsg=V_gfz06*diag(k_gfz2itsg.*exp(log_dd_gfz2itsg(i,:))')*V_gfz06';
Cov(:,:,i)=P_tj2itsg(i)*Cov_tj2itsg+P_gfz2itsg(i)*Cov_gfz2itsg;
end
end