% Optimize the SNR-constrained joint beamforming and reflection design with hybrid beamforming.
% This is used in the paper: R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, “SNR/CRB-constrained joint beamforming and reflection designs for RIS-ISAC systems,”IEEE Trans. Wireless Commun., to appear.
% Inputs: Prms: the structure of system parameters (requires Nrf and hybrid_maxiter);
%         Channel: the structure of the channels;
%         phi: the initial phi; F_RF/F_BB: the initial hybrid beamformers
% Outputs: W: hybrid beamforming; F_RF: analog beamformer; F_BB: digital beamformer
%          phi: RIS reflection coefficients; Vsr: the achieved sum-rate; gammat_r: the achieved radar SNR
function [W,F_RF,F_BB,phi,Vsr,gammat_r] = get_W_phi_SNR_hybrid(Prms,Channel,phi,F_RF,F_BB)

M = Prms.M; N = Prms.N; K = Prms.K; sigmar2 = Prms.sigmar2; sigmak2 = Prms.sigmak2;
sigmat2 = Prms.sigmat2; Nmax = Prms.Nmax; res_th = Prms.res_th; nL = Prms.L;
gammat = Prms.gammat; P = Prms.P; hdt = Channel.hdt;hrt = Channel.hrt; G = Channel.G; Hu = Channel.Hu; Hru = Channel.Hru;
Nrf = Prms.Nrf; hybrid_maxiter = Prms.hybrid_maxiter;

W = F_RF*F_BB;
if norm(W,'fro')^2 > P
    scale = sqrt(P)/norm(W,'fro');
    F_BB = scale*F_BB;
    W = F_RF*F_BB;
end

%%%% variable follow phi
Hk = Hu + Hru*diag(phi)*G;
Ht = (hdt + G.'*diag(phi)*hrt)*(hdt.' + hrt.'*diag(phi)*G);
r = zeros(K,1);
c = zeros(K,1);
for k = 1:1:K
    r(k) = abs(Hk(k,:)*W(:,k))^2/(norm(Hk(k,:)*W,2)^2-abs(Hk(k,:)*W(:,k))^2+sigmak2);
    c(k) = sqrt(1+r(k))*Hk(k,:)*W(:,k)/(norm(Hk(k,:)*W,2)^2+sigmak2);
end
u = kron(eye(K+M),Ht)*vec(W)/((vec(W))'*(kron(eye(K+M),Ht'*Ht))*vec(W));
e3 = sqrt(real(gammat*sigmar2*u'*u)/sigmat2/nL);
D = zeros(N,N);
g = zeros(N,1);
for k = 1:1:K
    g = g + 2*sqrt(1+r(k))*c(k)*diag(W(:,k)'*G')*Hru(k,:)';
    for j = 1:1:K+M
        temp = diag(W(:,j)'*G')*Hru(k,:)';
        D = D + abs(c(k))^2*(temp*temp');
        g = g - 2*abs(c(k))^2*diag(W(:,j)'*G')*Hru(k,:)'*Hu(k,:)*W(:,j);
    end
end

%%%% start iteration
iter = 1;
res = 1;
Vsr = zeros(1,Nmax);
Vres = zeros(1,Nmax);
Vgammat = zeros(1,Nmax);
mu = zeros(N,1);
varphi = phi;
rho = N/abs(phi'*D*phi-real(g'*phi));
while iter <= Nmax && res >= res_th
    %%%% update reflecting coefficients phi
    D = zeros(N,N);
    g = zeros(N,1);
    for k = 1:1:K
        g = g + 2*sqrt(1+r(k))*c(k)*diag(W(:,k)'*G')*Hru(k,:)';
        for j = 1:1:K+M
            temp = diag(W(:,j)'*G')*Hru(k,:)';
            D = D + abs(c(k))^2*(temp*temp');
            g = g - 2*abs(c(k))^2*diag(W(:,j)'*G')*Hru(k,:)'*Hu(k,:)*W(:,j);
        end
    end
    F = kron(W.'*hdt,G.'*diag(hrt)) + kron(W.'*G.'*diag(hrt),hdt);
    L = kron(W.'*G.'*diag(hrt),G.'*diag(hrt));
    Ltilde = reshape(L.'*conj(u),N,N);
    Lbar = [-real(Ltilde) imag(Ltilde);imag(Ltilde) real(Ltilde)];
    [~,bb] = eigsort(Lbar+Lbar.');
    lambda = bb(1,1);
    phibar = [real(phi);imag(phi)];
    uth = -u'*F + phibar.'*(Lbar+Lbar.'-lambda*eye(2*N))*[eye(N);-1i*eye(N)];
    e4 = -e3 + phibar.'*Lbar.'*phibar + real(u'*kron(eye(K+M),hdt*hdt.')*vec(W))-lambda*N;
    cvx_begin quiet
    cvx_solver SeDuMi
    variable phi(N,1) complex
    minimize real(phi'*D*phi)-real(g'*phi) + 0.5/rho*square_pos(norm(phi-varphi+mu*rho,2))
    subject to
    real(uth*phi) <= e4;
    abs(phi) <= 1;
    cvx_end

    %%%% update varphi
    varphi = exp(1i*angle(phi+mu*rho));

    %%% update mu
    mu = mu + (phi-varphi)/rho;

    %%%% variable follow phi
    Hk = Hu + Hru*diag(phi)*G;
    Ht = (hdt + G.'*diag(phi)*hrt)*(hdt.' + hrt.'*diag(phi)*G);

    %%%%% update transmit beamformer w (full-digital surrogate)
    a = zeros(M*(K+M),1);
    B = zeros(K*(K+M),M*(K+M));
    for k = 1:1:K
        a((k-1)*M+1:k*M) = 2*sqrt(1+r(k))*c(k)*Hk(k,:)';
        for j = 1:1:K+M
            Tj = zeros(M,M*(K+M));
            Tj(:,(j-1)*M+1:j*M) = eye(M);
            B((k-1)*(K+M)+j,:) = abs(c(k))*Tj.'*Hk(k,:).';
        end
    end
    kHt = kron(eye(K+M),Ht);
    cvx_begin quiet
    cvx_solver SeDuMi
    variable w(M*(K+M),1) complex
    minimize square_pos(norm(B*w,2))-real(a'*w)
    subject to
    real(u'*kHt*w) >= e3;
    norm(w,2) <= sqrt(P);
    cvx_end
    Wd = reshape(w,M,K+M);

    %%%%% hybrid decomposition
    [F_RF,F_BB,W] = hybrid_decompose(Wd,Nrf,hybrid_maxiter,F_RF);
    if norm(W,'fro')^2 > P
        scale = sqrt(P)/norm(W,'fro');
        F_BB = scale*F_BB;
        W = F_RF*F_BB;
    end

    %%%% update receive filter u
    u = kron(eye(K+M),Ht)*vec(W)/((vec(W))'*(kron(eye(K+M),Ht'*Ht))*vec(W));
    e3 = sqrt(real(gammat*sigmar2*u'*u)/sigmat2/nL);

    %%%% update auxiliary variables
    r = zeros(K,1);
    c = zeros(K,1);
    for k = 1:1:K
        r(k) = abs(Hk(k,:)*W(:,k))^2/(norm(Hk(k,:)*W,2)^2-abs(Hk(k,:)*W(:,k))^2+sigmak2);
        c(k) = sqrt(1+r(k))*Hk(k,:)*W(:,k)/(norm(Hk(k,:)*W,2)^2+sigmak2);
    end

    Vsr(iter) = sum(log2(1+r));
    Vres(iter) = norm(phi-varphi,2)^2;
    Vgammat(iter) = real(nL*sigmat2*abs(u'*kron(eye(K+M),Ht)*vec(W))^2/(sigmar2*u'*u));
    if iter > 10
        res = abs(Vsr(iter)-Vsr(iter-1))/Vsr(iter-1);
    end
    iter = iter + 1;
    rho = rho*0.8;
end
Vsr(iter:end) = [];
Vres(iter:end) = [];
Vgammat(iter:end) = [];
gammat_r = real(nL*sigmat2*abs(u'*kron(eye(K+M),Ht)*vec(W))^2/(sigmar2*u'*u));

end
