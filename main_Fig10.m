clear; clc;
rng(3);
tic;

% Parameters
NR = 4;
NT = 4;
Pt = 1;
SNR_dBs = -10:5:20;
s2 = Pt ./ db2pow(SNR_dBs);
nMonte = 1e4;
Y0 = 1/50;

Bs = [Inf,8,4,2];

SR_RZFBF = zeros(nMonte,length(Bs),length(SNR_dBs),2);



for iMonte = 1:nMonte
    H = sqrt(1/2) * (randn(NR,NT) + 1i * randn(NR,NT));

    for iB = 1:length(Bs)
        Hq_real = reshape(func_quant_Gauss(real(H(:)),Bs(iB)/2,sqrt(1/2)), NR, NT); % B/2 bits allocated to real part
        Hq_imag = reshape(func_quant_Gauss(imag(H(:)),Bs(iB)/2,sqrt(1/2)), NR, NT); % B/2 bits allocated to imag part
        Hq = Hq_real + 1i * Hq_imag;

        for iSNR = 1:length(SNR_dBs)
    
            % R-ZFBF: Digital
            F_RZFBF = Hq'/(Hq*Hq'+NR/db2pow(SNR_dBs(iSNR))*eye(NR));
            W_RZFBF = sqrt(Pt/NR) * F_RZFBF ./ vecnorm(F_RZFBF);
            %W_RZFBF = sqrt(Es) * F_RZFBF / norm(F_RZFBF,'fro');
            
            Pq = [[-NR/db2pow(SNR_dBs(iSNR))*eye(NR), Hq];[Hq', eye(NT)]];
            Yik_offdiag = - Y0 * Pq;
            Yik_diag = Y0 * sum(Pq) - Y0;
            Yik = Yik_offdiag - diag(diag(Yik_offdiag)) + diag(Yik_diag);
 
            Y_offdiag = - Yik;
            Y_diag = sum(Yik);
            Y = Y_offdiag - diag(diag(Y_offdiag)) + diag(Y_diag);

            % R-ZFBF: MiLAC (LMMSE)
            %Q_RZFBF = inv([[-NR/db2pow(SNR_dBs(iSNR))*eye(NR),Hq];[Hq',eye(NT)]]);
            Q_RZFBF = inv(Y/Y0 + eye(NT+NR));
            Q21_RZFBF = sqrt(Pt) * Q_RZFBF(NR+1:end,1:NR) / norm(Q_RZFBF(NR+1:end,1:NR),'fro'); % norm(Q_RZFBF(K+1:end,1:K),'fro') should be computed!
            
            % Sum Rate
            SR_RZFBF(iMonte,iB,iSNR,1) = func_sum_rate(W_RZFBF,H,s2(iSNR));
            SR_RZFBF(iMonte,iB,iSNR,2) = func_sum_rate(Q21_RZFBF,H,s2(iSNR));
        end
    end
end

toc;

SR_RZFBF_mean = squeeze(mean(SR_RZFBF));

%% Plot
figure('DefaultAxesFontSize',12);
LineW = 1.8; MarkS = 8;
plot(SNR_dBs,SR_RZFBF_mean(1,:,1),'-h','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Digital');
set(gca,'ColorOrderIndex',5)
hold on;
plot(SNR_dBs,SR_RZFBF_mean(1,:,2),'--^','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC, Continuous');
plot(SNR_dBs,SR_RZFBF_mean(2,:,2),'-->','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC,{\it B = 8}');
plot(SNR_dBs,SR_RZFBF_mean(3,:,2),'--<','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC,{\it B = 4}');
plot(SNR_dBs,SR_RZFBF_mean(4,:,2),'--v','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC,{\it B = 2}');

grid on;
set(gca,'GridLineStyle',':','GridAlpha',0.5,'LineWidth',1.2);
xlabel('SNR [dB]');
ylabel('Sum rate [bps/Hz]');
legend('Location','northwest','NumColumns',1,'FontSize',12);

set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1.5);