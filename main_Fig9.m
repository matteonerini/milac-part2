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

SNRH_dBs = [Inf,20,10,5];

SR_RZFBF = zeros(nMonte,length(SNRH_dBs),length(SNR_dBs),2);



for iMonte = 1:nMonte
    H = sqrt(1/2) * (randn(NR,NT) + 1i * randn(NR,NT));

    for iSNRH = 1:length(SNRH_dBs)
        Hn = H + sqrt(1/(2*db2pow(SNRH_dBs(iSNRH)))) * (randn(NR,NT) + 1i * randn(NR,NT));
    
        for iSNR = 1:length(SNR_dBs)
    
            % R-ZFBF: Digital
            F_RZFBF = Hn'/(Hn*Hn'+NR/db2pow(SNR_dBs(iSNR))*eye(NR));
            W_RZFBF = sqrt(Pt/NR) * F_RZFBF ./ vecnorm(F_RZFBF);
            %W_RZFBF = sqrt(Es) * F_RZFBF / norm(F_RZFBF,'fro');
    
            % R-ZFBF: MiLAC (LMMSE)
            Q_RZFBF = inv([[-NR/db2pow(SNR_dBs(iSNR))*eye(NR),Hn];[Hn',eye(NT)]]);
            Q21_RZFBF = sqrt(Pt) * Q_RZFBF(NR+1:end,1:NR) / norm(Q_RZFBF(NR+1:end,1:NR),'fro');
            
            % Sum Rate
            SR_RZFBF(iMonte,iSNRH,iSNR,1) = func_sum_rate(W_RZFBF,H,s2(iSNR));
            SR_RZFBF(iMonte,iSNRH,iSNR,2) = func_sum_rate(Q21_RZFBF,H,s2(iSNR));
        end
    end
end

toc;

SR_RZFBF_mean = squeeze(mean(SR_RZFBF));

%% Plot
figure('DefaultAxesFontSize',12);
LineW = 1.8; MarkS = 8;
plot(SNR_dBs,SR_RZFBF_mean(1,:,1),'-h','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Digital');
hold on;
plot(SNR_dBs,SR_RZFBF_mean(2,:,1),'-d','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Digital');
plot(SNR_dBs,SR_RZFBF_mean(3,:,1),'-s','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Digital');
plot(SNR_dBs,SR_RZFBF_mean(4,:,1),'-o','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Digital');

plot(SNR_dBs,SR_RZFBF_mean(1,:,2),'--^','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC, Perfect CSI');
plot(SNR_dBs,SR_RZFBF_mean(2,:,2),'-->','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC,{\it \rho = 20} dB');
plot(SNR_dBs,SR_RZFBF_mean(3,:,2),'--<','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC,{\it \rho = 10} dB');
plot(SNR_dBs,SR_RZFBF_mean(4,:,2),'--v','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC,{\it \rho = 5} dB');
grid on;
set(gca,'GridLineStyle',':','GridAlpha',0.5,'LineWidth',1.2);
xlabel('SNR [dB]');
ylabel('Sum rate [bps/Hz]');
legend('Location','northwest','NumColumns',2,'FontSize',12);

set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1.5);