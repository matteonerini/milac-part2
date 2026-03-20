clear; clc;
rng(3);
tic;

% Parameters
NR = 4;
NT = 4;
Pt = 1;
SNR_dBs = -10:5:20;
s2 = Pt ./ db2pow(SNR_dBs);
nMonte = 1e5;

SR_MBF = zeros(nMonte,length(SNR_dBs),3);
SR_ZFBF = zeros(nMonte,length(SNR_dBs),3);
SR_RZFBF = zeros(nMonte,length(SNR_dBs),3);

SR_ZFBF_teor = zeros(nMonte,length(SNR_dBs));

for iMonte = 1:nMonte
    H = sqrt(1/2) * (randn(NR,NT) + 1i * randn(NR,NT));

    % MBF: Digital
    F_MBF = H';
    W_MBF = sqrt(Pt/NR) * F_MBF ./ vecnorm(F_MBF);
    %W_MBF = sqrt(Pt) * F_MBF / norm(F_MBF,'fro');

    % MBF: MiLAC (arbit.)
    Q_MBF_arbit = inv([[eye(NR),zeros(NR,NT)];[-W_MBF,eye(NT)]]);
    Q21_MBF_arbit = Q_MBF_arbit(NR+1:end,1:NR);

    % MBF: MiLAC (LMMSE)
    Q_MBF_LMMSE = inv([[-eye(NR),zeros(NR,NT)];[H',eye(NT)]]);
    Q21_MBF_LMMSE = sqrt(Pt) * Q_MBF_LMMSE(NR+1:end,1:NR) / norm(Q_MBF_LMMSE(NR+1:end,1:NR),'fro');

    % ZFBF: Digital
    F_ZFBF = H'/(H*H');
    W_ZFBF = sqrt(Pt/NR) * F_ZFBF ./ vecnorm(F_ZFBF);
    %W_ZFBF = sqrt(Pt) * F_ZFBF / norm(F_ZFBF,'fro');

    % ZFBF: MiLAC (arbit.)
    Q_ZFBF_arbit = inv([[eye(NR),zeros(NR,NT)];[-W_ZFBF,eye(NT)]]);
    Q21_ZFBF_arbit = Q_ZFBF_arbit(NR+1:end,1:NR);

    % ZFBF: MiLAC (LMMSE)
    Q_ZFBF_LMMSE = inv([[zeros(NR),H];[H',eye(NT)]]);
    Q21_ZFBF_LMMSE = sqrt(Pt) * Q_ZFBF_LMMSE(NR+1:end,1:NR) / norm(Q_ZFBF_LMMSE(NR+1:end,1:NR),'fro');

    for iSNR = 1:length(SNR_dBs)

        % R-ZFBF: Digital
        F_RZFBF = H'/(H*H'+NR/db2pow(SNR_dBs(iSNR))*eye(NR));
        W_RZFBF = sqrt(Pt/NR) * F_RZFBF ./ vecnorm(F_RZFBF);
        %W_RZFBF = sqrt(Pt) * F_RZFBF / norm(F_RZFBF,'fro');

        % R-ZFBF: MiLAC (arbit.)
        Q_RZFBF_arbit = inv([[eye(NR),zeros(NR,NT)];[-W_RZFBF,eye(NT)]]);
        Q21_RZFBF_arbit = Q_RZFBF_arbit(NR+1:end,1:NR);

        % R-ZFBF: MiLAC (LMMSE)
        Q_RZFBF_LMMSE = inv([[-NR/db2pow(SNR_dBs(iSNR))*eye(NR),H];[H',eye(NT)]]);
        Q21_RZFBF_LMMSE = sqrt(Pt) * Q_RZFBF_LMMSE(NR+1:end,1:NR) / norm(Q_RZFBF_LMMSE(NR+1:end,1:NR),'fro');
        
        % Sum Rate
        SR_MBF(iMonte,iSNR,1) = func_sum_rate(W_MBF,H,s2(iSNR));
        SR_MBF(iMonte,iSNR,2) = func_sum_rate(Q21_MBF_arbit,H,s2(iSNR));
        SR_MBF(iMonte,iSNR,3) = func_sum_rate(Q21_MBF_LMMSE,H,s2(iSNR));

        SR_ZFBF(iMonte,iSNR,1) = func_sum_rate(W_ZFBF,H,s2(iSNR));
        SR_ZFBF(iMonte,iSNR,2) = func_sum_rate(Q21_ZFBF_arbit,H,s2(iSNR));
        SR_ZFBF(iMonte,iSNR,3) = func_sum_rate(Q21_ZFBF_LMMSE,H,s2(iSNR));

        SR_RZFBF(iMonte,iSNR,1) = func_sum_rate(W_RZFBF,H,s2(iSNR));
        SR_RZFBF(iMonte,iSNR,2) = func_sum_rate(Q21_RZFBF_arbit,H,s2(iSNR));
        SR_RZFBF(iMonte,iSNR,3) = func_sum_rate(Q21_RZFBF_LMMSE,H,s2(iSNR));

        %SR_ZFBF_teor(iMonte,iSNR) = func_sum_rate_ZFBF(H,db2pow(SNR_dBs(iSNR)));
    end
end

toc;

SR_MBF_mean = squeeze(mean(SR_MBF));
SR_ZFBF_mean = squeeze(mean(SR_ZFBF));
SR_RZFBF_mean = squeeze(mean(SR_RZFBF));

SR_ZFBF_teor_mean = mean(SR_ZFBF_teor);

%% Plot
figure('DefaultAxesFontSize',12);
LineW = 1.8; MarkS = 8;
plot(SNR_dBs,SR_RZFBF_mean(:,1),'-h','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Digital');
hold on;
plot(SNR_dBs,SR_ZFBF_mean(:,1),'-s','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Digital');
plot(SNR_dBs,SR_MBF_mean(:,1),'-v','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Digital');

plot(SNR_dBs,SR_RZFBF_mean(:,2),'--p','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC (arbit.)');
plot(SNR_dBs,SR_ZFBF_mean(:,2),'--o','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC (arbit.)');
plot(SNR_dBs,SR_MBF_mean(:,2),'--^','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC (arbit.)');

plot(SNR_dBs,SR_RZFBF_mean(:,3),':p','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC (LMMSE), R-ZFBF');
plot(SNR_dBs,SR_ZFBF_mean(:,3),':o','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC (LMMSE), ZFBF');
plot(SNR_dBs,SR_MBF_mean(:,3),':^','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC (LMMSE), MBF');

%plot(SNR_dBs,SR_ZFBF_teor_mean,':','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','ZFBF teor');
grid on;
set(gca,'GridLineStyle',':','GridAlpha',0.5,'LineWidth',1.2);
xlabel('SNR [dB]');
ylabel('Sum rate [bps/Hz]');
legend('Location','northwest','NumColumns',3,'FontSize',9);

set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1.5);