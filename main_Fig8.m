clear; clc;
rng(3);
tic;

% Parameters
NTNR = [[256,1024,1024,4096,4096,4096];
        [256, 256,1024, 256,1024,4096]];
Pt = 1;
SNR_dB = 10; % 0, 10, or 20
SNR = db2pow(SNR_dB);
s2 = Pt ./ SNR;
nMonte = 2;

tau = 1e2;

% Complexity
C_digital_RZFBF = nan(1,length(NTNR));
C_digital_MBF = nan(1,length(NTNR));
for iNTNR = 1:length(NTNR)
    NT = NTNR(1,iNTNR);
    NR = NTNR(2,iNTNR);
    C_digital_RZFBF(iNTNR) = 8*(NT*NR^2+NR^3/3) + 8*NT*NR*tau;
    C_digital_MBF(iNTNR) = 8*NT*NR*tau;
end
C_MiLAC_RZFBF = 6*prod(NTNR);

% Sum rate
SR_digital_RZFBF = nan(nMonte,length(NTNR));
SR_MiLAC_RZFBF = nan(nMonte,length(NTNR));
for iMonte = 1:nMonte
    for iNTNR = 1:length(NTNR)
        NT = NTNR(1,iNTNR);
        NR = NTNR(2,iNTNR);

        H = sqrt(1/2) * (randn(NR,NT) + 1i * randn(NR,NT));

        % R-ZFBF
        F_RZFBF = H'/(H*H'+NR/SNR*eye(NR));
        W_digital_RZFBF = sqrt(Pt/NR) * F_RZFBF ./ vecnorm(F_RZFBF);
        W_MiLAC_RZFBF = sqrt(Pt) * F_RZFBF / norm(F_RZFBF,'fro');

        % Sum Rate
        SR_digital_RZFBF(iMonte,iNTNR) = func_sum_rate(W_digital_RZFBF,H,s2);
        SR_MiLAC_RZFBF(iMonte,iNTNR) = func_sum_rate(W_MiLAC_RZFBF,H,s2);
    end
end

toc;

SR_digital_RZFBF_mean = mean(SR_digital_RZFBF);
SR_MiLAC_RZFBF_mean = mean(SR_MiLAC_RZFBF);

%% Plot
figure('DefaultAxesFontSize',12);
LineW = 1.8; MarkS = 8;
hold on;
semilogx(C_MiLAC_RZFBF(6),SR_MiLAC_RZFBF_mean(6),'h','color','#D95319','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','#D95319','HandleVisibility','off');
semilogx(C_MiLAC_RZFBF(5),SR_MiLAC_RZFBF_mean(5),'p','color','#D95319','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','#D95319','HandleVisibility','off');
semilogx(C_MiLAC_RZFBF(4),SR_MiLAC_RZFBF_mean(4),'^','color','#D95319','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','#D95319','HandleVisibility','off');
semilogx(C_MiLAC_RZFBF(3),SR_MiLAC_RZFBF_mean(3),'o','color','#D95319','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','#D95319','HandleVisibility','off');
semilogx(C_MiLAC_RZFBF(2),SR_MiLAC_RZFBF_mean(2),'x','color','#D95319','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','#D95319','HandleVisibility','off');
semilogx(C_MiLAC_RZFBF(1),SR_MiLAC_RZFBF_mean(1),'d','color','#D95319','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','#D95319','HandleVisibility','off');

semilogx(C_digital_RZFBF(6),SR_digital_RZFBF_mean(6),'h','color','#0072BD','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','#0072BD','HandleVisibility','off');
semilogx(C_digital_RZFBF(5),SR_digital_RZFBF_mean(5),'p','color','#0072BD','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','#0072BD','HandleVisibility','off');
semilogx(C_digital_RZFBF(4),SR_digital_RZFBF_mean(4),'^','color','#0072BD','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','#0072BD','HandleVisibility','off');
semilogx(C_digital_RZFBF(3),SR_digital_RZFBF_mean(3),'o','color','#0072BD','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','#0072BD','HandleVisibility','off');
semilogx(C_digital_RZFBF(2),SR_digital_RZFBF_mean(2),'x','color','#0072BD','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','#0072BD','HandleVisibility','off');
semilogx(C_digital_RZFBF(1),SR_digital_RZFBF_mean(1),'d','color','#0072BD','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','#0072BD','HandleVisibility','off');

semilogx(nan,nan,'hk','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','k','DisplayName','4096 x 4096');
semilogx(nan,nan,'pk','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','k','DisplayName','1024 x 4096');
semilogx(nan,nan,'^k','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','k','DisplayName','  256 x 4096');
semilogx(nan,nan,'ok','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','k','DisplayName','1024 x 1024');
semilogx(nan,nan,'xk','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','k','DisplayName','  256 x 1024');
semilogx(nan,nan,'dk','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','k','DisplayName','  256 x 256');
semilogx(nan,nan,'s','color','#0072BD','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','#0072BD','DisplayName','Digital');
semilogx(nan,nan,'s','color','#D95319','linewidth',LineW,'MarkerSize',MarkS,'MarkerFaceColor','#D95319','DisplayName','MiLAC');
grid on; box on;
set(gca,'xscale','log')
set(gca,'GridLineStyle',':','GridAlpha',0.5,'LineWidth',1.2);
xlabel('Computational complexity');
ylabel('Sum rate [bps/Hz]');
legend('Location','southeast','NumColumns',1,'FontSize',11);
ax = gca;
ax.XTick = 10.^(5:1:13);
ax.XLim = [1e5 1e13];
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1.5);