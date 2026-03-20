clear; clc;

% Parameters
Ns = [256,512,1024,1024:1024:1024*8];
tau = 1e2;

% Digital beamforming
DB_LS = 8*(Ns.^3+Ns.^3/3) + 8*Ns.^2 * tau; % Matrix-matrix product: 2N^3 complex ops. but only half is needed since hermitian
DB_MF = 8*Ns.^2 * tau; % Matrix-vector product: 2N^2 complex ops.
DB_DFT = 34/9*Ns.*log2(Ns) * tau; % https://ieeexplore.ieee.org/document/4034175

% Generalized analog beamforming
AB_LS = 6*Ns.^2; % Set 2N diagonal entries, each summing N values: 2N^2 complex adds. plus multiply by Y_0
AB_MF = 4*Ns.^2; % Set N diagonal entries, each summing N values: N^2 complex adds. plus multiply by Y_0
AB_DFT = ones(size(Ns));

%% Plot
figure('defaultaxesfontsize',12)
LineW = 1.8; MarkS = 8;
semilogy(Ns,DB_LS,'-h','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Digital')
hold on;
semilogy(Ns,DB_MF,'-s','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Digital')
semilogy(Ns,DB_DFT,'-v','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Digital')
set(gca,'ColorOrderIndex',1);
semilogy(Ns,AB_LS,'--p','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC, Zero-forcing')
semilogy(Ns,AB_MF,'--o','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC, Matched filtering')
semilogy(Ns,AB_DFT,'--^','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC, DFT')
grid on;
set(gca,'GridLineStyle',':','GridAlpha',0.5,'LineWidth',1.2);
xlabel('Number of antennas');
ylabel('Computational complexity');
legend('location','southeast','NumColumns',2,'FontSize',12);
ax = gca;
ax.XTick = 0:1024:1024*8;
ax.XLim = [0 max(Ns)];
ax.YTick = 10.^(0:2:13);
ax.YLim = [1e0 1e13];
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1.5);