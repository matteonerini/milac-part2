clear; clc;
rng(3);
tic;

% Parameters
NR = 4;
NT = 4;
Pt = 1;
SNR_dBs = -10:5:20;
s2 = Pt ./ db2pow(SNR_dBs);
nBits = 1e5;

bits = round(rand(nBits,1));
qpskSymbols = func_qpsk_mod(bits, 1);
s = reshape(qpskSymbols,NT,[]);

berMF = zeros(length(SNR_dBs),2);
berZF = zeros(length(SNR_dBs),2);
berMMSE = zeros(length(SNR_dBs),2);

for iSNR=1:length(SNR_dBs)
    
    y = zeros(NR,length(s));

    bitsMF = zeros(nBits,2);
    bitsZF = zeros(nBits,2);
    bitsMMSE = zeros(nBits,2);
    bitsML = zeros(nBits,1);

    for iSymbol = 1:length(s)
        H = sqrt(1 / 2) * (randn(NR, NT) + 1i * randn(NR, NT));
        n = sqrt(s2(iSNR) / 2) * (randn(NR,1) + 1i * randn(NR,1));

        % Spatial Multiplexing
        y(:,iSymbol) = sqrt(Pt/NT) * H * s(:,iSymbol) + n;
        
        % MF Receiver
        zMF_DB = sqrt(NT/Pt) * (db2pow(SNR_dBs(iSNR))/NT) * H' * y(:,iSymbol);
        bitsMF(2*NT*iSymbol-2*NT+1 : 2*NT*iSymbol,1) = func_qpsk_demod(zMF_DB);

        Q_MF = inv([[eye(NR),zeros(NR,NT)];[H',-NT/db2pow(SNR_dBs(iSNR))*eye(NT)]]);
        zMF_RF = Q_MF(NR+1:end,1:NR) * y(:,iSymbol);
        bitsMF(2*NT*iSymbol-2*NT+1 : 2*NT*iSymbol,2) = func_qpsk_demod(zMF_RF);

        % ZF Receiver
        zZF_DB = sqrt(NT/Pt) * ((H'*H) \ H') * y(:,iSymbol);
        bitsZF(2*NT*iSymbol-2*NT+1 : 2*NT*iSymbol,1) = func_qpsk_demod(zZF_DB);

        Q_ZF = inv([[eye(NR),H];[H',zeros(NT)]]);
        zZF_RF = Q_ZF(NR+1:end,1:NR) * y(:,iSymbol);
        bitsZF(2*NT*iSymbol-2*NT+1 : 2*NT*iSymbol,2) = func_qpsk_demod(zZF_RF);

        % MMSE Receiver
        zMMSE_DB = sqrt(NT/Pt) * ((H'*H + NT/db2pow(SNR_dBs(iSNR))*eye(NT)) \ H') * y(:,iSymbol);
        bitsMMSE(2*NT*iSymbol-2*NT+1 : 2*NT*iSymbol,1) = func_qpsk_demod(zMMSE_DB);

        Q_RZFBF = inv([[eye(NR),H];[H',-NT/db2pow(SNR_dBs(iSNR))*eye(NT)]]);
        zMMSE_RF = Q_RZFBF(NR+1:end,1:NR) * y(:,iSymbol);
        bitsMMSE(2*NT*iSymbol-2*NT+1 : 2*NT*iSymbol,2) = func_qpsk_demod(zMMSE_RF);

    end
    berMF(iSNR,:) = sum(xor([bits,bits], bitsMF)) / nBits;
    berZF(iSNR,:) = sum(xor([bits,bits], bitsZF)) / nBits;
    berMMSE(iSNR,:) = sum(xor([bits,bits], bitsMMSE)) / nBits;
end

toc;

%% Plot
figure('DefaultAxesFontSize',12);
LineW = 1.8; MarkS = 8;
semilogy(SNR_dBs,berMMSE(:,1),'-h','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Digital');
hold on;
semilogy(SNR_dBs,berZF(:,1),'-s','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Digital');
semilogy(SNR_dBs,berMF(:,1),'-v','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Digital');

semilogy(SNR_dBs,berMMSE(:,2),'--p','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC, MMSE');
semilogy(SNR_dBs,berZF(:,2),'--o','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC, ZF');
semilogy(SNR_dBs,berMF(:,2),'--^','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','MiLAC, MF');

grid on;
set(gca,'GridLineStyle',':','GridAlpha',0.5,'LineWidth',1.2);
xlabel('SNR [dB]');
ylabel('BER');
legend('Location','southwest','NumColumns',2,'FontSize',12);
ylim([1e-3 1])
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1.5);