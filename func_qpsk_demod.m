function [bits, SNR] = func_qpsk_demod(symbols)
    % ML detection for QPSK symbols
    % real > 0, imag > 0 -> close to +sqrt(E_s) + 1i * sqrt(E_s) -> [0, 0]
    % real > 0, imag < 0 -> close to +sqrt(E_s) - 1i * sqrt(E_s) -> [0, 1]
    % real < 0, imag > 0 -> close to -sqrt(E_s) + 1i * sqrt(E_s) -> [1, 0]
    % real < 0, imag < 0 -> close to -sqrt(E_s) - 1i * sqrt(E_s) -> [1, 1]
    
    bits(1:2:2 * length(symbols)-1) = 1 / 2 * (1 - sign(real(symbols)));
    bits(2:2:2 * length(symbols))   = 1 / 2 * (1 - sign(imag(symbols)));
    SNR = norm(symbols) ^ 2 / length(symbols);
end