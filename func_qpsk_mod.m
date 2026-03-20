function [symbols] = func_qpsk_mod(bits, Es)
    % [0, 0] -> +sqrt(Es/2) + 1i * sqrt(Es/2)
    % [0, 1] -> +sqrt(Es/2) - 1i * sqrt(Es/2)
    % [1, 0] -> -sqrt(Es/2) + 1i * sqrt(Es/2)
    % [1, 1] -> -sqrt(Es/2) - 1i * sqrt(Es/2)
    
    bits = 1 - 2 * reshape(bits, [2, length(bits) / 2]);
    symbols = sqrt(Es / 2) * (bits(1, :) + 1i * bits(2, :));
end