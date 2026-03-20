function quants = func_quant_Gauss(sig, B, sigma)
    % https://ieeexplore.ieee.org/document/1057548

    if B ~= Inf

        codebooks = cell(1,4);
        codebooks{1} = [0.7980];
        codebooks{2} = [0.4528,1.5100];
        codebooks{3} = [0.2451,0.7560,1.3440,2.1520];
        codebooks{4} = [0.1284,0.3881,0.6568,0.9424,1.2560,1.6180,2.0690,2.7330];
    
        codebook = sigma * [-flip(codebooks{B}), codebooks{B}];
    
        partition = (codebook(2:end) + codebook(1:end-1)) / 2;
    
        [~,quants] = quantiz(sig,partition,codebook);

    else

        quants = sig;

    end

end