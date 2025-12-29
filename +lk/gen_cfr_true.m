function cfr_true = gen_cfr_true(sysPar, carrier, CIR_cell)
%gen_cfr_true Generate true CFR from channel impulse response.
%
% Description:
%   Converts CIR to CFR using FFT and the same subcarrier mapping used in
%   OFDM demodulation.
%   Output: cfr_true: nSC * nRx * nTx * nRSslot * nRr * nTr
%
% Developer: Jia. Institution: PML. Date: 2024/06/02

nfft = carrier.Nfft;
nSC = carrier.NSizeGrid * 12;
nRr = sysPar.nRr;
nTr = sysPar.nTr;
nRSslot = sysPar.nRSslot;
nTx = sysPar.nTx;
nRx = sysPar.nRx;

cfr_true = zeros(nSC, nRx, nTx, nRSslot, nRr, nTr);

for iRr = 1 : nRr
    for iTr = 1 : nTr
        CIR = CIR_cell{iRr, iTr};
        for islot = 1 : nRSslot
            cir_slot = CIR(:, :, :, islot);
            hf = fft(cir_slot, nfft, 1);
            hf = circshift(hf, nSC / 2, 1);
            hf = hf(1:nSC, :, :);
            cfr_true(:, :, :, islot, iRr, iTr) = hf;
        end
    end
end
end
