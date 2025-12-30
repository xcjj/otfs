function [ddGrid, tfGrid] = OTFSDemodulate(carrier, rxWaveform, center_frequency, phasecomp)
%OTFSDemodulate Generate the OTFS demodulated grid.
%
% Description:
%   This function performs OFDM demodulation to recover the time-frequency
%   grid, then applies a 2-D SFFT (time-frequency -> delay-Doppler).
%   Output: ddGrid: nSC * nSym * nAnt, tfGrid: nSC * nSym * nAnt
%
% Developer: Jia. Institution: PML. Date: 2024/06/02

tfGrid = nr.OFDMDemodulate(carrier, rxWaveform, center_frequency, phasecomp);
ddGrid = otfs_sfft(tfGrid);
end

function ddGrid = otfs_sfft(tfGrid)
[nDelay, nDoppler, nAnt] = size(tfGrid);
ddGrid = zeros(nDelay, nDoppler, nAnt);
scale = sqrt(nDelay * nDoppler);
for iAnt = 1 : nAnt
    temp = ifft(tfGrid(:,:,iAnt), [], 1);
    temp = fft(temp, [], 2);
    ddGrid(:,:,iAnt) = temp * scale;
end
end
