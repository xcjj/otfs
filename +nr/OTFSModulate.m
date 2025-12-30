function waveform = OTFSModulate(carrier, ddGrid, center_frequency, phasecomp)
%OTFSModulate Generate the OTFS modulated waveform.
%
% Description:
%   This function performs a minimal OTFS modulation by applying a 2-D
%   ISFFT (delay-Doppler -> time-frequency) followed by OFDM modulation.
%   ddGrid: nSC * nSym * nAnt
%
% Developer: Jia. Institution: PML. Date: 2024/06/02

tfGrid = otfs_isfft(ddGrid);
waveform = nr.OFDMModulate(carrier, tfGrid, center_frequency, phasecomp);
end

function tfGrid = otfs_isfft(ddGrid)
[nDelay, nDoppler, nAnt] = size(ddGrid);
tfGrid = zeros(nDelay, nDoppler, nAnt);
scale = sqrt(nDelay * nDoppler);
for iAnt = 1 : nAnt
    temp = fft(ddGrid(:,:,iAnt), [], 1);
    temp = ifft(temp, [], 2);
    tfGrid(:,:,iAnt) = temp / scale;
end
end
