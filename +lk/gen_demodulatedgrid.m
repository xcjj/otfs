function [rxGrid, rxDDGrid] = gen_demodulatedgrid(sysPar,carrier,rxWaveform);
%gen_demodulatedgrid Generate receive resource grid (demodulation).
%
% Description:
% This function aims to perform demodulation at the receiver.
% Note: rxGrid: % nSC * nSyms * nRx * nRr
%       rxDDGrid: % nSC * nSyms * nRx * nRr (OTFS only)
%
% Developer: Jia. Institution: PML. Date: 2021/10/28

nRr = sysPar.nRr;
rxGrid = [];
rxDDGrid = [];
for iRr = 1 : nRr
    switch upper(sysPar.Modulation)
        case 'OFDM'
            rxGrid_r = nr.OFDMDemodulate( carrier, rxWaveform(:, :,iRr), ...
                sysPar.center_frequency, 0);
        case 'OTFS'
            [rxDDGrid_r, rxGrid_r] = nr.OTFSDemodulate( carrier, ...
                rxWaveform(:, :,iRr), sysPar.center_frequency, 0);
            rxDDGrid = cat(4, rxDDGrid, rxDDGrid_r);
        otherwise
            error('Unsupported modulation type: %s', sysPar.Modulation);
    end
    rxGrid = cat(4, rxGrid, rxGrid_r);
end
end
