function hdd_esti = gen_otfs_mimoestimate(sysPar, carrier, txDDGrid, rxDDGrid)
%gen_otfs_mimoestimate Minimal OTFS MIMO channel estimate from pilots.
%
% Description:
%   Estimates the delay-Doppler channel at pilot locations by dividing the
%   received DD grid by the transmitted pilot grid per (Tx,Tr).
%   Output: hdd_esti: nDelay * nDoppler * nRx * nTx * nRr * nTr
%
% Developer: Jia. Institution: PML. Date: 2024/06/02

nTx = sysPar.nTx;
nRx = sysPar.nRx;
nTr = sysPar.nTr;
nRr = sysPar.nRr;
nDelay = carrier.NSizeGrid * 12;
nDoppler = carrier.SymbolsPerSlot * sysPar.nRSslot;

hdd_esti = zeros(nDelay, nDoppler, nRx, nTx, nRr, nTr);

for iTr = 1 : nTr
    for iTx = 1 : nTx
        pilotIndex = (iTr - 1) * nTx + iTx - 1;
        delayIdx = mod(pilotIndex, nDelay) + 1;
        dopplerIdx = floor(pilotIndex / nDelay) + 1;
        pilotVal = txDDGrid(delayIdx, dopplerIdx, iTx, iTr);
        if pilotVal == 0
            continue;
        end
        for iRr = 1 : nRr
            for iRx = 1 : nRx
                rxVal = rxDDGrid(delayIdx, dopplerIdx, iRx, iRr);
                hdd_esti(delayIdx, dopplerIdx, iRx, iTx, iRr, iTr) = rxVal / pilotVal;
            end
        end
    end
end
end
