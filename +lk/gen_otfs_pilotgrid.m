function [txDDGrid, pilotMap] = gen_otfs_pilotgrid(sysPar, carrier)
%gen_otfs_pilotgrid Generate a minimal OTFS pilot grid for MIMO estimation.
%
% Description:
%   This function allocates one pilot per TX antenna in the delay-Doppler
%   grid. Pilot locations are unique across all Tx ports (and Trs).
%   Output:
%     txDDGrid: nDelay * nDoppler * nTx * nTr
%     pilotMap: struct with pilot delay/doppler indices per (iTx,iTr)
%
% Developer: Jia. Institution: PML. Date: 2024/06/02

nTx = sysPar.nTx;
nTr = sysPar.nTr;
nDelay = carrier.NSizeGrid * 12;
nDoppler = carrier.SymbolsPerSlot * sysPar.nRSslot;

totalPilots = nTx * nTr;
if totalPilots > nDelay * nDoppler
    error('OTFS pilot grid is too small for %d pilot ports.', totalPilots);
end

txDDGrid = zeros(nDelay, nDoppler, nTx, nTr);
pilotMap = struct('delay', [], 'doppler', []);
pilotMap.delay = zeros(nTx, nTr);
pilotMap.doppler = zeros(nTx, nTr);

for iTr = 1 : nTr
    for iTx = 1 : nTx
        pilotIndex = (iTr - 1) * nTx + iTx - 1;
        delayIdx = mod(pilotIndex, nDelay) + 1;
        dopplerIdx = floor(pilotIndex / nDelay) + 1;
        txDDGrid(delayIdx, dopplerIdx, iTx, iTr) = 1;
        pilotMap.delay(iTx, iTr) = delayIdx;
        pilotMap.doppler(iTx, iTr) = dopplerIdx;
    end
end
end
