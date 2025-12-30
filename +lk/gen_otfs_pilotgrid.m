function [txDDGrid, pilotMap] = gen_otfs_pilotgrid(sysPar, carrier, pilotCfg)
%gen_otfs_pilotgrid Generate a minimal OTFS pilot grid for MIMO estimation.
%
% Description:
%   This function allocates one pilot per TX antenna in the delay-Doppler
%   grid. Pilot locations are unique across all Tx ports (and Trs). Pilot
%   values follow OTFS_simulation-style embedded pilot settings.
%   Output:
%     txDDGrid: nDelay * nDoppler * nTx * nTr
%     pilotMap: struct with pilot delay/doppler indices per (iTx,iTr)
%
% Developer: Jia. Institution: PML. Date: 2024/06/02

if nargin < 3 || isempty(pilotCfg)
    if isfield(sysPar, 'otfsPilot')
        pilotCfg = sysPar.otfsPilot;
    else
        pilotCfg = struct();
    end
end

nTx = sysPar.nTx;
nTr = sysPar.nTr;
nDelay = carrier.NSizeGrid * 12;
nDoppler = carrier.SymbolsPerSlot * sysPar.nRSslot;
pilotCfg = apply_pilot_defaults(pilotCfg, nDoppler);

totalPilots = nTx * nTr;
if totalPilots > nDelay * nDoppler
    error('OTFS pilot grid is too small for %d pilot ports.', totalPilots);
end

txDDGrid = zeros(nDelay, nDoppler, nTx, nTr);
pilotMap = struct('delay', [], 'doppler', []);
pilotMap.delay = zeros(nTx, nTr);
pilotMap.doppler = zeros(nTx, nTr);
pilotMap.form = pilotCfg.form;
pilotMap.amp = pilotCfg.embeddedPilotAmp;
pilotMap.qamMode = pilotCfg.qamMode;

for iTr = 1 : nTr
    for iTx = 1 : nTx
        pilotIndex = (iTr - 1) * nTx + iTx - 1;
        delayIdx = mod(pilotIndex, nDelay) + 1;
        dopplerIdx = floor(pilotIndex / nDelay) + 1;
        txDDGrid(delayIdx, dopplerIdx, iTx, iTr) = ...
            gen_pilot_value(pilotCfg.form, pilotCfg.qamMode, pilotCfg.embeddedPilotAmp);
        pilotMap.delay(iTx, iTr) = delayIdx;
        pilotMap.doppler(iTx, iTr) = dopplerIdx;
    end
end
end

function pilotCfg = apply_pilot_defaults(pilotCfg, nDoppler)
if ~isfield(pilotCfg, 'form') || isempty(pilotCfg.form)
    pilotCfg.form = 'OTFS_delta';
end
if ~isfield(pilotCfg, 'qamMode') || isempty(pilotCfg.qamMode)
    pilotCfg.qamMode = 4;
end
if ~isfield(pilotCfg, 'embeddedPilotAmp') || isempty(pilotCfg.embeddedPilotAmp)
    pilotCfg.embeddedPilotAmp = nDoppler;
end
end

function pilotValue = gen_pilot_value(pilotForm, qamMode, embeddedPilotAmp)
switch pilotForm
    case 'OTFS_delta'
        pilotValue = qammod(qamMode - 1, qamMode, 'bin') * embeddedPilotAmp;
    case 'block'
        pilotValue = qammod(randi([0 qamMode - 1], 1, 1), qamMode);
    case 'comb_like'
        pilotValue = qammod(randi([0 qamMode - 1], 1, 1), qamMode);
    otherwise
        error(['Unsupported OTFS pilot form: ', pilotForm]);
end
end
