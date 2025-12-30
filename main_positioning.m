%--------------------------------------------------------------------------
% 5G localization link level simulator for 2-D positioning
%--------------------------------------------------------------------------
close all;
clear all;
clc;
rng('default');
data.LocErrall = [];
%
% please set ULA array at BS for this simulation, as 1D-angle-based positioning considered.
%
[sysPar, carrier, BeamSweep, Layout, RFI, PE, Chan] = lk.gen_sysconfig_pos;
sysPar.Modulation = 'OFDM';
% generate wireless channel information
[data.CIR_cell, data.hcoef, data.Hinfo] = lk.gen_channelcoeff(sysPar, ...
    carrier, Layout, Chan, RFI);
% generate RS symbols
[data.rsSymbols, data.rsIndices, data.txGrid] = lk.gen_rssymbol(sysPar, ...
    carrier,BeamSweep.IndBmSweep);
% OFDM modulation
[data.txWaveform] = lk.gen_transmitsignal(sysPar, carrier,data, RFI, BeamSweep);
% Channel filtering
[data.rxWaveform] = lk.gen_receivesignal(sysPar, carrier, data, RFI, BeamSweep);
% OFDM demodulation
[data.rxGrid] = lk.gen_demodulatedgrid(sysPar, carrier, data.rxWaveform);
% Channel estimation % works for SRS or CSIRS
[data.hcfr_esti] = lk.gen_estimated_cfr(sysPar, carrier, data);
data.cfr_true = lk.gen_cfr_true(sysPar, carrier, data.CIR_cell);
data.nmse_ofdm = lk.gen_nmse(data.hcfr_esti, data.cfr_true);
% %     angle estimation
[data.Angle_esti] = lk.gen_estimated_angle(sysPar, data.hcfr_esti,PE);
[data.Range_esti] = lk.gen_estimatedTOA(sysPar,data.hcfr_esti,PE);
%     % LS localization (multi-BS only)
if sysPar.nBS >= 2
    [data.EstiLoc, data.LocErr, data.BS_sel] = lk.gen_UElocation(sysPar, data,PE);
    data.LocErrall = cat(2, data.LocErrall, data.LocErr);
end
pf.plotSysLayout(sysPar,Layout, data);  % display Layout
pf.plotCIRMeas(sysPar, data);
if sysPar.nBS >= 2
    pf.plotCDF(sysPar, data);% display CDF
end
pf.plotPDP(sysPar, data);% display power-delay profile
pf.plotCIR(sysPar, data); % display perfect CIR
pf.plotPhaseDiff(sysPar, data); % display Phase difference
pf.plotChEstiCFR(sysPar, data); % display estimated CFR
pf.plotResources(sysPar, data); % display RS resources grid
% OTFS baseline (pilot-based MIMO estimation)
sysPar.Modulation = 'OTFS';
data_otfs = struct();
data_otfs.CIR_cell = data.CIR_cell;
data_otfs.hcoef = data.hcoef;
data_otfs.Hinfo = data.Hinfo;
[data_otfs.txGrid, data_otfs.pilotMap] = lk.gen_otfs_pilotgrid(sysPar, carrier);
[data_otfs.txWaveform] = lk.gen_transmitsignal(sysPar, carrier, data_otfs, RFI, BeamSweep);
[data_otfs.rxWaveform] = lk.gen_receivesignal(sysPar, carrier, data_otfs, RFI, BeamSweep);
[data_otfs.rxGrid, data_otfs.rxDDGrid] = lk.gen_demodulatedgrid(...
    sysPar, carrier, data_otfs.rxWaveform);
data_otfs.hdd_esti = lk.gen_otfs_mimoestimate(sysPar, carrier, data_otfs.txGrid, ...
    data_otfs.rxDDGrid);
data_otfs.htf_esti = lk.gen_otfs_tfestimate(data_otfs.hdd_esti);
symPerSlot = carrier.SymbolsPerSlot;
nRSslot = sysPar.nRSslot;
nSC = carrier.NSizeGrid * 12;
htf_slot = zeros(nSC, nRSslot, sysPar.nRx, sysPar.nTx, sysPar.nRr, sysPar.nTr);
for islot = 1 : nRSslot
    idxSym = (islot - 1) * symPerSlot + 1;
    htf_slot(:, islot, :, :, :, :) = data_otfs.htf_esti(:, idxSym, :, :, :, :);
end
data_otfs.nmse_otfs = lk.gen_nmse(permute(htf_slot, [1 3 4 2 5 6]), data.cfr_true);
disp(['OFDM NMSE (TF): ', num2str(data.nmse_ofdm)]);
disp(['OTFS NMSE (TF): ', num2str(data_otfs.nmse_otfs)]);
%%
