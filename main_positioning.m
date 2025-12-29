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
% OTFS baseline (minimal pilot-based estimation)
sysPar.Modulation = 'OTFS';
data_otfs = struct();
data_otfs.CIR_cell = data.CIR_cell;
data_otfs.hcoef = data.hcoef;
data_otfs.Hinfo = data.Hinfo;
[data_otfs.rsSymbols, data_otfs.rsIndices, data_otfs.txGrid] = lk.gen_rssymbol(...
    sysPar, carrier, BeamSweep.IndBmSweep);
[data_otfs.txWaveform] = lk.gen_transmitsignal(sysPar, carrier, data_otfs, RFI, BeamSweep);
[data_otfs.rxWaveform] = lk.gen_receivesignal(sysPar, carrier, data_otfs, RFI, BeamSweep);
[data_otfs.rxGrid, data_otfs.rxDDGrid] = lk.gen_demodulatedgrid(...
    sysPar, carrier, data_otfs.rxWaveform);
data_otfs.hdd_esti = lk.gen_otfs_pilotestimate(data_otfs.txGrid, data_otfs.rxDDGrid);
%%
