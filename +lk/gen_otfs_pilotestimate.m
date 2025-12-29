function hdd_esti = gen_otfs_pilotestimate(txDDGrid, rxDDGrid)
%gen_otfs_pilotestimate Minimal OTFS channel estimate using pilot division.
%
% Description:
%   This function provides a baseline OTFS estimate by dividing the received
%   delay-Doppler grid by the transmitted grid at pilot positions.
%   Note: This is a SISO baseline using the first Tx/Rx antenna only.
%
% Developer: Jia. Institution: PML. Date: 2024/06/02

txPilot = txDDGrid(:,:,1,1);
rxPilot = rxDDGrid(:,:,1,1);
mask = abs(txPilot) > 0;
hdd_esti = zeros(size(rxPilot));
hdd_esti(mask) = rxPilot(mask) ./ txPilot(mask);
end
