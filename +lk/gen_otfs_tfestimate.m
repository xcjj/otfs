function htf_esti = gen_otfs_tfestimate(hdd_esti)
%gen_otfs_tfestimate Convert OTFS DD estimate to TF estimate (SFFT).
%
% Description:
%   Applies SFFT over delay and Doppler dimensions for each antenna pair.
%   Input: hdd_esti: nDelay * nDoppler * nRx * nTx * nRr * nTr
%   Output: htf_esti: nDelay * nDoppler * nRx * nTx * nRr * nTr
%
% Developer: Jia. Institution: PML. Date: 2024/06/02

[nDelay, nDoppler, nRx, nTx, nRr, nTr] = size(hdd_esti);
htf_esti = zeros(size(hdd_esti));
scale = sqrt(nDelay * nDoppler);

for iTr = 1 : nTr
    for iRr = 1 : nRr
        for iTx = 1 : nTx
            for iRx = 1 : nRx
                temp = ifft(hdd_esti(:, :, iRx, iTx, iRr, iTr), [], 1);
                temp = fft(temp, [], 2);
                htf_esti(:, :, iRx, iTx, iRr, iTr) = temp * scale;
            end
        end
    end
end
end
