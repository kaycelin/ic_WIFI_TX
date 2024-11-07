%% 2024-05-29, draft
function [lut_matrix, lut] = paLookupTable_process(filename,sheetname,pinRange,fnum,isFitting)
    % filename = 'paData.xlsx';
    % sheetname = 'MEA2GCW';
    % pinRange = [-10 3];
    % fnum = 0529;
    
    %%
    if exist('pinRange','var')&&~isempty(pinRange)
        isExcludeNoise = 1;
    else
        isExcludeNoise = 0;
    end
    if exist('fnum','var')&&~isempty(fnum)
        isFnum = 1;
    else
        isFnum = 0;
    end
    if ~exist('isFitting','var')||isempty(isFitting)
        isFitting = 1;
    end
    % get excel
    lut = readtable(filename,'Sheet',sheetname);
    if isExcludeNoise
        [~, idx1] = min(abs(lut.PinDbm-pinRange(1)));
        [~, idxEnd] = min(abs(lut.PinDbm-pinRange(end)));
        lut = lut(idx1:idxEnd,:);
    end
    pin = lut.PinDbm;
    pout = lut.PoutDbm;
    phsShift = lut.PhaseShiftDeg;
    phsShift = phsShift - phsShift(1);
    gain = pout - pin;
    gainCompress = gain - gain(1);
    if isFitting
        p_ps = polyfit(pin,phsShift,7);
        phsShift = polyval(p_ps,pin);
        phsShift = phsShift - phsShift(1);

        p_gc = polyfit(pin,gainCompress,7);
        gainCompress = polyval(p_gc,pin);
        gainCompress = gainCompress - gainCompress(1);
        pout = pin + gainCompress + gain(1);

        % export to LUT
        lut.GC = gainCompress;
        lut.GainDb = gainCompress + gain(1);
        lut.PoutDbm = pout;
        lut.PhaseShiftDeg = phsShift;
    end

    if isFnum
        sheetname_2 = strrep(sheetname,'_','.');
        figure(fnum)
        subplot(1,3,1), plot(pin,pout,'DisplayName',sheetname_2), xlabel('Pin [dBm]'), ylabel('Pout [dBm]'), title('AMAM'), hold on; legend, grid minor
        subplot(1,3,2), plot(pin,gainCompress,'DisplayName',sheetname_2), xlabel('Pin [dBm]'), ylabel('Gain Compression [dB]'), title('Gain Compression'), hold on; legend, grid minor
        subplot(1,3,3), plot(pin,phsShift,'DisplayName',sheetname_2), xlabel('Pin [dBm]'), ylabel('Phase Shift [deg]'), title('AMPM'), hold on; legend, grid minor
    end

    lut_matrix = [pin, pout, phsShift];

end