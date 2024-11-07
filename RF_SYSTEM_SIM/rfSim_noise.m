%% 2024-04-13, draft
%% 2024-04-25, change fucntion name from noiseGen to noise and add input/ouput variable x/y

function [y, noiseOutput] = rfSim_noise(x, noiseMethod, noisePwr_nfDb_tempDeg, fs_nSamps_isCompx_nAnts_cell, Seed_Rohm_cell)
    if exist('x','var')&&~isempty(x)
        isX = 1;
        isCompx = ~isreal(x);
        nSamps = size(x,1);
        nAnts = size(x,2);
        sigPower = sum(abs(x(:)).^2)/numel(x); % linear
        sigPowerDbm = 10*log10(sigPower*1000);
    elseif exist('fs_nSamps_isCompx_nAnts_cell','var')&&~isempty(fs_nSamps_isCompx_nAnts_cell)
        isX = 0;
        nSamps = fs_nSamps_isCompx_nAnts_cell{2};
        isCompx = fs_nSamps_isCompx_nAnts_cell{3};
        nAnts = fs_nSamps_isCompx_nAnts_cell{4};
    else
        error('check input of fs_nSamps_isCompx_nAnts_cell')
    end
    if exist('fs_nSamps_isCompx_nAnts_cell','var')&&~isempty(fs_nSamps_isCompx_nAnts_cell)
        try
            fs = fs_nSamps_isCompx_nAnts_cell{1};
        catch
            fs = fs_nSamps_isCompx_nAnts_cell;
        end
    end
    if ~exist('Seed_Rohm_cell','var')||isempty(Seed_Rohm_cell)
        Rohm = 1;
        Seed = 'shuffle';
    else
        try
            Seed = Seed_Rohm_cell{1};
        catch
            Seed = Seed_Rohm_cell;
        end
        try
            Rohm = Seed_Rohm_cell{2};
        catch
            Rohm = 1;
        end
    end

    kfactor = 1.380649e-23; % Boltzmann constant
    T0 = 290;
    switch noiseMethod
        case 'NoisePowerDbmHz'
            if fs==1
                error('fs should not unit for NoisePowerDbmHz')
            end
            noisePowerDbmHz = noisePwr_nfDb_tempDeg;
            noisePowerDbm = noisePowerDbmHz + 10*log10(fs);
        case 'NoisePowerDbm'
            fs = 1;
            noisePowerDbm = noisePwr_nfDb_tempDeg;
        case 'NoiseTemperature'
            Tequivalent = 273.15 + noisePwr_nfDb_tempDeg; % Converting Celsius to Fahrenhit Temperature
            noisePowerDbm = 10*log10(1000*kfactor*Tequivalent*fs); % calculate noise power from KTB
        case 'NoiseFigure'
            nfDb = noisePwr_nfDb_tempDeg;
            nf = 10^(nfDb/10);
           if nf > 1
                nAddDbm = 10*log10(1000*(nf-1)*kfactor*T0*fs);
            else
                nAddDbm = 10*log10(1000*kfactor*T0*fs);
            end
            noisePowerDbm = nAddDbm;
        case {'SNR'} % Signal-to-noise ratio (SNR) per sample
            if numel(noisePwr_nfDb_tempDeg)==2
                snrDb = noisePwr_nfDb_tempDeg(1);
                sigPowerDbm = noisePwr_nfDb_tempDeg(2);
                sigPower = 10^(sigPowerDbm/10)/1000;
            else
                snrDb = noisePwr_nfDb_tempDeg(1);
                sigPowerDbm;
                sigPower;
            end
        case {'Eb/N0'} % Ratio of bit energy to noise power spectral density (EbN0).
            if isstruct(nsParms)
                numBitsPerSymb = nsParms.nBitsPerSymb;
                EbN0Db = nsParms.EbN0;
                osr = nsParms.osr;
            elseif iscell(nsParms)
                numBitsPerSymb = nsParms{1};
                EbN0Db = nsParms{2};
                osr = nsParms{3};
            end
            snrDb = EbN0Db - 10*log10(numBitsPerSymb) - 10*log10(osr);
        case {'Es/N0'} % Ratio of symbol energy to noise power spectral density (EsN0)
            if isstruct(nsParms)
                EsN0Db = nsParms.EsN0;
                osr = nsParms.osr;
            elseif iscell(nsParms)
                EsN0Db = nsParms{1};
                osr = nsParms{2};
            end
            snrDb = EsN0Db - 10*log10(osr);
    end

    % calculate nosie mangitude
    switch noiseMethod
        case {'NoisePowerDbmHz','NoisePowerDbm','NoiseTemperature'}
            noisePower = 0.001 * 10^((noisePowerDbm)/10) / 1;
        case {'Eb/N0','Es/N0','SNR'}
            noisePower = sigPower/10^(snrDb/10); % calculate noise power from SNR
            noisePowerDbm = 10*log10(noisePower*1000);
    end
    noiseVolt = sqrt(0.001 * 10^((noisePowerDbm)/10) / 1);
    try
        pwrNoiseDbmPsd = 10*log10(noisePower) - 10*log10(fs) +30;
        snrDb_cal =10*log10(sigPower / noisePower);
        pwrSigDbmPsd = 10*log10(sigPower) - 10*log10(fs) +30;
        snrDbPerHz_cal = pwrSigDbmPsd - pwrNoiseDbmPsd;
    end

    noise_total = zeros(nSamps,nAnts);
    kEnd = 20;
    if 0
        rng('shuffle')
    else
        rng(Seed)
    end
    for k=1:kEnd
        noise = zeros(nSamps,nAnts);
        try
            dlSlots = logical(dlSlots);
            nSampsDl = sum(dlSlots);
        catch
            dlSlots = 1:nSamps;
            nSampsDl =nSamps;
        end
        if 0
            scale_dlSlots = sqrt(nSamps/nSampsDl);
        else
            scale_dlSlots = 1; % remove scale
        end
        if ~isCompx
            noise(dlSlots,:) = scale_dlSlots * noiseVolt * randn(nSampsDl,nAnts);
        else
            noise(dlSlots,:) = scale_dlSlots * noiseVolt * (randn(nSampsDl,nAnts) + 1i*randn(nSampsDl,nAnts))/sqrt(2);
        end
        noise_total = noise_total + noise;
    end

    noiseOutput = noise_total/sqrt(kEnd);
    if 0
        powerDbm(noiseOutput)
    end

    if isX
        y = x + noiseOutput;
    else
        y = [];
    end

    if 0
    end
end
