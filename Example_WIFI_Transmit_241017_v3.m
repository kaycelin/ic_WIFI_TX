%% WIFI Transmitter
% History
%% 
% * 2024-09-18, Draft
% * 2024-09-19, Release V1.0
% * 2024-09-25, Release V2.0, Apply CFR - No need for CFR before the PA because 
% the DPD still expands the signal.
% * 2024-10-17, Draft V3.0
% Set WIFI line-up

% close all
% clear all
% clc
clear paInSignal paSignal_org paSignal; % clear signal

Rohm        = 1 % set resistor in ohm
RbwKHz      = 10 % set resolution bw for aclr plot
folderPath  = 'C:\Users\GrantLin\OneDrive - iCana Limited\Documents\MATLAB\waveform';
isPaModel   = 0; % turn on/off pa model

format shortG

simulation_trx = 'tx' % set simulation case as 'tx'
switch simulation_trx
    case {'tx','trx'}
        paModel_LUT_sheet   = '0527SIMCW'  % set LUT sheet
        paModel_methology   = 'RFSIMPA' % set pa model function
        paModel_type        = 'LUT' % set pa's amam/ampm from LUT
        isPaMultiStage      = 1 % set pa as multi-stages
        isPaFlatness        = 1
        isPaFlatnessCompensateLoss = 0
        isMIMO = 0
end

if 1 % set spec. using competitors
    competitorsCase = 'ICWF5519' % iC
    competitorsCase = 'ICWF5971' % iC
    % competitorsCase = 'KCT8570HE' % KCT8570HE
    % competitorsCase = 'SKY85797.11'
    % competitorsCase = 'OTHERS'
end

switch competitorsCase
    case "ICWF5519"
        signalFormat_vec    = ["EHT160","HE160","HE160","VHT80","HT40","HT20"];
        signalFormat_vec    = ["EHT160","HE160","HE160","VHT80","HT40","HT20"]; % test
        MCS_index_vec       = [13 11 11 9 7 0];
        MCS_index_vec       = [13 11 11 9 7 0]; % test
        paOutputPower_target_dBm_vec = [15 16 18 23 24 27] + 0.8;
        % paOutputPower_target_dBm_vec = [15 16 18 23 24 27] + 0.8; % test
        specEVM_vec         = [-47 -43 -43 -35 -30 nan];
    case "ICWF5971"
        signalFormat_vec = ["EHT320","HE160","VHT80","HT40","HT20"];
        MCS_index_vec = [13 11 9 7 0];
        paOutputPower_target_dBm_vec = [15 18 23 24 27] + 0.8;
        % paOutputPower_target_dBm_vec = [13 18 23 24 27] + 0.8;
        specEVM_vec = [-47 -43 -35 -30 nan];
    case "QPF4702"
        signalFormat_vec = ["EHT320","HE160","VHT80","HT40","HT20"];
        MCS_index_vec = [13 11 9 7 0];
        paOutputPower_target_dBm_vec = [18 21 24 24 27] + 0;
        specEVM_vec = [-47 -43 -35 -30 nan];
    case "KCT8570N"
        signalFormat_vec = ["EHT160","HE160","VHT80","HT40","HT20"];
        MCS_index_vec = [13 11 9 7 0];
        paOutputPower_target_dBm_vec = [21.5 22.5 24 25 27] + 0;
        specEVM_vec = [-47 -43 -35 -30 nan];
    otherwise % user-defined
        signalFormat_vec = [...
            "EHT20","EHT40","EHT80","EHT160","EHT320";...
            % "EHT320"
            ];
        MCS_index_vec = [...
            % 0 0 0 0 0;...
            % 1 1 1 1 1;...
            % 2 2 2 2 2;...
            % 3 3 3 3 3;...
            % 4 4 4 4 4;...
            % 5 5 5 5 5;...
            % 6 6 6 6 6;...
            % 7 7 7 7 7;...
            % 8 8 8 8 8;...
            % 9 9 9 9 9;...
            % 10 10 10 10 10;...
            % 11 11 11 11 11;...
            % 12 12 12 12 12;...
            13 13 13 13 13;...
            % 13
            ];
        paOutputPower_target_dBm_vec = [...
            25 24 22.5 21.5 20.5...
            % 20.5
            ];
end

if size(signalFormat_vec,1)~=size(MCS_index_vec,1) % check
    error('!')
end

% Set WIFI Waveform Configurations

if 0
    k1=k;
    kEnd=k
else
    k1=1;
    kEnd = numel(signalFormat_vec)
end

for k=k1:kEnd
    close all

    if 1 % initialization
        clear wv
        signalFormat = signalFormat_vec(k)
        MCS_index = MCS_index_vec(k)
        paOutputPower_target_dBm = paOutputPower_target_dBm_vec(k)
        dispLegd_format = ("MCS"+string(MCS_index)+"."+signalFormat);
        specEVM = specEVM_vec(k)
    end

    % set wv configurations
    configParamsWG                  = struct();
    configParamsWG.wlanPHY          = erase(signalFormat, ["20","40","80","160","320"]); % wv format
    configParamsWG.bw               = erase(signalFormat, ["HT","VHT","EHT","HE"])+"MHz"; % wv bandwidth
    configParamsWG.NumUsers         = 1; % MIMO number of users
    configParamsWG.MCS              = MCS_index; % mcs index
    if isMIMO
        configParamsWG.NumSpaceTimeStreams  = 2; %[1 1];
        configParamsWG.NumTransmitAntennas  = 2; %2;
    else
        configParamsWG.NumSpaceTimeStreams  = 1; %[1 1];
        configParamsWG.NumTransmitAntennas  = 1; %2;
    end
    configParamsWG.APEPLength       = 2048/1; % 8x and APEPLength <= PSDU Length
    configParamsWG.save_tvs         = 0; % if 1, output files as stored
    configParamsWG.show_figs        = 0;
    configParamsWG.windowingRatio   = 0;
    configParamsWG.bbFilter_en      = 1;
    configParamsWG.bbFilterSbAtten  = 112;
    configParamsWG.resampleFilterN  = 20*4; % increase this parameter to sharpen resample filter transition band (improve ACLR of ideal waveform)
    configParamsWG.resampleFilterB  = 13.5;
    configParamsWG.outputLevelDbfs  = -15;
    configParamsWG.outputSampleRate = [];
    configParamsWG.floatPoint_en    = 0;

    if 1 % set wv guard interval
        configParamsWG.GuardInterval = 0.8;
        disp('2024-04-17, If GuardInterval < 3.2us, the bbFilter significatly affects the band edge response and EVM')
        % configParamsWG.bbFilter_en = 0;5
    else
        configParamsWG.GuardInterval = 3.2;
    end

    if strcmpi(signalFormat, 'HT20')
        configParamsWG.bbFilter_en  = 1;
    end

    if 0 % set wv over sampling ratio
        configParamsWG.osr          = [];
    else
        configParamsWG.osr          = 4;
        % configParamsWG.osr = 6;
        % configParamsWG.osr = 8;
        % configParamsWG.osr = 10;
    end

    if strcmpi(simulation_trx, 'rx') && isAdjacentChannelRejection_sim
        configParamsWG.osr          = 4;
    end

    if 1
        configParamsWG.resampleFilter_method    = 'waveformGen';
        configParamsWG.resampleFilter_en        = 0;
    else
        configParamsWG.resampleFilter_method    = 'resampleFun'; % bad
        configParamsWG.resampleFilter_en        = 1;
    end

    if 0 % set wv number of packets and
        configParamsWG.NumPackets = 10;
    else
        configParamsWG.NumPackets = 4;
    end


    if 1*0 % set wv idle time
        disp('If IdleTime = 40us, papr is higher but less noise ?')
        configParamsWG.IdleTime = 40e-6;
    elseif 2*0
        disp('If IdleTime = 5us, papr is lower but more noise ?')
        configParamsWG.IdleTime = 5e-6;
    else
        configParamsWG.IdleTime = 10e-6; % user-defined
    end

    if 0 % save waveform
        configParamsWG.save_tvs     = 1;
        configParamsWG.fileFormat   = "TXT_ADS";
        configParamsWG.wlanBand     = "5GHz";
        configParamsWG.note         = "";
    end

% Generate WIFI Waveform 

    % get wv configurations
    configWG = ovwrConfig(waveformGenWLAN,configParamsWG);

    % get wv parameters
    if 1
        passParams = waveformGenWLAN.runWG(configWG);
    else
        passParams = wvParams; % load wv
        passParams.cfgFormat = wvParams.configWG;
    end
    % export wv
    wv = passParams.finalSignal;
    fs = passParams.outputSampleRate;
    fsMHz = fs / 1e6;
    cfg = passParams.cfgFormat;
    dlSlots = passParams.usefulWvfm;
    cbw = passParams.currentBw;

    % set transmitter power
    wv = powerDbm(wv,'set',paOutputPower_target_dBm,'dBm',dlSlots);
    transmitPwr_dBm = powerDbm(wv,'rms',[],'dBm',dlSlots);

    if 1 % store original wv
        wv_org = wv;
        wv = wv_org;
    end

    if 1 % initialize PLOT_WLAN
        rbwKHz = 10 % set plt rbwKHz
        pltWlan = PLOT_WLAN(cfg);
        pltWlan.RbwKHz = rbwKHz;
        pltComm = PLOT_COMM([], [], fsMHz, rbwKHz);
    else
        pltWlan = wvParams.PlotACLR;
        pltComm = wvParams.PlotComm;
    end

    if 1 % plot wv aclr / filter sa
        aclrVal_wv = pltWlan.aclr(pltWlan, wv(dlSlots), [], fsMHz, 10, {091201,'wv(dlSlots)'})
        pltComm = PLOT_COMM('aclr', aclrVal_wv, [], RbwKHz, Rohm); % set plt object !!
        pltComm.sa(pltComm, wv(dlSlots), {091801,'wv(dlSlots)'});
        pltComm.sa(pltComm, passParams.bbFilterParams.bbFilterCoeffs, {091801, 'bbFilter'});
    end

    if 0 % 2024-04-16, get wv packet timing line
        timingInfo = pltWlan.timing(pltWlan, wv, fsMHz, [], {091202,dispLegd_format})
        durationData = timingInfo.tDATA;
        durationData_start_us = (timingInfo.tPREAMBLE_NonHT+timingInfo.tPREAMBLE)
        durationData_end_us = timingInfo.txTIME
        dlSlots_tmp = [durationData_start_us*1e-6*fs:durationData_end_us*1e-6*fs]';1
        pltComm.TIME(wv, {123}, [], 'usec')
    end

    if 1 % plot wv ccdf
        ccdfVal_wv = PLOT_COMM.ccdf([], wv, 0.01, {[091203], [dispLegd_format], 'CCDF'})
        ccdfVal_wv_dlSlots = PLOT_COMM.ccdf([], wv(dlSlots), 0.01, {[091203], [dispLegd_format+".dlSlots"], 'CCDF'})
        ccdfVal_wv_userDefine = PLOT_COMM.ccdf([], wv(1:16e-6*fs), 0.01, {[091203], dispLegd_format+".userDefine", 'CCDF'})
    end

    if 0 % evaluate wv
        close all
        demod_wv = wlanDemodulation(wv,passParams,'tx')
        demod_wv.evmRMS
    end

    if 1 % export - wvParams
        passParams.finalSignal = wv;
        clear wvParams
        % -
        wvParams = passParams;
        wvParams.configWG = configWG;
        % -
        wvParams.SignalFormat = signalFormat;
        wvParams.MCS = MCS_index;
        wvParams.NumPackets = configWG.NumPackets;
        wvParams.SampleRateMHz = wvParams.outputSampleRate/1e6;
        wvParams.PlotACLR = pltWlan;
        wvParams.PlotComm = pltComm;
    end

    if ~isPaModel
        fileName = sprintf('WIFI_%s_MCS%d_%.2fdBm_NumPkts%d_%dMHz.mat',...
            signalFormat, MCS_index, transmitPwr_dBm,...
            configWG.NumPackets, wvParams.SampleRateMHz)
        fullFilePath = fullfile(folderPath, fileName);
        save(fullFilePath, 'wvParams');
        continue
    end
% Transmitter -  Apply Noise - Before PA

    isNoiseFloorBeforePA = 0; % set noise floor to signal before PA
    if 1 * isNoiseFloorBeforePA
        n0_dBmHz_set = -150 % set noise psd in dBm/Hz
        [wvNoise, noise0] = rfSim_noise(wv, 'NoisePowerDbmHz', n0_dBmHz_set, fs); % create noise wv
        noise0_dBmHz = powerDbm(noise0,'rms',[],fs,dlSlots)

        % export
        wv = wvNoise;
        if 0 % evaluate wv
            passParams.finalSignal = wvNoise;
            demod_txNoise = wlanDemodulation(wvNoise, passParams, 'tx')
            demod_txNoise.evmRMS
        end
    end
% Transmitter - Characterize PA and Signal

    if 1
        paLinearGain_dB = 34.8; % set pa linear gain in dB
    else
        paLinearGain_dB = lut_mx3(1,2) - lut_mx1(1,1);
    end

    if ~isPaMultiStage
        clear paModel
        paModel.Method = "Cubic polynomial";
        paModel.AMPMConversion = 0.0;
        paModel.LinearGain = paLinearGain_dB;

        switch paModel_type
            case 'LUT'
                paModel.Method = 'Lookup table'; % set pa model method to 'Lookup table'
                dispLegd_LUT = paModel_LUT_sheet

                % set pa model LUT from paData.xlsx
                [paModel.LookupTable, LookupTable] = paLookupTable_process('paData.xlsx',dispLegd_LUT,[-30,2],091204,1);
                paLinearGain_dB = paModel.LookupTable(1,2) - paModel.LookupTable(1,1);

            case 'LUT_userDefine' % set LUT by user define
                clear setLUT
                setLUT.PinLinearDbm = [-40:5:-15];
                setLUT.GainLinearDb = 34.5;
                if 1 % set lutCase
                    lutCase = {'exponentExpension', 'exponent', 'lutCase1'};
                    lutCase = {'exponentExpension', 'linear', 'lutCase2'};
                    % lutCase = {'linearExpension', 'exponent', 'lutCase3'};
                    % lutCase = {'linearExpension', 'linear', 'lutCase4'};
                    % lutCase = {'compression', 'exponent', 'lutCase5'};
                    % lutCase = {'compression', 'linear', 'lutCase6'};
                end
                lutGainType = lutCase{1};
                lutPhasType = lutCase{2};
                dispLegd_LUT = ['Gain: ',lutGainType, ', Phase: ',lutPhasType]

                if 1 % reference of pout dbm
                    setLUT.PoutDbm = [setLUT.PinLinearDbm + setLUT.GainLinearDb, 22.5  27.6 30.3 31.0 31.7 32.5 33.1]';
                end
                switch lutGainType
                    case 'exponentExpension'
                        setLUT.GainCompressionDb = [zeros(size(setLUT.PinLinearDbm)), 0.1 0.2 -0.15 -0.4 -0.7 -0.9 -1.2]';
                    case 'linearExpension'
                        setLUT.GainCompressionDb = [0 0.01 0.02 0.04 0.06 0.1, 0.05 0. -0.15 -0.4 -0.7 -0.9 -1.2]';
                    case 'compression'
                        setLUT.GainCompressionDb = [zeros(size(setLUT.PinLinearDbm)), -0.05 -0.18 -0.3 -0.4 -0.5 -0.8 -1]';
                end
                switch lutPhasType
                    case 'exponent'
                        setLUT.PhaseShiftDeg = [zeros(size(setLUT.PinLinearDbm)), 0.04 0.22 0.5 0.6 1 0.4 -0.5]';
                    case 'linear'
                        setLUT.PhaseShiftDeg = [0 0 0.01 0.02 0.04 0.08, 0.16 0.32 0.2 0.15 0.1 0 -0.5]';
                        setLUT.PhaseShiftDeg = [0 0 0.01 0.02 0.04 0.08, 0.12 0.24 0.36 0.15 0.1 0 -0.5]';
                end
                if  numel(setLUT.PoutDbm)~=numel(setLUT.GainCompressionDb)|| numel(setLUT.GainCompressionDb)~=numel(setLUT.PhaseShiftDeg)
                    error('check the setLUT, length should be equal')
                end
        end

    else % multi-stages pa
        clear paModel paMd1 paMd2 paMd3 paMd4
        switch paModel_type
            case 'OP1dB'
                paMd1.Method = "Cubic polynomial"; % pa model 1
                paMd1.LinearGain = 12;
                paMd1.TOISpecification = "OP1dB";
                paMd1.OP1dB = 14;
                paMd2.Method = "Cubic polynomial";
                paMd2.LinearGain = 11;
                paMd2.TOISpecification = "OP1dB";
                paMd2.OP1dB = 30;
                paMd3.Method = "Cubic polynomial";
                paMd3.LinearGain = 10;
                paMd3.TOISpecification = "OP1dB";
                paMd3.OP1dB = 38;
                paMd3.AMPMConversion = 0;
                paMd3.PowerLowerLimit = 10.8000;
                paMd3.PowerUpperLimit = 20.8000;
                paModel = {paMd1, paMd2, paMd3};
            case 'OIP3'
                paMd1.Method = "Cubic polynomial"; % pa model 1
                paMd1.LinearGain = 12;
                paMd1.TOISpecification = "OIP3";
                paMd1.OIP3 = 14+10;
                paMd2.Method = "Cubic polynomial";
                paMd2.LinearGain = 11;
                paMd2.TOISpecification = "OIP3";
                paMd2.OIP3 = 30+11;
                paMd3.Method = "Cubic polynomial";
                paMd3.LinearGain = 10;
                paMd3.TOISpecification = "OIP3";
                paMd3.OIP3 = 38+11;
                paMd3.AMPMConversion = 0;
                paMd3.PowerLowerLimit = 10.8000;
                paMd3.PowerUpperLimit = 20.8000;
                paModel = {paMd1, paMd2, paMd3};
            case 'LUT'
                [lut_mx1, LookupTable1] = paLookupTable_process('paData.xlsx','0614SIMCW1M',[-50,50],[],1); % get pa1 lut from excel
                [lut_mx2, LookupTable2] = paLookupTable_process('paData.xlsx','0614SIMCW2M',[-50,50],[],1); % get pa2 lut from excel
                [lut_mx3, LookupTable3] = paLookupTable_process('paData.xlsx','0614SIMCW3M',[-50,50],[],1); % get pa1 lut from excel
                paMd1.Method = "Lookup table"; % pa model 1
                paMd2.Method = "Lookup table";
                paMd3.Method = "Lookup table";
                if ~strcmpi(paModel_methology,'RFSIMPA')
                    paMd1.LookupTable = lut_mx1;
                    paMd2.LookupTable = lut_mx2;
                    paMd3.LookupTable = lut_mx3;
                else
                    paMd1.LookupTable = LookupTable1;
                    paMd2.LookupTable = LookupTable2;
                    paMd3.LookupTable = LookupTable3;
                end
                paModel = {paMd1, paMd2, paMd3};
                if 0
                    figure(0530), subplot(1,3,1), plot(LookupTable.PinDbm, LookupTable.PoutDbm,'DisplayName',dispLegd_LUT, 'LineWidth',2), title('AMAM'), hold on, grid minor, legend
                    figure(0530), subplot(1,3,2), plot(LookupTable.PoutDbm, LookupTable.GainDb-0*LookupTable.GainDb(1),'DisplayName',dispLegd_LUT, 'LineWidth',2), title('GainCompression'), hold on, grid minor, legend
                    figure(0530), subplot(1,3,3), plot(LookupTable.PoutDbm, LookupTable.PhaseShiftDeg,'DisplayName',dispLegd_LUT, 'LineWidth',2), title('AMPM'), hold on, grid minor, legend
                end
                if 0
                    paMd4.Method = "Lookup table";
                    paMd4.LookupTable = rfSim_pa_amam_model(LookupTable3.PoutDbm, LookupTable3.GainDb, LookupTable3.GainDb(1), 0.5, 0619)
                    paModel{end} = paMd4
                end
            case {'OP1dB+LUT','mixing case'}
                paMd1.Method = "Cubic polynomial"; % pa model 1
                paMd1.LinearGain = 12;
                paMd1.TOISpecification = "OP1dB";
                paMd1.OP1dB = 14;
                paMd2.Method = "Cubic polynomial";
                paMd2.LinearGain = 11;
                paMd2.TOISpecification = "OP1dB";
                paMd2.OP1dB = 30;
                paMd4.Method = "Lookup table";
                [lut_mx4, LookupTable4] = paLookupTable_process('paData.xlsx','0614SIMCW3M',[-50,50],[],1);
                if ~strcmpi(paModel_methology,'RFSIMPA')
                    paMd4.LookupTable = lut_mx4;
                else
                    paMd4.LookupTable = LookupTable4;
                end
                paModel = {paMd1, paMd2, paMd4};
        end
    end

% Transmitter - CFR - Release V2.0 - 2024-09-25

    isCfrBeforePA = 0 % set cfr before PA
    ccdfVal_wv = PLOT_COMM.ccdf([], wv, 0.01) % check ccdf

    if isCfrBeforePA
        cfr_paprDb = 8 % set cfr limitor
        cfr_nIteration = 20 % set iteration times

        wv_cfr = dfe_cfr_nco_2(wv, cfr_paprDb, [], [], cfr_nIteration, []);
        powerDbm(wv_cfr) % check pwr
        ccdfVal_wvCfr = PLOT_COMM.ccdf([], wv_cfr, 0.01) % check ccdf
        demod_wvCfr = wlanDemodulation(wv_cfr,passParams,'tx') % check demod.
        demod_wv = wlanDemodulation(wv,{passParams,'wv'},'tx')
        wv = wv_cfr;
    end
% Transmitter -  Generate PA Signal and Tuning to Target Output Power

    paOutputPower_tolerance_dB = 0.2; % set pa output power tolerance
    deltaTarget_dB = paOutputPower_tolerance_dB;
    ii = 1;
    paInputPowerDbm_set = paOutputPower_target_dBm - paLinearGain_dB; % get initial pa input power
    SIMULATING_tuningPower = 1;
    paOutputPower_dBm_vec = -1e10; % initial pa output power vector

    while SIMULATING_tuningPower
        clear paInSignal paSignal
        close all

        if 1 % set pa input power
            paInSignal = wv./rms(wv(dlSlots)).*10^(paInputPowerDbm_set/20)*sqrt(1/1000);
            paInSignal_dBm = 10*log10(mean(abs(paInSignal(dlSlots)).^2)) + 30
        else
            paInSignal = powerDbm(wv,'set',paInputPowerDbm_set,'dBm',dlSlots);
            paInSignal_dBm = powerDbm(paInSignal,'rms',[],[],dlSlots)
        end

        if 0 % check demod.
            demod_paIn = wlanDemodulation(paInSignal,passParams,'tx')
        end

        if isfield(paModel,'AMPMConversion') % set pa input boundary for AMPMConversion
            paModel.PowerLowerLimit = paInSignal_dBm - 5; % 5 is margin
            paModel.PowerUpperLimit = paInSignal_dBm + 11; % 11 is PAPR
        end

        % apply pa model
        switch paModel_methology
            case 'commMemorylessNonlinearity'
                [paSignal, paSignal_info_stages, paSignal_stages] = rfSim_pa_stages(paInSignal,'COMM',paModel,dlSlots);
            case 'RFSIMPA'
                [paSignal, paSignal_info_stages, paSignal_stages] = rfSim_pa_stages(paInSignal,'RFSIM',paModel,dlSlots); % generate pa output wv
        end

        if 1 % evaluate - pa output
            demod_paSignal = wlanDemodulation(paSignal, passParams, 'tx')
            paSignal_dBm = powerDbm(paSignal, 'rms', [], 'dBm', dlSlots)
            PLOT_COMM.time([], paInSignal, fs/1e6, 'us', {061203, ['pa input']});
            if isPaMultiStage
                PLOT_COMM.time([], paSignal_stages(:,1), fs/1e6, 'us', {061203, ['pa stage 1']});
                PLOT_COMM.time([], paSignal_stages(:,2), fs/1e6, 'us', {061203, ['pa stage 2']});
                PLOT_COMM.time([], paSignal_stages(:,3), fs/1e6, 'us', {061203, ['pa stage 3']});
                pltComm.AMAM([paInSignal(dlSlots), paSignal_stages(dlSlots,1)], {061901, ['pa stage 1']}, [], [], lut_mx1(1,2)-lut_mx1(1,1));
                pltComm.AMAM([paInSignal(dlSlots), paSignal_stages(dlSlots,2)], {061901, ['pa stage 2']}, [], [], lut_mx2(1,2)-lut_mx1(1,1));
                pltComm.AMAM([paInSignal(dlSlots), paSignal_stages(dlSlots,3)], {061901, ['pa stage 3']}, [], [], lut_mx3(1,2)-lut_mx1(1,1));
            end
        end

% Transmitter - Flatness - After PA

        paSignal_org = paSignal; % store org

        if isPaFlatness
            paFlatnessMax_dB    = 2 % set max flatness
            paFlatFreqMHz       = cbw/2/1e6 * [0:1:configParamsWG.osr]; % implement flat. freq. vector in MHz
            paFlatMagnDb        = -paFlatnessMax_dB * [0:1:configParamsWG.osr] + paFlatnessMax_dB/2; % implement flat. mag. vector in Db
            nTapsRange_paFlat   = 10e3; % set range of taps for flat. fir
            [paSignalFlat, b_paFlat] = rfSim_fltaness(paSignal, fs, paFlatMagnDb, paFlatFreqMHz, nTapsRange_paFlat); % set wv with flatness
            if 0 % debug
                paSignalFlat_tmp = conv(paSignal, b_paFlat, 'same');
                ccdfVal_paSignal = PLOT_COMM.ccdf([], paSignal, 0.01, {[1021], ['paSignal'], 'CCDF'})
                ccdfVal_paSignalFlat = PLOT_COMM.ccdf([], paSignalFlat, 0.01, {[1021], ['paSignalFlat'], 'CCDF'})
            end

            if isPaFlatnessCompensateLoss % compensate flatness loss
                paInSignal_dBm = powerDbm(paInSignal,'rms',[],[],dlSlots)
                paSignal_dBm = powerDbm(paSignal,'rms',[],[],dlSlots)
                paSignalFlat_dBm = powerDbm(paSignalFlat,'rms',[],[],dlSlots)
                paSignalFlat = powerDbm(paSignalFlat, 'set', paSignal_dBm, 'dBm', dlSlots);
            end

            % export
            paSignal = paSignalFlat;

            if 0 % evaluate - pa flatness
                aclrVal_paOut = pltComm.ACLR(paSignal_org(dlSlots), {091901, 'paSignal'},[],[],[],[],10,[],'dBm'); % plot aclr
                aclrVal_paIn = pltComm.ACLR(paInSignal(dlSlots), {091901, 'paInSignal'},[],[],[],[],10,[],'dBm'); % plot aclr
                aclrVal_paOutFlat = pltComm.ACLR(paSignalFlat(dlSlots), {091901, 'paSignalFlat'},[],[],[],[],10,[],'dBm'); % plot aclr
                powerDbm(bb)
                aclrVal_b_paFlat = pltComm.ACLR(b_paFlat, {091901, 'b_paFlat'},[],[],[],[],10,[],'dBr'); % plot aclr
            end

        else
            paFlatnessMax_dB = 0;
            b_paFlat = 1;
        end

        % check pa output power
        paSignal_dBm = powerDbm(paSignal,'rms',[],[],dlSlots)
        paOutputPower_dBm_vec = [paOutputPower_dBm_vec, paSignal_dBm];
        ii=ii+1;

        % 2024-03-21, To get the fixed output power
        deltaTarget_dB = round(paSignal_dBm - paOutputPower_target_dBm, 3)
        if deltaTarget_dB > paOutputPower_tolerance_dB
            paInputPowerDbm_set = paInputPowerDbm_set - paOutputPower_tolerance_dB/2; % decrease the boostDb
            paSignal_tmp = paSignal;
            continue
        elseif any(diff(paOutputPower_dBm_vec)<0) || ii>20 % saturation
            SIMULATING_tuningPower = 0;
            paSignal = paSignal_tmp;
            paSignal_dBm = paOutputPower_dBm_vec(end-1);
        elseif deltaTarget_dB < 0
            paInputPowerDbm_set = paInputPowerDbm_set - deltaTarget_dB; % increase the boostDb
            paSignal_tmp = paSignal;
            continue
        else
            SIMULATING_tuningPower = 0;
        end
    end % SIMULATING_tuningPower

    % export - tx final signal
    passParams.finalSignal = paSignal;
    
    if isPaMultiStage % evaluate - pa output stages
        demod_casd_pa_stage_1 = wlanDemodulation(paSignal_stages(:,1), passParams, 'tx');
        demod_casd_pa_stage_2 = wlanDemodulation(paSignal_stages(:,2), passParams, 'tx');
        demod_casd_pa_stage_3 = wlanDemodulation(paSignal_stages(:,3), passParams, 'tx');
    end
% Transmitter - wlanDemodulation and Measurements

    if deltaTarget_dB < 0 || deltaTarget_dB > paOutputPower_tolerance_dB + 0.1 % check deltaTarget_dB
        error('pa output power could not met the target power')
    end
    demod_finalSignal = wlanDemodulation(passParams.finalSignal, passParams, 'tx') % demodulate wv
% Transmitter - Summary

    mappingPass = containers.Map({1,0},{"Pass","Fail"}); % set mapping to pass or fail
    clear summary

    % signal format and mcs index
    summary.Format = signalFormat_vec(k);
    summary.MCS = MCS_index;

    if 1 % pa impairments
        summary.PaFlatnessMaxDbPerCBW = paFlatnessMax_dB;
    end

    % summary power
    summary.TxPowerDbm_Spec = paOutputPower_target_dBm;
    summary.TxPowerDbm_Final = round(paSignal_dBm, 2);
    summary.TxPower_Pass = string(arrayfun(@(x)mappingPass(x),...
        (summary.TxPowerDbm_Final>summary.TxPowerDbm_Spec-paOutputPower_tolerance_dB), 'UniformOutput', false));

    if isPaMultiStage
        paInPowerDbm = round(paInSignal_dBm)
        pa1PowerDbm = round(powerDbm(paSignal_stages(dlSlots,1),'rms'),2)
        pa2PowerDbm = round(powerDbm(paSignal_stages(dlSlots,2),'rms'),2)
        pa3PowerDbm = round(powerDbm(paSignal_stages(dlSlots,3),'rms'),2)
        summary.TxPowerDbm_PaStage =  {[paInPowerDbm, pa1PowerDbm, pa2PowerDbm, pa3PowerDbm]};
    end

    % summary demodulation evm
    summary.EVM_RMS_Spec =  specEVM;
    if isPaMultiStage
        evm_rms_casd_pa_stage = round([max(demod_casd_pa_stage_1.evmRMS),max(demod_casd_pa_stage_2.evmRMS),max(demod_casd_pa_stage_3.evmRMS)], 2);
        summary.EVM_RMS_Demod_Casd_PaStage = {evm_rms_casd_pa_stage};
    else
        evm_rms_casd_pa_stage = round(max(demod_finalSignal.evmRMS), 2);
        summary.EVM_RMS_Demod_Casd_PaStage = evm_rms_casd_pa_stage;
    end
    summary.EVM_RMS_Pass = string(arrayfun(@(x)mappingPass(x),...
        (evm_rms_casd_pa_stage(end)<summary.EVM_RMS_Spec), 'UniformOutput', false));

    % summary raw evm
    summary.EVM_RMS_Raw_Casd_PaStage = [];
    for nStage=1:numel(paModel)
        evm_rms_raw_casd_paStage = round(evm(paInSignal(dlSlots,1),paSignal_stages(dlSlots,nStage),'dB',1,[fs,-cbw/2+2e6,cbw/2-2e6]), 2)
        summary.EVM_RMS_Raw_Casd_PaStage = [summary.EVM_RMS_Raw_Casd_PaStage, evm_rms_raw_casd_paStage];
    end
    summary.EVM_RMS_Raw_Casd_PaStage = {summary.EVM_RMS_Raw_Casd_PaStage}

    % summary flatness
    summary.SpectralFlatness_Pass = string(arrayfun(@(x)mappingPass(x),...
        all(demod_finalSignal.isFlatPass), 'UniformOutput', false));

    % summary mask
    summary.SpectralMask_Pass = string(arrayfun(@(x)mappingPass(x),...
        all(demod_finalSignal.isMaskPass), 'UniformOutput', false));

    summary_tab_tmp = struct2table(summary); % convert struct to table
    if exist("summary_tab",'var') % vertcat
        summary_tab = vertcat(summary_tab,summary_tab_tmp)
    else
        summary_tab = summary_tab_tmp
    end

    % save
    if 1
        % clear wvParams fileName
        % % - 
        % wvParams = passParams;
        % wvParams.configWG = configWG;
        % % -
        wvParams.SignalFormat = signalFormat;
        wvParams.MCS = MCS_index; 
        wvParams.NumPackets = configWG.NumPackets;
        wvParams.SampleRateMHz = wvParams.outputSampleRate/1e6;
        wvParams.PaInSignal = paInSignal; % export
        wvParams.PaOutSignal_ideal = paSignal_org; % export
        wvParams.PaOutSignal = paSignal; % export
        wvParams.PlotACLR = pltWlan;
        wvParams.PlotComm = pltComm;
        if 1 % pa model
            wvParams.PaModel = paModel;
            wvParams.PaModelFlatnessPerCBW = paFlatnessMax_dB;
            wvParams.PaModelFlatFirCoefs = b_paFlat;
            wvParams.PaModelFlatCompensateLoss = logical(isPaFlatnessCompensateLoss);
        end
        % -
        if isPaFlatness & isPaFlatnessCompensateLoss
            fileName = sprintf('WIFI_%s_MCS%d_%.2fdBm_NumPkts%d_%dMHz_%s_Cp%ddB.mat',...
                signalFormat, MCS_index, paOutputPower_target_dBm,...
                configWG.NumPackets, wvParams.SampleRateMHz, paModel_LUT_sheet, paFlatnessMax_dB)
        else
            fileName = sprintf('WIFI_%s_MCS%d_%.2fdBm_NumPkts%d_%dMHz_%s_%ddB.mat',...
                signalFormat, MCS_index, paOutputPower_target_dBm,...
                configWG.NumPackets, wvParams.SampleRateMHz, paModel_LUT_sheet, paFlatnessMax_dB)
        end
        fullFilePath = fullfile(folderPath, fileName);
        save(fullFilePath, 'wvParams');
    end
end % k

% Package 

if 0
    varPackage = packageVars('paInSignal', 'paSignal','summary_tab')
end
% Fnunctions

function selectedVars = packageVars(varargin)
% packageVars packages the selected variables into a structure.
% Input: Variable names (strings)
% Output: Structure containing the selected variables

% Initialize an empty structure
selectedVars = struct();

% Loop through each input variable name
for i = 1:length(varargin)
    varName = varargin{i};  % Get the variable name

    if 1 % upper case
        newFieldName = [upper(varName(1)), varName(2:end)];  % Capitalize first letter
    else
        newFieldName = varName;
    end

    if evalin('caller', ['exist(''', varName, ''', ''var'')'])  % Check if variable exists in the caller workspace
        selectedVars.(newFieldName) = evalin('caller', varName);  % Package the variable into the structure
    else
        warning(['Variable ', varName, ' does not exist.']);
    end
end
end

function [SNROut_dB, NFcal_dB] = cal_snr(ref,mea_output,mea_input,bwInband,fs)
isFFTFreqSelect = 1;
if ~exist('bwInband','var')||isempty(bwInband)
    isFFTFreqSelect = 0;
elseif isscalar(bwInband)
    freq_1 = -bwInband/2;
    freq_end = bwInband/2;
else
    freq_1 = bwInband(1);
    freq_end = bwInband(2);
end

if isFFTFreqSelect
    freq_vec = fs/numel(ref) * [-numel(ref)/2+1:numel(ref)/2];
    [~,idx_cbw_1] = min(abs(freq_vec-freq_1));
    [~,idx_cbw_end] = min(abs(freq_vec-freq_end));
    idx_cbw = [idx_cbw_1:idx_cbw_end]';
    REF = fftshift(fft(ref));
    MEA_OUT = fftshift(fft(mea_output));
    MEA_IN = fftshift(fft(mea_input));
    ref = ifft(fftshift(REF(idx_cbw)));
    mea_output = ifft(fftshift(MEA_OUT(idx_cbw)));
    mea_input = ifft(fftshift(MEA_IN(idx_cbw)));
end

[~, SNRIn_dB] = evm(ref, mea_input);
[~, SNROut_dB] = evm(ref, mea_output);
NFcal_dB = SNRIn_dB - SNROut_dB;
end