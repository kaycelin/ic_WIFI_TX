% 2024-02-26, Draft, Referring to Matlab example "802.11be Transmitter Measurements"
% 2024-02-27, Add measurementSummary
% 2024-03-05, Use plot_wlan instead of ehtSpectralMaskTest
% 2024-03-06, Use plot_wlan instead of ehtPlotTxSpectralFlatness
% 2024-03-07, Use plot_wlan instead of ehtTxEVMConstellationPlots
% 2024-03-14, timingSync: overSamplingFreq, use over sampling freq to improve the accuracy
% 2024-03-15, Update 802.11ax (HE) demodulation
% 2024-03-25, Issue, users and signal puncturing ???
% 2024-03-25, add isTxMeasurement
% 2024-04-16, calculate RAW EVM
% 2024-04-16, Optimized the decimation filter by a1 and enable_cbw_tolearnce
% 2024-04-23, Add subcarier allocation related constants for decimation fir design
% 2024-04-24, Add dfe_decimate
% 2024-05-08, frequency selected by fft
% 2024-05-13, update numTxAnt for MIMO structure
% 2024-06-19, change function from wlanReceiver to wlanDemodulation, and
% add txWaveform input parameters
%% 2024-09-16, apply PLOT_WLAN instead of plot_wlan


function output = wlanDemodulation(txWaveform, passParams_txWaveform_dispLgndCell, isTxMea_sampleRate, cfgFormat, packetParams, inBits, isTxMea)

% initialize
if ~exist('txWaveform','var')||isempty(txWaveform)
    txWaveform = []; % 2024-06-19
end
dispLgnd = ""; % init. empty
if isstruct(passParams_txWaveform_dispLgndCell)
    parms = passParams_txWaveform_dispLgndCell;
    if isempty(txWaveform)
        txWaveform = parms.finalSignal;
    end
    fs_osf = parms.outputSampleRate;
    cfg = parms.cfgFormat;
    packetParams = parms.packetParams;
    txPSDU = {parms.inBits};

    if exist('isTxMea_sampleRate','var')&&~isempty(isTxMea_sampleRate)
        isTxMeasurement = isTxMea_sampleRate; % 2024-03-25, add isTxMeasurement
    end
elseif iscell(passParams_txWaveform_dispLgndCell)
    parms = passParams_txWaveform_dispLgndCell{1};
    if isempty(txWaveform)
        txWaveform = parms.finalSignal;
    end
    fs_osf = parms.outputSampleRate;
    cfg = parms.cfgFormat;
    packetParams = parms.packetParams;
    txPSDU = {parms.inBits};
    if exist('isTxMea_sampleRate','var')&&~isempty(isTxMea_sampleRate)
        isTxMeasurement = isTxMea_sampleRate; % 2024-03-25, add isTxMeasurement
    end

   dispLgnd = "."+ passParams_txWaveform_dispLgndCell{end};
else
    if isempty(txWaveform)
        txWaveform = passParams_txWaveform_dispLgndCell;
    end
    fs_osf = isTxMea_sampleRate;
    cfg = cfgFormat;
    packetParams;
    if iscell(inBits)
        txPSDU = inBits;
    else
        txPSDU = {inBits};
    end
    isTxMeasurement = isTxMea;
end
switch upper(isTxMeasurement)
    case {'TX'}
        isTxMeasurement = 1;
    case {'RX'}
        isTxMeasurement = 0;
    otherwise
        error('isTxMeasurement ?')
end

% get numPackets/numUsers/psduLen/ppduType
if isstruct(packetParams)&&isfield(packetParams,'NumPackets')
    numPackets = packetParams.NumPackets;
else
    numPackets = packetParams;
end
% Get indices for accessing each field within the time-domain packet.
ind = wlanFieldIndices(cfg);
cbw = cfg.ChannelBandwidth;
switch class(cfg)
    case {'wlanEHTMUConfig'}
        signalFormat = "EHT"+erase(cbw, ["CBW"]);
        psduLen = psduLength(cfg).*8;
        ofdmInfo = wlanEHTOFDMInfo("EHT-Data",cfg,1);
        pktLength = double(ind.EHTData(2));
        numUsers = cfg.NumUsers;
    case {'wlanHESUConfig','wlanHEMUConfig'}
        signalFormat = "HE"+erase(cbw, ["CBW"]);
        psduLen = getPSDULength(cfg).*8;
        ofdmInfo = wlanHEOFDMInfo("HE-Data",cfg,1);
        pktLength = double(ind.HEData(2));
        try
            numUsers = cfg.NumUsers;
        catch
            numUsers = 1;
        end
    case {'wlanVHTConfig'}
        signalFormat = "VHT"+erase(cbw, ["CBW"]);
        psduLen = cfg.PSDULength.*8;
        ofdmInfo = wlanVHTOFDMInfo("VHT-Data",cfg);
        pktLength = double(ind.VHTData(2));
        numUsers = 1;
    case {'wlanHTConfig'}
        signalFormat = "HT"+erase(cbw, ["CBW"]);
        psduLen = cfg.PSDULength.*8;
        ofdmInfo = wlanHTOFDMInfo("HT-Data",cfg);
        pktLength = double(ind.HTData(2));
        numUsers = 1;
end

if numUsers > 1
    ppduType = "OFDMA";
    error('2024-05-13, not support ofdma!')
else
    ppduType = "non-OFDMA";
end

numTxAnt = size(txWaveform,2); % 2024-05-13, update numTxAnt for MIMO structure

fs = ofdmInfo.SampleRate;
osf = fs_osf / fs;
rxWaveform = [];
for k=1:numTxAnt
    if osf==1
        % remove filter
        rxWaveform(:,k) = txWaveform(1:osf:end,k);
        delay_firdec = 0;
    elseif 1*1
        % 2024-04-23, Add subcarier allocation related constants for
        % decimation fir design
        NSR = max(ofdmInfo.ActiveFrequencyIndices);
        scs = fs / ofdmInfo.FFTLength;
        switch parms.cfgFormat.GuardInterval
            case {0.8,'short','long'}
                fstart = (fix(NSR+fix(fs/2/scs))/2)*scs;
            otherwise
                fstart = (NSR+0)*scs;
        end
        ftrans = (fs/2-fstart)/1;
        % 2024-04-24, Add dfe_decimate
        [rxWaveform(:,k), b_cell] = dfe_decimate(txWaveform(:,k),fs_osf,fs,fstart,ftrans);
        delay_firdec = 0;
        if 0
            pltParms.rbw = scs;
            plot_comm(b_cell{:}, fs_osf, 'spectra', pltParms, {041801, 'decimate filter'}); hold on
            plot_comm(txWaveform(:,k), fs_osf, 'spectra', pltParms, {041801, 'tx wv'}); hold on
            plot_comm(rxWaveform(:,k), fs_osf, 'spectra', pltParms, {041801, 'tx wv + filter'}); hold on
            pltParms.unit = 'us';
            plot_comm(txWaveform(:,k), fs_osf, 'time', pltParms, {041801, 'tx wv time'}); hold on
        end
    elseif 2*0
        % 2024-05-08, frequency selected by fft, however the method
        % contribute the fft noise
        TXWV = fftshift(fft(txWaveform(:,k)));
        FREQ_OSF = fs_osf/numel(TXWV)*[-numel(TXWV)/2+1:numel(TXWV)/2]';
        % figure, plot(FREQ_OSF,abs(TXWV))
        [~,idx_cbw_1] = min(abs(FREQ_OSF+parms.currentBw/2));
        [~,idx_cbw_end] = min(abs(FREQ_OSF-parms.currentBw/2));
        TXWV_CBW = TXWV(idx_cbw_1:idx_cbw_end);
        % figure, plot(FREQ_OSF(idx_cbw_1:idx_cbw_end),abs(TXWV_CBW))
        TXWV_CBW_OSF = zeros(numel(TXWV),1);
        TXWV_CBW_OSF(idx_cbw_1:idx_cbw_end,:) = TXWV_CBW;
        % figure, plot(FREQ_OSF,abs(TXWV_CBW_OSF))
        txwv_cbw_osf = ifft(fftshift(TXWV_CBW_OSF));
        % figure, plot(abs(TXWV_CBW_OSF))
        rxWaveform(:,k) = txwv_cbw_osf(1:osf:end);
    else
        % Matlab: designMultirateFIR
        if 1 % 2024-04-16, Optimized the decimation filter by a1 and enable_cbw_tolearnce
            a1 = 2;
            enable_cbw_tolearnce = 1;
        else
            a1 = 0;
            enable_cbw_tolearnce = 0;
        end
        % Design resampling filter.
        aStop = 20 * a1; % Stopband attenuation
        SCS = fs/ofdmInfo.FFTLength; % Subcarrier spacing
        txbw = max(abs(ofdmInfo.ActiveFrequencyIndices))*2*SCS; % Occupied bandwidth
        txbw = txbw + enable_cbw_tolearnce*(fs - txbw)/2;
        [L,M] = rat(1/osf);
        maxLM = max([L M]);
        R = (fs-txbw)/fs;
        TW = 2*R/maxLM; % Transition width
        firdec = designMultirateFIR(L,M,TW,aStop,'SystemObject',true);
        % Resample the waveform to baseband
        rxWaveform(:,k) = firdec(txWaveform(:,k));
        delay_firdec = grpdelay(firdec,1);          % Group delay of downsampling filter
        if 0
            pltParms.rbw = 10e3;
            pltParms.winType = 'hann';
            plot_comm(txWaveform(:,k), fs_osf, 'spectra', pltParms, {041801, 'txWaveform'}); hold on
            plot_comm(firdec.Numerator(:), fs_osf, 'spectra', pltParms, {041801, 'decimate filter'}); hold on
            wvInput = conv(txWaveform(:,k), firdec.Numerator, 'same');
            plot_comm(wvInput, fs_osf, 'spectra', pltParms, {041801, 'after decimate filter'}); hold on
            plot_comm(rxWaveform(:,k), fs, 'spectra', pltParms, {041801, 'rxWaveform'}); hold on
        end
    end
end
rxWaveform_ant = rxWaveform;
if ~isTxMeasurement * 0
    agc = comm.AGC; % Automatic gain control
    % Pass each channel through AGC
    for ic = 1:size(rxWaveform_ant,2) % cfg.NumSpaceTimeStreams or cfg.NumTransmitAntennas ???
        rxWaveform(:,ic) = agc(rxWaveform_ant(:,ic));
        reset(agc);
    end
end
%% Receiver Processing
% Define the minimum detectable length of data, in samples.
minPktLen = double(ind.LSTF(2)-ind.LSTF(1))+1;

% Detect and process packets within the received waveform by using a while
% loop, which performs these steps.
% 1. Detect a packet by indexing into rxWaveform with the sample offset, searchOffset
% 2. Detect and process the first packet within rxWaveform
% 3. Detect and process the next packet by incrementing the sample index offset
% 4. Repeat until no further packets are detected
rxWaveformLength = size(rxWaveform,1);
if 0
    rmsEVM = zeros(numUsers,numPackets);
    eqSym = cell(numUsers,numPackets);
    evmPerSC = cell(numUsers,numPackets);
    pktOffsetStore = zeros(numPackets);
    isEVMPass = false(numUsers,numPackets);
    isFlatPass = false(numUsers,numPackets);
    isMaskPass = false(1,numPackets);
    isDecodePass = false(numUsers,numPackets);
    flatDeviation = cell(numUsers,numPackets);
else
    % 2024-05-13, update numTxAnt for MIMO structure
    % row is antenna no./ column is package no./ layer is user no.
    rmsEVM = zeros(numTxAnt,numPackets,numUsers);
    eqSym = cell(numUsers,numPackets);
    evmPerSC = cell(numUsers,numPackets);
    pktOffsetStore = zeros(1,numPackets);
    isEVMPass = false(numTxAnt,numPackets,numUsers);
    isFlatPass = false(numTxAnt,numPackets,numUsers);
    isMaskPass = false(numTxAnt,numPackets,numUsers);
    isDecodePass = false(numTxAnt,numPackets,numUsers);
    flatDeviation = cell(numUsers,numPackets);
end

nPkt = 1; % Start from No.1
searchOffset = 0; % Start at first sample (no offset)
numPacketErrors = 0; % init.
packetErrorRatePercent = zeros(numUsers,1);
offsetDebug = [];
nUser = 0;
while (searchOffset+minPktLen)<=rxWaveformLength % ......................nPkt
    [rxPacket, pktOffset, coarsefreqOff, fineFreqOff] = timingSync(rxWaveform,cfg,fs_osf,searchOffset);

    for nUser = 1:numUsers %.............................................nUser
        try
            numSTSs = cfg.User{nUser}.NumSpaceTimeStreams;
        catch
            numSTSs = cfg.NumSpaceTimeStreams; % set NumSpaceTimeStreams=1
        end

        for nSTS = 1:numSTSs %...........................................nSTS
            if isempty(rxPacket)
                if nUser==0
                    error('Synchronization Failed in First Packet')
                end
                break;
            end
            if 0 % debug
                nPkt
                offsetDebug = [offsetDebug; pktOffset, coarsefreqOff, fineFreqOff, searchOffset];
            end

            if isTxMeasurement % Spectral Mask Measurement
                switch class(cfg)
                    case {'wlanEHTMUConfig'}
                        idx_Data_start = osf*(ind.EHTData(1)-1)+1; % Upsampled start of EHT Data
                        idx_Data_end = osf*(ind.EHTData(2));
                    case {'wlanHESUConfig','wlanHEMUConfig'}
                        idx_Data_start = osf*(ind.HEData(1)-1)+1; % Upsampled start of HE Data
                        idx_Data_end = osf*(ind.HEData(2));
                    case {'wlanVHTConfig'}
                        idx_Data_start = osf*(ind.VHTData(1)-1)+1; % Upsampled start of VHT Data
                        idx_Data_end = osf*(ind.VHTData(2));
                    case {'wlanHTConfig'}
                        idx_Data_start = osf*(ind.HTData(1)-1)+1; % Upsampled start of HT Data
                        idx_Data_end = osf*(ind.HTData(2));
                end

                % Start of packet in txWaveform
                pktOffset_osf = round(osf*pktOffset) - delay_firdec;
                % Indices of EHT Data in txWaveform
                idxSignalFormatData(:,nPkt) = pktOffset_osf+(idx_Data_start:idx_Data_end);
                if nPkt == numPackets || 1 % plot spectral mask based on all packets
                    gatedEHTTData = txWaveform(idxSignalFormatData(:,nPkt));
                    if 0
                        output_mask.isMaskPass = ehtSpectralMaskTest(gatedEHTTData,fs,osf);
                    else % 2024-03-05, Use plot_wlan instead of ehtSpectralMaskTest
                        if 0
                            output_mask = plot_wlan(gatedEHTTData, fs_osf, 'spectral mask', cfg, numPackets, {0304,[signalFormat+dispLgnd],'Spectral Mask'});
                        else
                            if 0
                                output_mask = PLOT_WLAN.spectralMask([], gatedEHTTData, signalFormat, fs_osf/1e6,...
                                    {[0304, numSTSs, numPackets, nPkt], [signalFormat+dispLgnd], 'Spectral Mask'});
                            else
                                output_mask = PLOT_WLAN.spectralMask([], gatedEHTTData, signalFormat, fs_osf/1e6,...
                                    {[0304, numSTSs, numPackets, nPkt], [signalFormat+dispLgnd], sprintf('Spectral Mask'+".User%d.Ant%d.Pkt%d", nUser, nSTS, nPkt)});
                            end
                        end
                        isMaskPass(:,nPkt,nUser) = output_mask.isMaskPass;
                        outputPowerPerPktsDbm(:,nPkt,nUser) = output_mask.outputPowerPerPktsDbm;
                    end
                end
            end

            % Extract demoulated LTF field and perform and the channel
            % estimation
            switch class(cfg)
                case {'wlanEHTMUConfig'}
                    ehtLTF = rxPacket(ind.EHTLTF(1):ind.EHTLTF(2),:); % extract LTF signal
                    ehtLTFDemod = wlanEHTDemodulate(ehtLTF,'EHT-LTF',cfg,nUser); % LTF demodulation
                    [chanEst,pilotEst] = wlanEHTLTFChannelEstimate(ehtLTFDemod,cfg,nUser); % channel estimation
                    ofdmInfo_i = wlanEHTOFDMInfo('EHT-Data',cfg,nUser); % OFDM parameters
                case {'wlanHESUConfig','wlanHEMUConfig'}
                    heLTF = rxPacket(ind.HELTF(1):ind.HELTF(2),:);
                    heLTFDemod = wlanHEDemodulate(heLTF,'HE-LTF',cfg,nUser); % LTF demodulation
                    [chanEst,pilotEst] = wlanHELTFChannelEstimate(heLTFDemod,cfg,nUser); % channel estimation
                    if 0
                        rxHELTF = rxPacket((ind.HELTF(1):ind.HELTF(2)),:);
                        % HE-LTF demodulation
                        try
                            heltfDemod = wlanHEDemodulate(rxHELTF,'HE-LTF',cfg(1).ChannelBandwidth,cfg(1).GuardInterval, ...
                                cfg(1).HELTFType,[allocInfo.RUSizes(r) allocInfo.RUIndices(r)]);
                        catch
                            heltfDemod = wlanHEDemodulate(rxHELTF,'HE-LTF',cfg(1).ChannelBandwidth,cfg(1).GuardInterval, ...
                                cfg(1).HELTFType);
                        end
                        % Channel estimate
                        [heltfRU(r).ChanEst,heltfRU(r).PilotEst] = wlanHELTFChannelEstimate(heltfDemod,cfgUser(uidx));
                    end
                    % end
                    ofdmInfo_i = wlanHEOFDMInfo('HE-Data',cfg,nUser); % OFDM parameters
                case {'wlanVHTConfig'}
                    % It is used for MIMO channel estimation and pilot subcarrier tracking.
                    % The VHT-LTF includes one VHT long training symbol for each spatial stream indicated by the selected modulation and coding scheme (MCS).
                    % Each symbol is 4 μs long. A maximum of eight symbols are permitted in the VHT-LTF.
                    vhtLTF = rxPacket((ind.VHTLTF(1):ind.VHTLTF(2)), :);
                    vhtLTFDemod = wlanVHTLTFDemodulate(vhtLTF, cfg);
                    [chanEst, pilotEst] = wlanVHTLTFChannelEstimate(vhtLTFDemod, cfg);
                    ofdmInfo_i = wlanVHTOFDMInfo('VHT-Data',cfg); % OFDM parameters
                case {'wlanHTConfig'}
                    htLTF = rxPacket((ind.HTLTF(1):ind.HTLTF(2)), :);
                    htLTFDemod = wlanHTLTFDemodulate(htLTF,cfg);
                    chanEst = wlanHTLTFChannelEstimate(htLTFDemod, cfg);
                    ofdmInfo_i = wlanHTOFDMInfo('HT-Data',cfg); % OFDM parameters
            end

            % Measure spectral flatness for non-OFDMA PPDY type
            if strcmpi(ppduType,"Non-OFDMA")
                if isTxMeasurement
                    % Plot deviation against limits
                    % 2024-03-06, Use plot_wlan instead of ehtPlotTxSpectralFlatness
                    if 0
                        output_flat = plot_wlan(chanEst, [], 'spectral flatness', cfg, [], {0306, [signalFormat+", Pkt"+string(nPkt)+", User"+string(nUser)], 'Spectral Flatness'});
                    else
                        output_flat = PLOT_WLAN.spectralFlatness([], chanEst, signalFormat,...
                            {[0306, numSTSs, numPackets, nPkt], [signalFormat+dispLgnd], sprintf('Spectral Flatness'+".User%d.Ant%d.Pkt%d", nUser, nSTS, nPkt)});
                    end
                    isFlatPass(:,nPkt,nUser) = output_flat.isFlatPass;
                    flatDeviation{:,nPkt,nUser} = output_flat.flatDeviation;
                end

                % Extract Data field and perform the demodulation/noise
                % estimation and
                switch class(cfg)
                    case {'wlanEHTMUConfig'}
                        ehtData = rxPacket(ind.EHTData(1):ind.EHTData(2),:); % extract data signal
                        ehtDataDemod = wlanEHTDemodulate(ehtData,'EHT-Data',cfg,nUser); % demodulated data signal for noise estimation
                        ehtDataDemod_phsPilotTracking = wlanEHTTrackPilotError(ehtDataDemod,chanEst,cfg,'EHT-Data',nUser); % Perform pilot phase tracking
                        nVarEst = wlanEHTDataNoiseEstimate(ehtDataDemod_phsPilotTracking(ofdmInfo_i.PilotIndices,:,:),pilotEst,cfg,nUser); % Perform noise estimation
                        ehtDataDemod_scs = ehtDataDemod_phsPilotTracking(ofdmInfo_i.DataIndices,:,:); % extract data subcarrier form demodulation data
                        chanEstData = chanEst(ofdmInfo_i.DataIndices,:,:); % get channel estimation data subcarriers
                        [eqSym{nUser,nPkt},csi] = wlanEHTEqualize(ehtDataDemod_scs,chanEstData,nVarEst,cfg,'EHT-Data',nUser); % perform the equalization for equalized symbols
                        rxPSDU = wlanEHTDataBitRecover(eqSym{nUser,nPkt},nVarEst,csi,cfg,nUser); % recover data bit (rxPSDU) from equalized symbols
                    case {'wlanHESUConfig','wlanHEMUConfig'}
                        heData = rxPacket(ind.HEData(1):ind.HEData(2),:); % extract data signal
                        heDataDemod = wlanHEDemodulate(heData,'HE-Data',cfg,nUser); % demodulated data signal for noise estimation
                        nVarEst = wlanHEDataNoiseEstimate(heDataDemod(ofdmInfo_i.PilotIndices,:,:),pilotEst,cfg,nUser); % Perform noise estimation
                        heDataDemod_phsPilotTracking = wlanHETrackPilotError(heDataDemod,chanEst,cfg,'HE-Data',nUser); % Perform pilot phase tracking
                        heDataDemod_scs = heDataDemod_phsPilotTracking(ofdmInfo_i.DataIndices,:,:); % extract data subcarrier form demodulation data
                        chanEstData = chanEst(ofdmInfo_i.DataIndices,:,:); % get channel estimation data subcarriers
                        [eqSym{nUser,nPkt},csi] = wlanHEEqualize(heDataDemod_scs,chanEstData,nVarEst,cfg,'HE-Data',nUser); % perform the equalization for equalized symbols
                        if 0
                            PLOT_Constellation(eqSym{nUser,nPkt},[],12345)
                        end
                        if 0
                            rxPSDU = wlanHEDataBitRecover(eqSym{nUser,nPkt},nVarEst,csi,cfg,nUser); % recover data bit (rxPSDU) from equalized symbols
                        else
                            rxPSDU = wlanHEDataBitRecover(eqSym{nUser,nPkt},nVarEst,csi,cfg); % recover data bit (rxPSDU) from equalized symbols
                        end
                    case {'wlanVHTConfig'}
                        rxData = rxPacket(ind.VHTData(1):ind.VHTData(2),:);
                        % dataDemod = wlanVHTLTFDemodulate(rxData,cfg);
                        nVarEst = vhtNoiseEstimate(rxData,pilotEst,cfg);

                        % recover data bit (rxPSDU) from Data signal
                        % EqualizationMethod, Equalization method:
                        % 'MMSE' — The receiver uses a minimum mean-square error equalizer.
                        % 'ZF' — The receiver uses a zero-forcing equalizer.
                        % PilotPhaseTracking, Pilot phase tracking:
                        % 'PreEQ' — Enable pilot phase tracking, which the function performs before any equalization operation.
                        % 'None' — Disable pilot phase tracking.
                        % LDPCDecodingMethod, LDPC decoding algorithm, specified as the comma-separated pair consisting of 'LDPCDecodingMethod' and one of these values.
                        % 'bp' — Use the belief propagation (BP) decoding algorithm. For more information, see Belief Propagation Decoding.
                        % 'layered-bp' — Use the layered BP decoding algorithm, suitable for quasi-cyclic parity check matrices (PCMs). For more information, see Layered Belief Propagation Decoding.
                        % 'norm-min-sum' — Use the layered BP decoding algorithm with the normalized min-sum approximation. For more information, see Normalized Min-Sum Decoding.
                        % 'offset-min-sum' — Use the layered BP decoding algorithm with the offset min-sum approximation. For more information, see Offset Min-Sum Decoding.
                        if 1 % tx
                            [rxPSDU, rxSIGBCRC, eqSym{nUser,nPkt}] = wlanVHTDataRecover(rxData,chanEst,nVarEst,cfg,...
                                'EqualizationMethod','ZF','PilotPhaseTracking','PreEQ','LDPCDecodingMethod','norm-min-sum');
                        else % rx
                            [rxPSDU, rxSIGBCRC, eqSym{nUser,nPkt}] = wlanVHTDataRecover(rxData,chanEst,nVarEst,cfg,...
                                'LDPCDecodingMethod','norm-min-sum');
                        end
                    case {'wlanHTConfig'}
                        rxData = rxPacket(ind.HTData(1):ind.HTData(2),:);
                        nVarEst = htNoiseEstimate(rxData,chanEst,cfg);
                        if 1 % tx
                            [rxPSDU, eqSym{nUser,nPkt}] = wlanHTDataRecover(rxData,chanEst,nVarEst,cfg,...
                                'EqualizationMethod','ZF','PilotPhaseTracking','PreEQ','LDPCDecodingMethod','norm-min-sum'); % Use zero forcing algorithm for equalization
                        else % rx
                            [rxPSDU, eqSym{nUser,nPkt}] = wlanHTDataRecover(rxData,chanEst,nVarEst,cfg,...
                                'LDPCDecodingMethod','norm-min-sum'); % Use zero forcing algorithm for equalization
                        end
                end

                if 1  % 2024-03-07, Use plot_wlan instead of ehtTxEVMConstellationPlots
                    if 0
                        output_evm = plot_wlan(eqSym{nUser,nPkt}, [], 'evm', cfg, [] , {0307,[signalFormat+", Pkt"+string(nPkt)+", User"+string(nUser)]});
                    else
                        output_evm = PLOT_WLAN.evm([],eqSym{nUser,nPkt},cfg,...
                            {[0307, numSTSs, numPackets, nPkt], [signalFormat+dispLgnd], sprintf('EVM'+".User%d.Ant%d.Pkt%d", nUser, nSTS, nPkt)});
                    end
                    rmsEVM(:,nPkt,nUser) = output_evm.evmRMS;
                    evmPerSC{:,nPkt,nUser} =  output_evm.evmPerSc;
                    isEVMPass(:,nPkt,nUser) = output_evm.isEVMPass;
                end

                % check decoded
                idx_psduLen = (1:psduLen(nUser))' + (nPkt-1) * psduLen(nUser);
                if nPkt<=numPackets
                    if 1
                        % 2024-03-25, Determine if any bits are in error, i.e. a packet error
                        packetError = any(biterr(txPSDU{nUser}(idx_psduLen),rxPSDU));
                        isDecodePass(nUser,nPkt) = ~packetError;
                        numPacketErrors = numPacketErrors + packetError;
                        packetErrorRatePercent(nUser) = numPacketErrors/(numPackets)*100; % 2024-03-25, Calculate packet error rate (PER)
                    end
                end
            end
        end % nSTS
    end % nUser end

    % Store the offset of each packet within the waveform
    pktOffsetStore(nPkt) = pktOffset;

    % Increment waveform offset and search remaining waveform for a packet
    searchOffset = pktOffset+pktLength+minPktLen;
    nPkt = nPkt+1;
end % nPkt end

if isTxMeasurement % export - output
    output.isEVMPass = isEVMPass;
    output.isFlatPass = isFlatPass;
    output.isMaskPass = isMaskPass;
    output.evmRMS = rmsEVM;
    output.evmPerSc = evmPerSC;
    output.flatDeviation = flatDeviation;
    output.evmSpec = output_evm.evmSpec;
    output.flatSpec = output_flat.flatSpec;
    output.maskSpec = output_mask.maskSpec;
    output.outputPowerPerPktsDbm = outputPowerPerPktsDbm;
    % output.aclr = output_mask.aclr;
else
    % Calculate average input power per antenna
    output = plot_wlan(txWaveform(parms.usefulWvfm,:), [], 'sensitivity', cfg, [packetErrorRatePercent], {0326, [signalFormat+ ", User"+string(nUser)]});
    output.evmRMS = rmsEVM;
    output.evmPerSc = evmPerSC;
end
end

function [rxPacket, pktOffset, coarseFreqOff, fineFreqOff] = timingSync(rxWaveform,cfg,overSamplingFreq,searchOffset)
isSetToZero = 0;
% initialize
rxPacket = [];
pktOffset =[];
coarseFreqOff = [];
fineFreqOff = [];
cbw = cfg.ChannelBandwidth;
ind = wlanFieldIndices(cfg);
switch class(cfg)
    case {'wlanEHTMUConfig'}
        pktLength = double(ind.EHTData(2));
    case {'wlanHESUConfig','wlanHEMUConfig'}
        pktLength = double(ind.HEData(2));
    case {'wlanVHTConfig'}
        pktLength = double(ind.VHTData(2));
    case {'wlanHTConfig'}
        pktLength = double(ind.HTData(2));
end
if 0
    fs = ofdmInfo.SampleRate; % bug
else
    fs = overSamplingFreq; % 2024-03-14, timingSync: overSamplingFreq, use over sampling freq to improve the accuracy
end
rxWaveformLength = size(rxWaveform,1);

if ~exist('searchOffset','var')||isempty(searchOffset)
    searchOffset = 0;
end

% Detect packet and determine coarse packet offset
if 1
    pktOffset = wlanPacketDetect(rxWaveform,cbw,searchOffset);
else % 2024-03-13, modify wlanPacketDetect
    pktOffset = wlanPacketDetect(rxWaveform,cbw);
end
pktOffset = searchOffset+pktOffset; % Packet offset from start of the waveform
% Skip packet if L-STF is empty
if isempty(pktOffset) || (pktOffset<0) || ...
        ((pktOffset+ind.LSIG(2))>rxWaveformLength)
    return;
end

if 0
    % Scale the waveform based on L-STF power (AGC)
    gain = 1./(sqrt(mean(rxLSTF.*conj(rxLSTF))));
    rx = rx.*gain;
end

% Extract non-HT field(L-STF, L-LTF and L-SIG) and perform coarse frequency offset correction
% L-STF has good correlation properties, it is used for start-of-packet detection,
% for coarse frequency correction, and for setting the AGC.
nonht = rxWaveform(pktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
coarseFreqOff = wlanCoarseCFOEstimate(nonht,cbw)*(~isSetToZero);
% coarseFreqOff = 0;
nonht = frequencyOffset(nonht,fs,-coarseFreqOff);

% Extract the L-LTF of non-HT field and determine fine packet offset
lltfOffset = wlanSymbolTimingEstimate(nonht,cbw); % Channel estimation, fine frequency offset estimation, and fine symbol timing offset estimation rely on the L-LTF.
pktOffset = pktOffset+lltfOffset; % Determin packet offset

% If offset is outside the bounds of the waveform, then skip samples
% and continue searching within remainder of the waveform
if (pktOffset<0) || ((pktOffset+pktLength)>rxWaveformLength)
    searchOffset = pktOffset+double(ind.LSTF(2))+1;
    return;
end

% Apply coarse frequency correction to the extracted packed
rxPacket = rxWaveform(pktOffset+(1:pktLength),:);
rxPacket = frequencyOffset(rxPacket,fs,-coarseFreqOff);

% Extract L-LTF field and perform fine frequency offset correction
lltf = rxPacket(ind.LLTF(1):ind.LLTF(2),:); % Extract L-LTF
fineFreqOff = wlanFineCFOEstimate(lltf,cbw)*(~isSetToZero);; % Channel estimation, fine frequency offset estimation, and fine symbol timing offset estimation rely on the L-LTF.
% fineFreqOff = 0;
rxPacket = frequencyOffset(rxPacket,fs,-fineFreqOff);
end


