%% 2024-09-12, Draft


classdef PLOT_WLAN < handle
    properties
        CarrierConfig = [];
        ReferenceImpedance = 1;
        SignalFormat = string;
        MCS = nan;
        NumPackets = 1;
        ScsKHz = nan;
        RbwKHz = nan;
    end

    methods (Access=public)
        function obj=PLOT_WLAN(carrierConfig, numPackets)
            % set parameters
            if nargin > 0 && ~isempty(carrierConfig)
                obj.CarrierConfig = carrierConfig;
            end
            if nargin > 1 && ~isempty(numPackets)
                obj.NumPackets = numPackets;
            end

            if ~isempty(obj.CarrierConfig)
                try % EHT, HE
                    numUsers = numel(carrierConfig.User);
                catch % HT, VHT
                    numUsers = 1;
                end

                for nUser=1:numUsers
                    switch class(obj.CarrierConfig)
                        case 'wlanEHTMUConfig'
                            obj.ScsKHz(nUser) = 78.125;
                            obj.SignalFormat(nUser) = strrep(carrierConfig.ChannelBandwidth, 'CBW', "EHT");
                            tDFT(nUser) = 12.8e-6;
                            ofdmInfoPerUser(nUser) = wlanEHTOFDMInfo('EHT-Data',carrierConfig,nUser);
                            obj.MCS(nUser) = carrierConfig.User{:}.MCS;
                        case {'wlanHESUConfig','wlanHEMUConfig'}
                            obj.ScsKHz(nUser) = 78.125;
                            obj.SignalFormat(nUser) = strrep(carrierConfig.ChannelBandwidth, 'CBW', "HE");
                            tDFT(nUser) = 12.8e-6;
                            ofdmInfoPerUser(nUser) = wlanHEOFDMInfo('HE-Data',carrierConfig,nUser);
                            try
                                MCSIndex(nUser) = carrierConfig.MCS;
                            catch
                                MCSIndex(nUser) = carrierConfig.User{:}.MCS;
                            end
                        case 'wlanVHTConfig'
                            obj.ScsKHz(nUser) = 312.5;
                            obj.SignalFormat(nUser) = strrep(carrierConfig.ChannelBandwidth, 'CBW', "VHT");
                            ofdmInfoPerUser(nUser) = wlanVHTOFDMInfo('VHT-Data',carrierConfig);
                            tDFT(nUser) = 3.2e-6;
                            MCSIndex(nUser) = carrierConfig.MCS;
                        case 'wlanHTConfig'
                            obj.ScsKHz(nUser) = 312.5;
                            obj.SignalFormat(nUser) = strrep(carrierConfig.ChannelBandwidth, 'CBW', "HT");
                            ofdmInfoPerUser(nUser) = wlanHTOFDMInfo('HT-Data',carrierConfig);
                            tDFT(nUser) = 3.2e-6;
                        otherwise
                            error('obj.CarrierConfig?!')
                    end
                end % nUser
            end

        end
    end % method - public

    methods (Static)
        %% Plot - Timing
        function varargout = timing(obj, data, sampleRateMHz, carrierConfig, fnum_lgnd_titl_color_note_cell)
            if ~exist('carrierConfig','var')||isempty(carrierConfig)
                cfg = obj.CarrierConfig;
            else
                cfg = carrierConfig;
            end
            PLOT_COMM.time([], data, sampleRateMHz, 'usec', fnum_lgnd_titl_color_note_cell, 'Abs');
            maxData = ceil(max(abs(data))/1) * 1 ;

            if ~isempty(cfg)
                % get config.
                switch class(cfg)
                    case 'wlanEHTMUConfig'
                        [PSDU_LENGTH,TXTIME,~,~,trc,info] = ehtPLMETxTimePrimative_ica(cfg);
                    case {'wlanHESUConfig','wlanHEMUConfig'}
                        [PSDU_LENGTH,TXTIME,~,trc,info] = hePLMETxTimePrimative_ica(cfg);
                    case 'wlanVHTConfig'
                        [PSDU_LENGTH,TXTIME,~,trc,info] = hePLMETxTimePrimative_ica(cfg);
                    case 'wlanHTConfig'
                    otherwise
                        error('cfg!?')
                end

                % get symbols
                NSYM = info.commonCodingParams.NSYM;
                try
                    NEHTSIG = info.NEHTSIG;
                    NEHTLTF = info.NEHTLTF;
                end
                try
                    NHELTF = info.NHELTF;
                    NHESIGB = info.NHESIGB;
                    Nma = info.Nma;
                end
                SignalExtension = info.SignalExtension;

                % calculate timing
                u1 = 1e3;
                switch u1
                    case 1e3
                        t.unit = "us";
                end
                tLEG_PREAMBLE = (trc.TLSTF + trc.TLLTF + trc.TLSIG)/u1;
                t_DATA = NSYM*trc.TSYM/u1;
                switch packetFormat(cfg)
                    case 'HE-MU'
                        tHE_PREAMBLE = (trc.TRLSIG + trc.THESIGA + NHESIGB * trc.THESIGB + trc.THESTFNT + NHELTF * trc.THELTFSYM) / u1;
                        trc_series = [trc.TLSTF, trc.TLLTF, trc.TLSIG,...
                            trc.TRLSIG, trc.THESIGA, repmat(trc.THESIGB, 1, NHESIGB),...
                            trc.THESTFNT, repmat(trc.THELTFSYM, NHELTF),...
                            t_DATA*u1, repmat(trc.THELTFSYM, Nma * NHELTF), trc.TPE] / u1;
                    case 'HE-SU'
                        tHE_PREAMBLE = (trc.TRLSIG + trc.THESIGA + trc.THESTFNT + NHELTF * trc.THELTFSYM) / u1;
                        trc_series = [trc.TLSTF, trc.TLLTF, trc.TLSIG,...
                            trc.TRLSIG, trc.THESIGA, [],...
                            trc.THESTFNT, repmat(trc.THELTFSYM, NHELTF),...
                            t_DATA*u1, repmat(trc.THELTFSYM, Nma * NHELTF), trc.TPE] / u1;
                        txt_format = ["L-STF", "L-LTF", "L-SIG", "RL-SIG", "HE-SIG-A", [],...
                            "HE-STF", repmat("HE-LTF", NHELTF), "HE-Data", "PE"];
                    case 'HE-TB'
                        tHE_PREAMBLE = (trc.TRLSIG + trc.THESIGA + trc.THESTFT + NHELTF * trc.THELTFSYM) / u1;
                    case 'HE-EXT-SU'
                        tHE_PREAMBLE = (trc.TRLSIG + trc.THESIGAR + trc.THESTFNT + NHELTF * trc.THELTFSYM) / u1;
                    case 'EHT-MU'
                        tEHT_PREAMBLE = (trc.TRLSIG + trc.TUSIG + NEHTSIG * trc.TEHTSIG + trc.TEHTSTFNT + NEHTLTF * trc.TEHTLTFSYM) / u1; % Equation 36-97
                        trc_series = [trc.TLSTF, trc.TLLTF, trc.TLSIG, trc.TRLSIG, trc.TUSIG, repmat(trc.TEHTSIG, 1, NEHTSIG),...
                            trc.TEHTSTFNT, repmat(trc.TEHTLTFSYM, NEHTLTF), t_DATA*u1, trc.TPE] / u1;
                        txt_format = ["L-STF", "L-LTF", "L-SIG", "RL-SIG", "U-SIG", repmat("EHT-SIG", 1, NEHTSIG),...
                            "EHT-STF", repmat("EHT-LTF", NEHTLTF), "EHT-Data", "PE"];
                    case 'EHT-TB'
                        tEHT_PREAMBLE = (trc.TRLSIG + trc.TUSIG + trc.TEHTSTFT + NEHTLTF * trc.TEHTLTFSYM) / u1; % Equation 36-97
                        trc_series = [trc.TLSTF, trc.TLLTF, trc.TLSIG, trc.TRLSIG, trc.TUSIG, [],...
                            trc.TEHTSTFT, repmat(trc.TEHTLTFSYM, NEHTLTF), t_DATA*u1, trc.TPE] / u1;
                        txt_format = ["L-STF", "L-LTF", "L-SIG", "RL-SIG", "U-SIG", [],...
                            "EHT-STF", repmat("EHT-LTF", NEHTLTF), "EHT-Data", "PE"];
                end

                % get txTIME
                switch class(cfg)
                    case 'wlanEHTMUConfig'
                        txTIME = tLEG_PREAMBLE + tEHT_PREAMBLE + t_DATA + trc.TPE/u1 + SignalExtension; % TXTIME in ns. Equation 36-110
                    case {'wlanHESUConfig','wlanHEMUConfig'}
                        txTIME = tLEG_PREAMBLE + tHE_PREAMBLE + t_DATA + Nma * NHELTF * trc.THELTFSYM/u1 + trc.TPE/u1 + SignalExtension; % TXTIME in ns, IEEE P802.11ax/D4.1, Section 27.4.3, Equation 27-135
                end

                % plot timing format
                if trc.TPE==0
                    trc_series(end) = [];
                    txt_format(end) = [];
                end
                txt_series = txt_format + newline + string(trc_series) + t.unit;
                xIdx = [];
                yIdx = [];
                zIdx = [];
                xIdx_txt = [];
                yIdx_txt = 0.8;
                xIdx_xtick = 0;
                for k = 1:numel(trc_series)
                    if k==1
                        xIdx_xtick = 0;
                        tStart = 0;
                        tEnd = tStart + trc_series(k);
                    else
                        tStart = tEnd + 1 / (sampleRateMHz*1e6);
                        tEnd = tStart + trc_series(k);
                    end
                    xIdx = [xIdx, [tStart; tStart; tEnd; tEnd]];
                    yIdx = [yIdx, [0; 1; 1; 0]];
                    zIdx = [zIdx, [1; 1; 1; 1]];
                    xIdx_xtick = [xIdx_xtick, xIdx_xtick(end)+trc_series(k)];
                end

                patch(xIdx, maxData*yIdx, 0.1*zIdx, 'FaceAlpha', 0.1, 'DisplayName', "PPDU Format")
                xticks(xIdx_xtick);

                xIdx_txt = xIdx_xtick(1:end-1) + diff(xIdx_xtick)/2;
                text(xIdx_txt, repmat(0.8,size(xIdx_txt)), txt_series(1:numel(trc_series)), "FontWeight", 'normal');

                % return
                t.tPREAMBLE_NonHT   = tLEG_PREAMBLE;
                switch class(cfg)
                    case 'wlanEHTMUConfig'
                        t.tPREAMBLE    = tEHT_PREAMBLE;
                        t.tDATA        = t_DATA;
                    case {'wlanHESUConfig','wlanHEMUConfig'}
                        t.tPREAMBLE    = tHE_PREAMBLE;
                        t.tDATA        = t_DATA;
                end
                t.txTIME           = txTIME;
                varargout{1} = t;
            end
        end
        %% Plot - EVM
        function varargout = evm(obj, eqSymbPerPktUserSTS, carrierConfig, fnum_lgnd_titl_color_note_cell)
            [fnum, dispLgnd, dispTitle, color, pltType] = PLOT_WLAN.getPlotInfo(fnum_lgnd_titl_color_note_cell);

            % get carrier configuration
            cfg = PLOT_WLAN.getGarrierConfiguration(carrierConfig);

            % plot constellation
            mcs_spec_tab = wlan_spec('mcs',[],{"MCS"},{cfg.MCS});
            evm_spec_tab = wlan_spec('evm',[],{"MCS"},{cfg.MCS});
            evmVal = PLOT_COMM.evmConstellation(PLOT_COMM, eqSymbPerPktUserSTS, mcs_spec_tab.Modulation, 'dB', fnum_lgnd_titl_color_note_cell);
            isEVMRmsPass = (evm_spec_tab.EVM > evmVal.evmRMS);
            if isEVMRmsPass
                dispLgnd_pass = sprintf(".EVM %.2fdB.Pass", round(evmVal.evmRMS,2));
            else
                dispLgnd_pass = sprintf(".EVM %.2fdB.Fail", round(evmVal.evmRMS,2));
            end
            lgd = legend;
            lgd.String{end-1} = lgd.String{end-1}+dispLgnd_pass; % add pass info.
            if 0
                legend('Location', 'eastoutside');
            elseif 1*1
                legend('Location', 'northoutside');
            end

            % plot evm vs subcarrier
            evmPerSc = sqrt(sum(abs(evmVal.errorVec).^2, 2)) / size(evmVal.errorVec, 2);
            unit = 'dB';
            switch unit
                case '%'
                    evmPerSc = evmPerSc * 100;
                case 'dB'
                    evmPerSc = 20*log10(evmPerSc);
            end

            scsIdx = cfg.OfdmInfoPerUser.ActiveFrequencyIndices(cfg.OfdmInfoPerUser.DataIndices);
            nFFT = cfg.OfdmInfoPerUser.FFTLength;
            scsIdx_nFFT = (1:nFFT)' - (nFFT/2+1);
            
            figure(fnum(1)*10)
            if numel(fnum)==4
                hAx = subplot(fnum(2),fnum(3),fnum(4));
            else
                hAx = subplot(fnum(2),fnum(3),fnum(4));
            end
            hAx.YLimMode = 'auto';
            evmPerSc_nFFT = NaN(numel(scsIdx_nFFT),1);
            evmPerSc_nFFT(scsIdx+nFFT/2+1,1) = evmPerSc; % Convert the unit of EVM to dB for plotting
            plot(scsIdx_nFFT,evmPerSc_nFFT,'.-','DisplayName',dispLgnd), hold on, grid minor
            xlabel('Subcarrier Index'), ylabel(['EVM (', unit,')']), legend
            title(dispTitle)
            switch get(gca, 'XMinorGrid')
                case 'on'
                case 'off'
                    grid minor
            end

            if 1 % align the YLim
                allAxes = findobj(gcf, 'Type', 'axes');
                setYLim(1) = min([allAxes(end).YLim, hAx.YLim]);
                setYLim(2) = max([allAxes(end).YLim, hAx.YLim]);
                set(allAxes, 'YLim', setYLim);
            end

            if 0
                legend('Location', 'eastoutside');
            elseif 1*1
                legend('Location', 'northoutside');
                % legend('Location', 'northoutside');
            end
            
            % return
            varargout{1}.evmRMS = evmVal.evmRMS;
            varargout{1}.evmPerSc = evmPerSc;
            varargout{1}.isEVMPass = isEVMRmsPass;
            varargout{1}.evmSpec = evm_spec_tab;

        end

        %% Plot - Spectral Flatness
        function varargout = spectralFlatness(obj, chEstDataPerPktUserSTS, signalFormat_str, fnum_lgnd_titl_color_note_cell)
            [fnum, dispLgnd, dispTitle, color, pltType] = PLOT_WLAN.getPlotInfo(fnum_lgnd_titl_color_note_cell);

            if ~exist('signalFormat_str','var')||isempty(signalFormat_str)
                signalFormat_str = obj.SignalFormat;
            end
            % Validate number of subcarriers
            chanEstData = chEstDataPerPktUserSTS;
            [numSCs,numSTSs,numRxAnts] = size(chanEstData);
            if numRxAnts<numSTSs
                numSTSs = numRxAnts;
            end
            chanEstData = chanEstData(:,1:numSTSs,1:numSTSs);
            t = reshape(repelem(eye(numSTSs),numSCs,1),numSCs,numSTSs,[]);
            estUse = reshape(chanEstData(t==1),numSCs,[]);

            % get spec.
            flat_spec_tab = wlan_spec('spectral flatness', string(signalFormat_str));
            % get requirement
            activeScsIdx = sort([flat_spec_tab.AvgScsIdx{:}';flat_spec_tab.TestScsIdxEdge{:}'],1,'ascend');
            activeScsIdx_vec = [];
            for k=1:(numel(activeScsIdx)/2)
                activeScsIdx_vec = [activeScsIdx_vec; [activeScsIdx(2*k-1):activeScsIdx(2*k)]'];
            end
            avgScsIdx_cent = sort([flat_spec_tab.AvgScsIdx{:}'],1,'ascend');
            avgScsIdx_cent_vec = [];
            for k=1:(numel(avgScsIdx_cent)/2)
                avgScsIdx_cent_vec = [avgScsIdx_cent_vec; [avgScsIdx_cent(2*k-1):avgScsIdx_cent(2*k)]'];
            end
            testScsIdx_cent = sort([flat_spec_tab.TestScsIdx{:}'],1,'ascend');
            testScsIdx_cent_vec = [];
            for k=1:(numel(testScsIdx_cent)/2)
                testScsIdx_cent_vec = [testScsIdx_cent_vec; [testScsIdx_cent(2*k-1):testScsIdx_cent(2*k)]'];
            end
            testScsIdx_edge = sort([flat_spec_tab.TestScsIdxEdge{:}'],1,'ascend');
            testScsIdx_edge_vec = [];
            for k=1:(numel(testScsIdx_edge)/2)
                testScsIdx_edge_vec = [testScsIdx_edge_vec; [testScsIdx_edge(2*k-1):testScsIdx_edge(2*k)]'];
            end

            % get Ndft
            Ndft = flat_spec_tab.ChannelBandwidth_MHz*1e6 / (flat_spec_tab.SubcarrierSpacing_kHz*1e3);

            % Store channel estimates in subcarrier locations with nulls
            estFullFFT = complex(zeros(Ndft,numSTSs));
            estFullFFT(activeScsIdx_vec+Ndft/2+1,:) = estUse;

            % Calculate average power for the Average Subcarrier Indices of channel estimate
            avgChanEst_avgScs = 20*log10(mean(abs(estFullFFT(avgScsIdx_cent_vec+Ndft/2+1,:))));

            % Calculate psd power of Test Subcarrier Indices of channel estimate
            avgChanEst_testScs_1 = 20*log10((abs(estFullFFT(testScsIdx_cent_vec+Ndft/2+1,:))));
            avgChanEst_testScs_2 = 20*log10((abs(estFullFFT(testScsIdx_edge_vec+Ndft/2+1,:))));
            chanEstDeviation{1} = minus(avgChanEst_testScs_1, avgChanEst_avgScs);
            chanEstDeviation{2} = minus(avgChanEst_testScs_2, avgChanEst_avgScs);

            % check
            isFlatPass_l1 = ~any(chanEstDeviation{1}(:)<(min(flat_spec_tab.MaxDev_dB)));
            isFlatPass_u1 = ~any(chanEstDeviation{1}(:)>(max(flat_spec_tab.MaxDev_dB)));
            isFlatPass_l2 = ~any(chanEstDeviation{2}(:)<(min(flat_spec_tab.MaxDev_Edg_dB)));
            isFlatPass_u2 = ~any(chanEstDeviation{2}(:)>(max(flat_spec_tab.MaxDev_Edg_dB)));
            isFlatPass = all([isFlatPass_l1,isFlatPass_u1,isFlatPass_l2,isFlatPass_u2]);

            % plot flatness
            if ~isempty(fnum)
                chanEstDeviation_vec = nan(Ndft,numSTSs);
                chanEstDeviation_vec(testScsIdx_cent_vec+Ndft/2+1,:) = chanEstDeviation{1};
                chanEstDeviation_vec(testScsIdx_edge_vec+Ndft/2+1,:) = chanEstDeviation{2};
                chanEstDeviation_vec = rmmissing(chanEstDeviation_vec);
                figure(fnum(1))
                if numel(fnum)==4
                    hAx = subplot(fnum(2),fnum(3),fnum(4));
                else
                    hAx = subplot(1,1,1);
                end
                % get gca
                % h = gca(fnum(1));
                isCombinedMask = 1;
                hChid = hAx.Children;
                if numel(hChid)>=2 % delete mask
                    delete(hChid(1));
                    if ~isCombinedMask
                        delete(hChid(2));
                    end
                end

                if isFlatPass
                    dispLgnd_pass = '.Pass';
                else
                    dispLgnd_pass = '.Fail';
                end
                plot(activeScsIdx_vec, chanEstDeviation_vec, 'o', 'DisplayName', [dispLgnd+dispLgnd_pass]), hold on, legend
                title(dispTitle)
                xlabel('Subcarrier Index'), ylabel('Deviation (dB)')

                % plot mask
                dispLgnd_mask_l = ['Spectral Flatness Low'];
                specDeviation_low_vec = nan(Ndft,1);
                specDeviation_low_vec(testScsIdx_cent_vec+Ndft/2+1) = min(flat_spec_tab.MaxDev_dB);
                specDeviation_low_vec(testScsIdx_edge_vec+Ndft/2+1) = min(flat_spec_tab.MaxDev_Edg_dB);
                if 1
                    specDeviation_low_vec = rmmissing(specDeviation_low_vec);
                else
                    specDeviation_low_vec = specDeviation_low_vec(~isnan(specDeviation_low_vec));
                end
                if ~isCombinedMask
                    plot(activeScsIdx_vec, specDeviation_low_vec, 'r.', 'DisplayName', dispLgnd_mask_l), hold on
                end
                dispLgnd_mask_u = ['Spectral Flatness Upper'];
                specDeviation_upp_vec = nan*ones(Ndft,1);
                specDeviation_upp_vec(testScsIdx_cent_vec+Ndft/2+1) = max(flat_spec_tab.MaxDev_dB);
                specDeviation_upp_vec(testScsIdx_edge_vec+Ndft/2+1) = max(flat_spec_tab.MaxDev_Edg_dB);
                if 1
                    specDeviation_upp_vec = rmmissing(specDeviation_upp_vec);
                else
                    specDeviation_upp_vec = specDeviation_upp_vec(~isnan(specDeviation_upp_vec));
                end
                if ~isCombinedMask
                    plot(activeScsIdx_vec, specDeviation_upp_vec, 'r.', 'DisplayName', dispLgnd_mask_u), hold on
                else
                    % scatter combined mask
                    scatter([activeScsIdx_vec(:);activeScsIdx_vec(:)], [specDeviation_low_vec(:);specDeviation_upp_vec(:)], 'r.', 'DisplayName', 'Mask'), hold on
                end
                switch get(gca, 'XMinorGrid')
                    case 'on'
                    case 'off'
                        grid minor
                end

                if 0
                    legend('Location', 'eastoutside');
                else
                    legend('Location', 'northoutside');
                end

                % return
                output.isFlatPass = isFlatPass;
                output.flatDeviation = chanEstDeviation;
                output.flatSpec = flat_spec_tab;
            end
            % end % s

            % return
            varargout{1} = output;


        end

        %% Plot - Spectral Mask
        function varargout = spectralMask(obj, dataPerPktUserSTS, signalFormat_str, sampleRateMHz, fnum_lgnd_titl_color_note_cell, numPackets_numUsers_numSSsPerUser, nPkt_nUser_nSSsPerUser)
            if ~exist('signalFormat_str','var')||isempty(signalFormat_str)
                signalFormat_str = obj.SignalFormat;
            end
            % if ~exist('numPackets','var')||isempty(numPackets)
            %     numPackets = obj.NumPackets;
            % end

            if ~exist('numPackets_numUsers_numSSsPerUser','var')||isempty(numPackets_numUsers_numSSsPerUser)
                numPackets = [];
            else
                numPackets = numPackets_numUsers_numSSsPerUser(1);
                numUsers = numPackets_numUsers_numSSsPerUser(2);
                numSSsPerUser = numPackets_numUsers_numSSsPerUser(3);
            end

            if ~exist('nPkt_nUser_nSSsPerUser','var')||isempty(nPkt_nUser_nSSsPerUser)
                nPkt = [];
            else
                nPkt = nPkt_nUser_nSSsPerUser(1); % Nth packet in NumPackets
                nUser = nPkt_nUser_nSSsPerUser(2); % Nth user in NumUsers
                nSTS = nPkt_nUser_nSSsPerUser(3); % Nth in NumSpaceTimeStreams
            end

            % n=0;
            isMaskPass = 1;
            % for u=1:numUsers
            signalFormat_chr = char(signalFormat_str);
            switch signalFormat_chr(1:2)
                case {'EH','HE'}
                    rbwKHz = 100;
                    vbwKHz = 7.5;
                case {'VH','HT'}
                    rbwKHz = 100;
                    vbwKHz = 30;
                otherwise
                    error('signalFormat_str !?')
            end
            [fnum, dispLgnd, dispTitle, color, pltType] = PLOT_WLAN.getPlotInfo(fnum_lgnd_titl_color_note_cell);
            fnum_lgnd_titl_color_note_cell_tmp = {fnum, dispLgnd, dispTitle, color, pltType};
            if numel(fnum)~=4
                fnum_tmp = [fnum(1), numSSsPerUser, numPackets, nPkt]; % add subplot
                fnum_lgnd_titl_color_note_cell_tmp{1} = fnum_tmp;
                fnum_lgnd_titl_color_note_cell_tmp{3} = sprintf(fnum_lgnd_titl_color_note_cell{3} +".User%d.Ant%d.Pkt%d", nUser, nSTS, nPkt);
                figure(fnum(1))
                hAx = subplot(fnum(2),fnum(3),fnum(4));
            elseif 0
                nPkt = fnum_lgnd_titl_color_note_cell_tmp{1}(end);
            else
                figure(fnum(1))
                hAx = subplot(fnum(2),fnum(3),fnum(4));
            end

            % get gca
            hChid = hAx.Children; % get hAx children
            if numel(hChid)>=2 % delete mask
                delete(hChid(1)); % no.1 is the mask
            end

            psd_data_dBr = [];
            if 0
                data = reshape(data(:,nSTS), [], numPackets);
                saVal = PLOT_COMM.sa([], data(:,p), fnum_lgnd_titl_color_note_cell_tmp, sampleRateMHz, {rbwKHz,vbwKHz}, [], []);
                powerDataPktDbm(s,p) = 10*log10(mean(abs(data(:,p)).^2)) + 30;
            else
                data = dataPerPktUserSTS;
                saVal = PLOT_COMM.sa([], data, fnum_lgnd_titl_color_note_cell_tmp, sampleRateMHz, {rbwKHz,vbwKHz}, [], 'dBr');
                powerDataPktDbm = 10*log10(mean(abs(data).^2)) + 30;
            end
            psd_data_dBr = [saVal.Plot.PlotData, psd_data_dBr];

            % plot mask
            % try
            %     signalFormat_str_p = signalFormat_str(p);
            % catch
            %     signalFormat_str_p = signalFormat_str;
            % end
            mask_spec_tab = wlan_spec('spectral mask', signalFormat_str);
            freqs = saVal.Plot.FreqData;
            mask_data = zeros(size(saVal.Plot.PsdData));
            mask_vec = [0 mask_spec_tab.SpectralMask_dBr];
            freqs_mask_vec = [mask_spec_tab.FreqOffsetRange_MHz] * 1e6;
            for m=1:numel(mask_vec)
                if m~=numel(mask_vec)
                    freqs_mask_lp = freqs_mask_vec(m) * 1;
                    freqs_mask_hp = freqs_mask_vec(m+1) * 1;
                    freqs_mask_ln = - freqs_mask_lp;
                    freqs_mask_hn = - freqs_mask_hp;
                    mask_h = mask_vec(m);
                    mask_l = mask_vec(m);
                    idx_p = find(freqs > freqs_mask_lp & freqs <= freqs_mask_hp);
                    idx_n = find(freqs < freqs_mask_ln & freqs >= freqs_mask_hn);
                    mask_data(idx_p) = linspace(mask_h, mask_l, numel(idx_p));
                    mask_data(idx_n) = linspace(mask_l, mask_h, numel(idx_n));
                else
                    freqs_mask_p = freqs_mask_vec(m) * 1;
                    freqs_mask_n = - freqs_mask_p;
                    idx_lim = find(freqs > freqs_mask_p | freqs < freqs_mask_n);
                    mask_lim = mask_vec(m);
                    mask_data(idx_lim) = mask_lim;
                end
            end

            hChid = hAx.Children; % update hAx children
            isMaskPass = all(mask_data >= psd_data_dBr);
            if isMaskPass
                hChid(1).DisplayName = fnum_lgnd_titl_color_note_cell_tmp{2} + ".Pass";
            else
                hChid(1).DisplayName = fnum_lgnd_titl_color_note_cell_tmp{2} + ".Fail";
            end
            plot(freqs/1e6, mask_data, 'LineWidth', 1.5, 'Color', 'r', 'DisplayName', 'Mask'); grid minor
            switch get(gca, 'XMinorGrid')
                case 'on'
                case 'off'
                    grid minor
            end

            if 1 % align the YLim
                allAxes = findobj(gcf, 'Type', 'axes');
                if diff(allAxes(1).YLim) > diff(allAxes(end).YLim)
                    setYLim = allAxes(1).YLim;
                else
                    setYLim = allAxes(end).YLim;
                end
                set(allAxes, 'YLim', setYLim);
            end

            if 0
                legend('Location', 'eastoutside');
            else
                legend('Location', 'northoutside');
            end

            % return
            output.isMaskPass = isMaskPass;
            output.outputPowerPerPktsDbm = powerDataPktDbm;
            output.maskSpec = mask_spec_tab;
            varargout{1} = output;

        end

        %% Plot - ACLR
        function varargout = aclr(obj, data, signalFormat_str, sampleRateMHz, rbwKHz_vbwKHz_cell, fnum_lgnd_titl_color_note_cell)
            if ~exist('signalFormat_str','var')||isempty(signalFormat_str)
                signalFormat_str = obj.SignalFormat;
            end
            vbwKHz = [];
            if isempty(rbwKHz_vbwKHz_cell) && ~isempty(obj.ScsKHz)
                rbwKHz = obj.ScsKHz;
            elseif iscell(rbwKHz_vbwKHz_cell)
                rbwKHz = rbwKHz_vbwKHz_cell{1};
                vbwKHz = rbwKHz_vbwKHz_cell{2};
            else
                rbwKHz = 10;
            end
            data_pkt = reshape(data(:,1), [], obj.NumPackets);

            spectral_spec_tab = wlan_spec('spectral mask', signalFormat_str);
            clear caPm
            caPm.carrierBwMHz   = spectral_spec_tab.FreqOffsetRange_MHz(1)*2;
            caPm.channelBwMHz   = spectral_spec_tab.ChannelBandwidth_MHz;

            aclrVal = PLOT_COMM.aclr([], data_pkt, fnum_lgnd_titl_color_note_cell, sampleRateMHz, [], {caPm.carrierBwMHz, caPm.channelBwMHz}, {1}, {rbwKHz,vbwKHz});
            varargout{1} = aclrVal;

        end

    end % methods - Static

    methods (Static)

        function cfg = getGarrierConfiguration(carrierConfig)
            try % EHT, HE
                numUsers = numel(carrierConfig.User);
            catch % HT, VHT
                numUsers = 1;
            end

            for nUser=1:numUsers
                switch class(carrierConfig)
                    case 'wlanEHTMUConfig'
                        cfg.ScsKHz(nUser) = 78.125;
                        cfg.SignalFormat(nUser) = strrep(carrierConfig.ChannelBandwidth, 'CBW', "EHT");
                        tDFT(nUser) = 12.8e-6;
                        cfg.OfdmInfoPerUser(nUser) = wlanEHTOFDMInfo('EHT-Data',carrierConfig,nUser);
                        cfg.MCS(nUser) = carrierConfig.User{:}.MCS;
                    case {'wlanHESUConfig','wlanHEMUConfig'}
                        cfg.ScsKHz(nUser) = 78.125;
                        cfg.SignalFormat(nUser) = strrep(carrierConfig.ChannelBandwidth, 'CBW', "HE");
                        tDFT(nUser) = 12.8e-6;
                        cfg.OfdmInfoPerUser(nUser) = wlanHEOFDMInfo('HE-Data',carrierConfig,nUser);
                        try
                            cfg.MCS(nUser) = carrierConfig.MCS;
                        catch
                            cfg.MCS(nUser) = carrierConfig.User{:}.MCS;
                        end
                    case 'wlanVHTConfig'
                        cfg.ScsKHz(nUser) = 312.5;
                        cfg.SignalFormat(nUser) = strrep(carrierConfig.ChannelBandwidth, 'CBW', "VHT");
                        cfg.OfdmInfoPerUser(nUser) = wlanVHTOFDMInfo('VHT-Data',carrierConfig);
                        tDFT(nUser) = 3.2e-6;
                        cfg.MCS(nUser) = carrierConfig.MCS;
                    case 'wlanHTConfig'
                        cfg.ScsKHz(nUser) = 312.5;
                        cfg.SignalFormat(nUser) = strrep(carrierConfig.ChannelBandwidth, 'CBW', "HT");
                        cfg.OfdmInfoPerUser(nUser) = wlanHTOFDMInfo('HT-Data',carrierConfig);
                        tDFT(nUser) = 3.2e-6;
                        cfg.MCS(nUser) = carrierConfig.MCS;
                end
            end % nUser

            % return
            cfg = struct(cfg);
        end

        function [fnum, dispLgnd, dispTitle, color, note] = getPlotInfo(fnum_lgnd_titl_color_note_cell)
            fnum = 0;
            dispLgnd = "";
            dispTitle = "";
            color = [];
            note = [];

            if nargin>0
                if ~iscell(fnum_lgnd_titl_color_note_cell)
                    fnum = fnum_lgnd_titl_color_note_cell;
                else
                    if 1
                        fnum = PLOT_WLAN.getFieldOrDefault(fnum_lgnd_titl_color_note_cell, 1, []);
                        dispLgnd = PLOT_WLAN.getFieldOrDefault(fnum_lgnd_titl_color_note_cell, 2, "");
                        if ~strcmpi(dispLgnd,"")
                            dispLgnd = strrep(string(fnum_lgnd_titl_color_note_cell{2}), "_", ".");
                        end
                        dispTitle = PLOT_WLAN.getFieldOrDefault(fnum_lgnd_titl_color_note_cell, 3, "");
                        if ~strcmpi(dispTitle,"")
                            dispTitle = strrep(string(fnum_lgnd_titl_color_note_cell{3}), "_", ".");
                        end
                        color = PLOT_WLAN.getFieldOrDefault(fnum_lgnd_titl_color_note_cell, 4, []);
                        note = PLOT_WLAN.getFieldOrDefault(fnum_lgnd_titl_color_note_cell, 5, []);
                    else
                        try
                            fnum = fnum_lgnd_titl_color_note_cell{1};
                        catch
                        end
                        try
                            dispLgnd = strrep(string(fnum_lgnd_titl_color_note_cell{2}), "_", ".");
                        catch
                        end
                        try
                            dispTitle = strrep(string(fnum_lgnd_titl_color_note_cell{3}), "_", ".");
                        catch
                        end
                        try
                            color = fnum_lgnd_titl_color_note_cell{4};
                        catch
                        end
                        try
                            note = fnum_lgnd_titl_color_note_cell{5};
                        catch
                        end
                    end
                end
            end
        end

        function value = getFieldOrDefault(input, fieldOrIndex, defaultValue)
            % getFieldOrDefault: Get value from a struct, vector, or cell array at specified field or index, or return default value
            %
            % Inputs:
            %   input - struct, vector, or cell array to extract value from
            %   fieldOrIndex - field name (for struct) or index (for vector/cell array)
            %   defaultValue - default value to return if field/index is not found
            %
            % Output:
            %   output - struct containing the field 'DEF' with the value from the input or default value

            try
                if isempty(input) || isempty(fieldOrIndex)
                    value = input;
                elseif isstruct(input)
                    % Handle struct input
                    if isfield(input, fieldOrIndex)
                        value = input.(fieldOrIndex);
                    else
                        value = defaultValue;
                    end
                elseif iscell(input)
                    % Handle cell array input
                    if fieldOrIndex > 0 && fieldOrIndex <= numel(input)
                        value = input{fieldOrIndex};
                    else
                        value = defaultValue;
                    end
                elseif isnumeric(input) || islogical(input)
                    % Handle numeric or logical vector input
                    if fieldOrIndex > 0 && fieldOrIndex <= numel(input)
                        value = input(fieldOrIndex);
                    else
                        value = defaultValue;
                    end
                else
                    % Unsupported input type
                    value = defaultValue;
                end

                % if isempty of value set to defaultValue
                if isempty(value)
                    value = defaultValue;
                end

            catch
                value = defaultValue;
            end
        end

    end % methods - Static (internal)
end % class

