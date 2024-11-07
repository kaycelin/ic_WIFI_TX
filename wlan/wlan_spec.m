%% 2024-01-18, draft
% 2024-03-26, add spec_val

function [spec_tab, spec_val] = wlan_spec(wlanCase, rowName_str, varNames_cell, varValues_cell)

    if 0 % initi.
        wlanCase;
        if ~exist('rowName_str','var')||~isstring(rowName_str)
            rowName_str = [];
        end
        if ~exist('varNames_cell','var')||~exist('varValues_cell','var')
            varNames_cell = [];
            varValues_cell = [];
        elseif iscell(varNames_cell)&&iscell(varValues_cell)
            varNames_cell;
            varValues_cell;
        elseif ~iscell(varNames_cell)&&~iscell(varValues_cell)
            varNames_cell = mat2cell(varNames_cell,1);
            varValues_cell = mat2cell(varValues_cell,1);
        else
            varNames_cell = [];
            varValues_cell = [];
        end
    else
        wlanCase;
        if ~exist('rowName_str','var')||isempty(rowName_str)
            rowName_str = [];
        elseif ~isstring(rowName_str)
            rowName_str = string(rowName_str);
        end
        if ~exist('varNames_cell','var')||isempty(varNames_cell)
            varNames_cell = [];
        elseif ~iscell(varNames_cell)
            varNames_cell = mat2cell(varNames_cell,1);
        end
        if ~exist('varValues_cell','var')||isempty(varValues_cell)
            varValues_cell = [];
        elseif ~iscell(varValues_cell)
            varValues_cell = mat2cell(varValues_cell,1);
        end
    end

    switch string(wlanCase) % generate specification
        case {"adjacentChannelRejection", "acr", "nonAdjacentChannelRejection", "nacr", "Table 36-68"}
            clear spec
            spec_mcs_tab = wlan_spec("mcs");
            spec_tab = spec_mcs_tab(:,{'MCS','Modulation','CodingRate'});
            if 1 % get signalFormat_str
                rowName_char = char(rowName_str);
                signalFormat_str = string(rowName_char(isstrprop(rowName_char, 'alpha')));
            end
            if strcmpi(signalFormat_str,"EHT")
                % Table 36-68
                spec_tab.AdjacentChannelRejection_dB = [16 13 11 8 4 0 -1 -2 -7 -9 -12 -14 -17 -20 16 16]';
                spec_tab.NonAdjacentChannelRejection_dB = [32 29 27 24 20 16 15 14 9 7 4 2 -1 -4 32 32]';
            elseif strcmpi(signalFormat_str,"HE")
                % Table 27-52
                spec_tab.AdjacentChannelRejection_dB = [16 13 11 8 4 0 -1 -2 -7 -9 -12 -14 nan nan nan nan]';
                spec_tab.AdjacentChannelRejection_80x2_dB = spec_tab.AdjacentChannelRejection_dB - 3;
                spec_tab.NonAdjacentChannelRejection_dB = [32 29 27 24 20 16 15 14 9 7 4 2 nan nan nan nan]';
                spec_tab.NonAdjacentChannelRejection_80x2_dB = spec_tab.NonAdjacentChannelRejection_dB - 3;
            elseif strcmpi(signalFormat_str,"VHT")
                % Table 21-26
                spec_tab.AdjacentChannelRejection_dB = [16 13 11 8 4 0 -1 -2 -7 -9 nan nan nan nan nan nan]';
                spec_tab.NonAdjacentChannelRejection_dB = [32 29 27 24 20 16 15 14 9 7 nan nan nan nan nan nan]';
            elseif strcmpi(signalFormat_str,"HT")
                % Table 19-23
                spec_tab.AdjacentChannelRejection_dB = [16 13 11 8 4 0 -1 -2 nan nan nan nan nan nan nan nan]';
                spec_tab.NonAdjacentChannelRejection_dB = [32 29 27 24 20 16 15 14 nan nan nan nan nan nan nan nan]';
            else
                error('wlan_spec.set rowName_str using signal format')
            end
            % spec = table2struct(spec_tab);
        case {"sensitivity", "Table 36-67"}
            clear spec
            spec_mcs_tab = wlan_spec("mcs");
            spec_tab = spec_mcs_tab(:,{'MCS','Modulation','CodingRate'});
            spec_tab.Sensitivity_20MHz_dBm = [-82 -79 -77 -74 -70 -66 -65 -64 -59 -57 -54 -52 -49 -46 nan -82]';
            spec_tab.Sensitivity_40MHz_dBm = spec_tab.Sensitivity_20MHz_dBm + 3;
            spec_tab.Sensitivity_80MHz_dBm = spec_tab.Sensitivity_40MHz_dBm + 3;  spec_tab.Sensitivity_80MHz_dBm(end-1) = -78;
            spec_tab.Sensitivity_160MHz_dBm = spec_tab.Sensitivity_80MHz_dBm + 3;
            spec_tab.Sensitivity_320MHz_dBm = spec_tab.Sensitivity_160MHz_dBm + 3;
            % spec = table2struct(spec_tab);
        case {"evm", "Table 36-65"}
            clear spec
            spec_mcs_tab = wlan_spec("mcs");
            spec_tab = spec_mcs_tab(:,{'MCS','Modulation','CodingRate'});
            spec_tab.EVM = [-27 -27 -27 -27 -27 -27 -27 -27 -30 -32 -35 -35 -38 -38 NaN -27]';
            % spec = table2struct(spec_tab);
        case {"mcs","Table 36-71"}
            clear spec
            spec.MCS = (0:15)';
            spec.MCS_HT = nan(size(spec.MCS));
            spec.MCS_VHT = nan(size(spec.MCS));
            spec.MCS_HE = nan(size(spec.MCS));
            spec.MCS_EHT = nan(size(spec.MCS));
            spec.MCS_HT((0:7)+1) = (0:7)';
            spec.MCS_VHT((0:9)+1) = (0:9)';
            spec.MCS_HE((0:11)+1) = (0:11)';
            spec.MCS_EHT((0:15)+1) = (0:15)';
            spec.NumSpatialStream = ones(size(spec.MCS));
            spec.Modulation = ["BPSK"; "QPSK"; "QPSK"; "16QAM"; "16QAM";...
                "64QAM"; "64QAM"; "64QAM"; "256QAM"; "256QAM";...
                "1024QAM"; "1024QAM"; "4096QAM"; "4096QAM"; "BPSK-DCM"; "BPSK-DCM"];
            spec.CodingRate = ["1/2"; "1/2"; "3/4"; "1/2"; "3/4";...
                "2/3"; "3/4"; "5/6"; "3/4"; "5/6";...
                "3/4"; "5/6"; "3/4"; "5/6"; "1/2"; "1/2"];
            spec.NumCodedBits = log2([2; 4; 4; 16; 16;...
                64; 64; 64; 256; 256;....
                1024; 1024; 4096; 4096; 2; 2]);
            numScsPerSymb_26RU = [24*ones(14,1); NaN; 12];
            spec.NumDataPerSymb_26RU = numScsPerSymb_26RU;
            spec.NumCodedBitsPerSymb = spec.NumDataPerSymb_26RU.*spec.NumCodedBits;
            codingRate = [1/2; 1/2; 3/4; 1/2; 3/4;...
                2/3; 3/4; 5/6; 3/4; 5/6;...
                3/4; 5/6; 3/4; 5/6; 1/2; 1/2];

            tDFT_us = 12.8;
            if 1 % equations (36-122, 36-123)
                numDataPerPerSymb_RU26 = [24*ones(13+1,1); NaN; 12];
                numCodedBitsPerSymb_RU26 = spec.NumCodedBits .* numDataPerPerSymb_RU26;
                numDataBitsPerSymb_RU26 = numCodedBitsPerSymb_RU26 .* codingRate;
                spec.DataRateMbs_GI0P8us_RU26 = round(10*(numDataBitsPerSymb_RU26 ./ (tDFT_us + 0.8)))/10;
                spec.DataRateMbs_GI1P6us_RU26 = round(10*(numDataBitsPerSymb_RU26 ./ (tDFT_us + 1.6)))/10;
                spec.DataRateMbs_GI3P2us_RU26 = round(10*(numDataBitsPerSymb_RU26 ./ (tDFT_us + 3.2)))/10;

                numDataPerPerSymb_RU52 = [48*ones(13+1,1); NaN; 24];
                numCodedBitsPerSymb_RU52 = spec.NumCodedBits .* numDataPerPerSymb_RU52;
                numDataBitsPerSymb_RU52 = numCodedBitsPerSymb_RU52 .* codingRate;
                spec.DataRateMbs_GI0P8us_RU52 = round(10*(numDataBitsPerSymb_RU52 ./ (tDFT_us + 0.8)))/10;
                spec.DataRateMbs_GI1P6us_RU52 = round(10*(numDataBitsPerSymb_RU52 ./ (tDFT_us + 1.6)))/10;
                spec.DataRateMbs_GI3P2us_RU52 = round(10*(numDataBitsPerSymb_RU52 ./ (tDFT_us + 3.2)))/10;

                numDataPerPerSymb_RU106 = [102*ones(13+1,1); NaN; 51];
                numCodedBitsPerSymb_RU106 = spec.NumCodedBits .* numDataPerPerSymb_RU106;
                numDataBitsPerSymb_RU106 = numCodedBitsPerSymb_RU106 .* codingRate;
                spec.DataRateMbs_GI0P8us_RU106 = round(10*(numDataBitsPerSymb_RU106 ./ (tDFT_us + 0.8)))/10;
                spec.DataRateMbs_GI1P6us_RU106 = round(10*(numDataBitsPerSymb_RU106 ./ (tDFT_us + 1.6)))/10;
                spec.DataRateMbs_GI3P2us_RU106 = round(10*(numDataBitsPerSymb_RU106 ./ (tDFT_us + 3.2)))/10;

                numDataPerPerSymb_RU242 = [234*ones(13+1,1); NaN; 117];
                numCodedBitsPerSymb_RU242 = spec.NumCodedBits .* numDataPerPerSymb_RU242;
                numDataBitsPerSymb_RU242 = numCodedBitsPerSymb_RU242 .* codingRate;
                spec.DataRateMbs_GI0P8us_RU242 = round(10*(numDataBitsPerSymb_RU242 ./ (tDFT_us + 0.8)))/10;
                spec.DataRateMbs_GI1P6us_RU242 = round(10*(numDataBitsPerSymb_RU242 ./ (tDFT_us + 1.6)))/10;
                spec.DataRateMbs_GI3P2us_RU242 = round(10*(numDataBitsPerSymb_RU242 ./ (tDFT_us + 3.2)))/10;

                numDataPerPerSymb_RU484 = [468*ones(13+1,1); NaN; 234];
                numCodedBitsPerSymb_RU484 = spec.NumCodedBits .* numDataPerPerSymb_RU484;
                numDataBitsPerSymb_RU484 = numCodedBitsPerSymb_RU484 .* codingRate;
                spec.DataRateMbs_GI0P8us_RU484 = round(10*(numDataBitsPerSymb_RU484 ./ (tDFT_us + 0.8)))/10;
                spec.DataRateMbs_GI1P6us_RU484 = round(10*(numDataBitsPerSymb_RU484 ./ (tDFT_us + 1.6)))/10;
                spec.DataRateMbs_GI3P2us_RU484 = round(10*(numDataBitsPerSymb_RU484 ./ (tDFT_us + 3.2)))/10;

                numDataPerPerSymb_RU996 = [980*ones(13+1,1); NaN; 490];
                numCodedBitsPerSymb_RU996 = spec.NumCodedBits .* numDataPerPerSymb_RU996;
                numDataBitsPerSymb_RU996 = numCodedBitsPerSymb_RU996 .* codingRate;
                spec.DataRateMbs_GI0P8us_RU996 = round(10*(numDataBitsPerSymb_RU996 ./ (tDFT_us + 0.8)))/10;
                spec.DataRateMbs_GI1P6us_RU996 = round(10*(numDataBitsPerSymb_RU996 ./ (tDFT_us + 1.6)))/10;
                spec.DataRateMbs_GI3P2us_RU996 = round(10*(numDataBitsPerSymb_RU996 ./ (tDFT_us + 3.2)))/10;

                numDataPerPerSymb_RU996_2 = [1960*ones(13+1,1); NaN; 980];
                numCodedBitsPerSymb_RU996_2 = spec.NumCodedBits .* numDataPerPerSymb_RU996_2;
                numDataBitsPerSymb_RU996_2 = numCodedBitsPerSymb_RU996_2 .* codingRate;
                spec.DataRateMbs_GI0P8us_RU996_2 = round(10*(numDataBitsPerSymb_RU996_2 ./ (tDFT_us + 0.8)))/10;
                spec.DataRateMbs_GI1P6us_RU996_2 = round(10*(numDataBitsPerSymb_RU996_2 ./ (tDFT_us + 1.6)))/10;
                spec.DataRateMbs_GI3P2us_RU996_2 = round(10*(numDataBitsPerSymb_RU996_2 ./ (tDFT_us + 3.2)))/10;

                numDataPerPerSymb_RU996_4 = [3920*ones(13+1,1); NaN; 1960];
                numCodedBitsPerSymb_RU996_4 = spec.NumCodedBits .* numDataPerPerSymb_RU996_4;
                numDataBitsPerSymb_RU996_4 = numCodedBitsPerSymb_RU996_4 .* codingRate;
                spec.DataRateMbs_GI0P8us_RU996_4 = round(10*(numDataBitsPerSymb_RU996_4 ./ (tDFT_us + 0.8)))/10;
                spec.DataRateMbs_GI1P6us_RU996_4 = round(10*(numDataBitsPerSymb_RU996_4 ./ (tDFT_us + 1.6)))/10;
                spec.DataRateMbs_GI3P2us_RU996_4 = round(10*(numDataBitsPerSymb_RU996_4 ./ (tDFT_us + 3.2)))/10;
            end
            % % struct to talbe
            if 1 % Table 36-21 Subcarrier allocation related constants for RUs in an OFDMA EHT PPDU
                rus.RU26 = [24 2 26]';
                rus.RU52 = [48 4 52]';
                rus.RU106 = [102 4 106]';
                rus.RU242 = [234 8 242]';
                rus.RU484 = [468 16 484]';
                rus.RU996 = [980 16 996]';
                rus.RU996_2 = [1960 32 1992]';
                rus.RU996_4 = [3920 64 3984]';
                rus.Description = ["Total number of data subcarriers per RU", "Number of pilot subcarriers per RU", 'Total number of subcarriers per RU']';
                rus_rowName = ["NSDtot", "NSP", "NST"]';
                rus_tab = struct2table(rus);
                rus_tab.Properties.RowNames = rus_rowName;
            end
        case {"mask","spectral mask"}
            clear spec
            % sigFormat
            spec_row = ["EHT20";"EHT40";"EHT80"; "EHT160"; "EHT320";...
                "HE20"; "HE40"; "HE80"; "HE160";...
                "VHT20"; "VHT40"; "VHT80"; "VHT160";...
                "HT20"; "HT40"; "HT20_5G"; "HT40_5G"];
            spec.Format...
                = ["EHT"; "EHT"; "EHT"; "EHT"; "EHT";...
                "HE"; "HE"; "HE"; "HE";...
                "VHT"; "VHT"; "VHT"; "VHT";...
                "HT"; "HT"; "HT"; "HT"];
            spec.ChannelBandwidth_MHz...
                = [20; 40; 80; 160; 320;...
                20; 40; 80; 160;...
                20; 40; 80; 160;...
                20; 40; 20; 40];
            spec.FreqOffsetRange_MHz...
                = [[9.75 10.5 20 30]; [19.5 20.5 40 60]; [39.5 40.5 80 120]; [79.5 80.5 160 240]; [159.5 160.5 320 480];...
                [9.75 10.5 20 30]; [19.5 20.5 40 60]; [39.5 40.5 80 120]; [79.5 80.5 160 240];...
                [9 11 20 30]; [19 21 40 60]; [39 41 80 120]; [79 81 160 240];...
                [9 11 20 30]; [19 21 40 60]; [9 11 20 30]; [19 21 40 60]];
            spec.SpectralMask_dBr...
                = [[-20 -28 -40]; [-20 -28 -40]; [-20 -28 -40]; [-20 -28 -40]; [-20 -28 -40];...
                [-20 -28 -40]; [-20 -28 -40]; [-20 -28 -40]; [-20 -28 -40];...
                [-20 -28 -40]; [-20 -28 -40]; [-20 -28 -40]; [-20 -28 -40];...
                [-20 -28 -45]; [-20 -28 -45]; [-20 -28 -40]; [-20 -28 -40];];
            spec.SpectralMaskLimit_dBmPerHz_vs_Freq...
                = {[-53, -39]; [-56, -39]; [-39, -39]; [-39, -39]; [-39, -39];...
                [-53 -53]; [-56 -56]; [-59 -59]; [-59 -59];...
                -53; -56; -59; -59;...
                -53; -53; -53; -53};
            spec.FreqBandGHz...
                = {[2.4 5]; [2.4 5]; [5 6]; [5 6]; [5 6];...
                [2.4 5]; [2.4 5]; [5 6]; [5 6];...
                5; 5; 5; 5;...
                2.4; 2.4; 5; 5};
            spec.Note...
                = ["none"; "none"; "none"; "none"; "none";...
                "none"; "none"; "none"; "80+80MHz";...
                "none"; "none"; "none"; "80+80MHz";...
                "none"; "none"; "none"; "none"];
        case {"flatness","spectral flatness"}
            clear spec
            % sigFormat
            spec_row = ["EHT20";"EHT40";"EHT80"; "EHT160"; "EHT320";...
                "HE20"; "HE40"; "HE80"; "HE160";...
                "VHT20"; "VHT40"; "VHT80"; "VHT160";...
                "nonHT40Duplicate"; "nonHT80Duplicate"; "nonHT160Duplicate"; "nonHT320Duplicate";...
                "HT20"; "HT40"];
            spec.Format...
                = ["EHT"; "EHT"; "EHT"; "EHT"; "EHT";...
                "HE"; "HE"; "HE"; "HE";...
                "VHT"; "VHT"; "VHT"; "VHT";...
                "nonHT"; "nonHT"; "nonHT"; "nonHT";...
                "HT"; "HT"];
            spec.ChannelBandwidth_MHz...
                = [20; 40; 80; 160; 320;...
                20; 40; 80; 160;...
                20; 40; 80; 160;...
                40; 80; 160; 320;...
                20; 40];
            spec.SubcarrierSpacing_kHz...
                = [repmat(78.125, 5, 1);...
                repmat(78.125, 4, 1);...
                repmat(312.5, 4, 1);...
                repmat(312.5, 4, 1);...
                repmat(312.5, 2, 1)];

            avgScsIdx_cbw20_scs78 = [-84 -2 2 84];
            avgScsIdx_cbw40_scs78 = [-168 -3 3 168];
            avgScsIdx_cbw80_scs78 = [-344 -3 3 344];
            avgScsIdx_cbw160_scs78 = [-696 -515 -509 -12 12 509 515 696];
            avgScsIdx_cbw320_scs78 = [-1400 -1036, -1012 -515, -509 -12 12 509 515 1012 1036 1400];
            avgScsIdx_cbw20_scs312 = [-16 -1 1 16];
            avgScsIdx_cbw40_scs312 = [-42 -2 2 42];
            avgScsIdx_cbw80_scs312 = [-84 -2 2 84];
            avgScsIdx_cbw160_scs312 = [-172 -130 -126 -44 44 126 130 172];
            avgScsIdx_cbw40_scs312_nonHTdup = [-42 -33 -31 -6 6 31 33 42];
            avgScsIdx_cbw80_scs312_nonHTdup = [-84 -70 -58 -33 -31 -6 6 31 33 58 70 84];
            avgScsIdx_cbw160_scs312_nonHTdup = [-172 -161 -159 -134 -122 -97 -95 -70 -58 -44 44 58 70 95 97 122 134 159 161 172];
            avgScsIdx_cbw320_scs312_nonHTdup = [-348 -326 -314 -300 -212 -198 -186 -161 -159 -134 -122 -97 -95 -84,...
                84 95 97 122 134 159 161 186 198 212 300 314 326 348];

            spec.AvgScsIdx...
                = {avgScsIdx_cbw20_scs78; avgScsIdx_cbw40_scs78; avgScsIdx_cbw80_scs78; avgScsIdx_cbw160_scs78; avgScsIdx_cbw320_scs78;...
                avgScsIdx_cbw20_scs78; avgScsIdx_cbw40_scs78; avgScsIdx_cbw80_scs78; avgScsIdx_cbw160_scs78;...
                avgScsIdx_cbw20_scs312; avgScsIdx_cbw40_scs312; avgScsIdx_cbw80_scs312; avgScsIdx_cbw160_scs312;...
                avgScsIdx_cbw40_scs312_nonHTdup; avgScsIdx_cbw80_scs312_nonHTdup; avgScsIdx_cbw160_scs312_nonHTdup; avgScsIdx_cbw320_scs312_nonHTdup;...
                avgScsIdx_cbw20_scs312; avgScsIdx_cbw40_scs312};
            spec.TestScsIdx...
                = {avgScsIdx_cbw20_scs78; avgScsIdx_cbw40_scs78; avgScsIdx_cbw80_scs78; avgScsIdx_cbw160_scs78; avgScsIdx_cbw320_scs78;...
                avgScsIdx_cbw20_scs78; avgScsIdx_cbw40_scs78; avgScsIdx_cbw80_scs78; avgScsIdx_cbw160_scs78;...
                avgScsIdx_cbw20_scs312; avgScsIdx_cbw40_scs312; avgScsIdx_cbw80_scs312; avgScsIdx_cbw160_scs312;...
                avgScsIdx_cbw40_scs312_nonHTdup; avgScsIdx_cbw80_scs312_nonHTdup; avgScsIdx_cbw160_scs312_nonHTdup; avgScsIdx_cbw320_scs312_nonHTdup;...
                avgScsIdx_cbw20_scs312; avgScsIdx_cbw40_scs312};

            spec.MaxDev_dB...
                = [[-4 4]; [-4 4]; [-4 4]; [-4 4]; [-4 4];...
                [-4 4]; [-4 4]; [-4 4]; [-4 4];...
                [-4 4]; [-4 4]; [-4 4]; [-4 4];...
                [-4 4]; [-4 4]; [-4 4]; [-4 4];...
                [-4 4]; [-4 4];...
                ];
            spec.MaxDev_wiPreamblePunct_dB...
                = [[-6 4]; [-6 4]; [-6 4]; [-6 4]; [-6 4];...
                [nan nan]; [nan nan]; [nan nan]; [nan nan];...
                [nan nan]; [nan nan]; [nan nan]; [nan nan];...
                [nan nan]; [nan nan]; [-6 4]; [-6 4];...
                [nan nan]; [nan nan];...
                ];

            testScsIdxEdg_cbw20_scs78 = [-122 -85 85 122];
            testScsIdxEdg_cbw40_scs78 = [-244 -169 169 244];
            testScsIdxEdg_cbw80_scs78 = [-500 -345 345 500];
            testScsIdxEdg_cbw160_scs78 = [-1012 -697 697 1012];
            testScsIdxEdg_cbw320_scs78 = [-2036 -1539 -1533 -1401 1401 1533 1539 2036];
            testScsIdxEdg_cbw20_scs312 = [-28 -17 17 28];
            testScsIdxEdg_cbw40_scs312 = [-58 -43 43 58];
            testScsIdxEdg_cbw80_scs312 = [-122 -85 85 122];
            testScsIdxEdg_cbw160_scs312 = [-250 -173 -43 -6 6 43 173 250];

            testScsIdxEdg_cbw40_nonHTdup = testScsIdxEdg_cbw40_scs312;
            testScsIdxEdg_cbw80_nonHTdup = [-122 -97 -95 -85 85 95 97 122];
            testScsIdxEdg_cbw160_nonHTdup = [-250 -225 -223 -198 -186 -173 -43 -33 -31 -6 6 31 33 43 173 186 198 223 225 250];
            testScsIdxEdg_cbw320_nonHTdup =...
                [-506 -481 -479 -454 -442 -417 -415 -390 -378 -353 -351 -349 -299 -289 -287 -262,...
                -250 -225 -223 -213 -83 -70 -58 -33 -31 -6,...
                6 31 33 58 70 83 213 223 225 250,...
                262 287 289 299 349 351 353 378 390 415 417 442 454 479 481 506];

            spec.TestScsIdxEdge...
                = {testScsIdxEdg_cbw20_scs78; testScsIdxEdg_cbw40_scs78; testScsIdxEdg_cbw80_scs78; testScsIdxEdg_cbw160_scs78; testScsIdxEdg_cbw320_scs78;...
                testScsIdxEdg_cbw20_scs78; testScsIdxEdg_cbw40_scs78; testScsIdxEdg_cbw80_scs78; testScsIdxEdg_cbw160_scs78;...
                testScsIdxEdg_cbw20_scs312; testScsIdxEdg_cbw40_scs312; testScsIdxEdg_cbw80_scs312; testScsIdxEdg_cbw160_scs312;...
                testScsIdxEdg_cbw40_nonHTdup; testScsIdxEdg_cbw80_nonHTdup; testScsIdxEdg_cbw160_nonHTdup; testScsIdxEdg_cbw320_nonHTdup;...
                testScsIdxEdg_cbw20_scs312; testScsIdxEdg_cbw40_scs312};

            spec.MaxDev_Edg_dB...
                = [[-6 4]; [-6 4]; [-6 4]; [-6 4]; [-6 4];...
                [-6 4]; [-6 4]; [-6 4]; [-6 4];...
                [-6 4]; [-6 4]; [-6 4]; [-6 4];...
                [-6 4]; [-6 4]; [-6 4]; [-6 4];...
                [-6 4]; [-6 4];...
                ];
            spec.MaxDev_wiPreamblePunctEdg_dB...
                = [[-6 4]; [-6 4]; [-6 4]; [-6 4]; [-6 4];...
                [nan nan]; [nan nan]; [nan nan]; [nan nan];...
                [nan nan]; [nan nan]; [nan nan]; [nan nan];...
                [nan nan]; [-6 4]; [-6 4]; [-6 4];...
                [nan nan]; [nan nan];...
                ];
            spec.Default...
                = repmat({"woPreamablePuncturing", "Center"}, numel(spec_row), 1);
        otherwise
            error('check the wlanCase ?')
    end

    % struct to talbe
    try
        spec_tab;
    catch
        spec_tab = struct2table(spec);
    end
    try
        spec_tab.Properties.RowNames = spec_row;
    catch
        rowName_str = [];
    end
    [spec_tab, spec_val] = tabFilter(spec_tab, rowName_str, varNames_cell, varValues_cell);
end

%% functions
function [output_tab, output_val] = tabFilter(input_tab, rowName_str, varNames_cell, varValues_cell)

    if exist('rowName_str','var')&&~isempty(rowName_str)
        isRowFilter = 1;
    else
        isRowFilter = 0;
    end

    if exist('varNames_cell','var')&&~isempty(varNames_cell)
        isVarFilter = 1;
    else
        isVarFilter = 0;
    end
    if exist('varValues_cell','var')&&~isempty(varValues_cell)
        isVarSelect = 1;
    else
        isVarSelect = 0;
    end

    if isRowFilter
        output_tab = input_tab(rowName_str,:);
    else
        output_tab = input_tab;
    end

    varNames_cell_str = [];
    output_val = [];
    output_tab_2 = [];
    if isVarFilter
        for k=1:numel(varNames_cell)
            varNames_cell_str = [varNames_cell_str, string(varNames_cell{k})];
            % filter
            try
                output_tab_tmp = output_tab(:, varNames_cell_str(k));
            catch
                output_tab_tmp = output_tab(:, char(varNames_cell_str(k)));
            end
            % try
            %     output_tab_2 = horzcat(output_tab_2,output_tab_tmp);
            % catch
            %     output_tab_2 = [];
            % end
            % get value
            if isVarSelect && k<=numel(varValues_cell)
                if ~isempty(varValues_cell{k})
                    output_tab = output_tab(ismember(output_tab_tmp.Variables, varValues_cell{k}), :);
                else
                    output_val = [output_val, output_tab(:,varNames_cell_str(k)).Variables];
                end
            else
                output_tab_2 = horzcat(output_tab_2,output_tab_tmp);
                output_val = output_tab_tmp.Variables;
            end
        end
        
        if ~isempty(output_tab_2)
            % export output table
            output_tab = output_tab_2;
        end
    end


end

% function output_tab = tabFilter(input_tab, rowName_str, varNames_cell, varValues_cell)
%
%     if exist('rowName_str','var')&&~isempty(rowName_str)&&isstring(rowName_str)
%         isRowSelect = 1;
%     else
%         isRowSelect = 0;
%     end
%
%     if exist('varNames_cell','var')&&~isempty(varNames_cell)&&iscell(varNames_cell)...
%             &exist('varValues_cell','var')&&~isempty(varValues_cell)&&iscell(varValues_cell)
%         isVarSelect = 1;
%     else
%         isVarSelect = 0;
%     end
%
%     if isRowSelect
%         output_tab = input_tab(rowName_str,:);
%     else
%         output_tab = input_tab;
%     end
%
%     if isVarSelect
%         for k=1:numel(varNames_cell)
%             try
%                 output_tab_tmp = output_tab(:, varNames_cell{k});
%                 output_tab = output_tab(ismember(output_tab_tmp.Variables, varValues_cell{k}), :);
%             end
%         end
%     end
%
%     if isempty(output_tab)
%         error('tabFilter, check the inputs')
%     end
%
% end

