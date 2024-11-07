%% 2024-04-17, draft

function [phy_tab, phy_val] = wlan_phy(phyCase, signalFormat_str, rowName_str, varNames_cell, varValues_cell)

    if 1 % initi.
        phyCase;
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
        if ~exist('signalFormat_str','var')||isempty('signalFormat_str')
            signalFormat_str = [];
        elseif ~isstring(signalFormat_str)
            signalFormat_str = string(signalFormat_str);
        end
    end

    switch string(phyCase) % generate specification
        case {"sa", "Subcarrier allocation related constants"}
            if strcmpi(signalFormat_str,"EHT")
                % Table 36-19
                sa_rowName = ["CBW20","CBW40","CBW80(non-EHT-MCS14)","CBW80(EHT-MCS14)","CBW160","CBW320"];
                sa.NSD_total = [234,468,980,936,1960,3920]';
                sa.NSP = [8,16,16,32,32,64]';
                sa.NST = sa.NSD_total + sa.NSP;
                sa.NSR = [122,244,500,500,1012,2036]';
                sa.NDC = [3,5,5,23,23,23]';
                sa.NGuardLeft = [6,12,12,12,12,12]';
                sa.NGuardRight = [5,11,11,11,11,11]';
            elseif strcmpi(signalFormat_str,"HE")
                % Table 27-13
                sa_rowName = ["CBW20","CBW40","CBW80","CBW80+80","CBW160"];
                sa.NSD = ["See 27.5","See 27.5","See 27.5","See 27.5","See 27.5"]';
                sa.NSP = [8,16,16,32,32]';
                sa.NST = [242,484,996,996,1992]';
                sa.NSR = [122,244,500,500,1012]';
                sa.NSeg = [1,1,1,2,1]';
                sa.NDC = [3,5,5,5,23]';
            elseif strcmpi(signalFormat_str,"VHT")
                % Table 19-6
                sa_rowName = ["NonHT_CBW20","HT_CBW20","HT_CBW40","NonHTDuplicate_MCS32_CBW40"]';
                sa.NSD          = [48,52,108,48]';
                sa.NSP          = [4,4,6,4]';
                sa.NST          = [52,56,114,104]';
                sa.NSR          = [26,28,58,58]';
                sa.deltaF_kHz   = [312.5,312.5,312.5,312.5]';
                sa.tDFT_us      = 1./sa.deltaF_kHz*1e3;
                sa.tGI_us       = sa.tDFT_us/4;
                sa.tGI2_us      = sa.tDFT_us/2;
                sa.tGIS_us      = [nan,sa.tDFT_us(2)/8,sa.tDFT_us(3)/8,nan]';
                sa.tLSTF_us     = 10*sa.tDFT_us;
                sa.tHTGFSTF_us  = [nan,10*sa.tDFT_us(2)/4,10*sa.tDFT_us(3)/4,nan]';
                sa.tLLTF_us     = 2*sa.tDFT_us+sa.tGI2_us;
                sa.tSYM_us      = sa.tDFT_us + sa.tGI_us;
                sa.tSYMS_us     = sa.tDFT_us+sa.tGIS_us;
                sa.tLSIG_us     = sa.tSYM_us;
                sa.tHTSIG_us    = 2*sa.tSYM_us;
                sa.tHTSTF_us    = [nan,4,4,nan]';
                sa.tHTLTF1_us   = [nan,4,4,nan]';
                sa.tHTLTFs_us   = [nan,4,4,nan]';
            elseif strcmpi(signalFormat_str,"HT")
                % Table 17-5
                sa_rowName      = ["CBW20","CBW10","CBW5"];
                sa.NSD          = [48,48,48]';
                sa.NSP          = [4,4,4]';
                sa.NST          = sa.NSD + sa.NSP;
                sa.deltaF_kHz   = [20,10,5]'./64*1e3;
                sa.tFFT_us      = 1./sa.deltaF_kHz*1e3;
                sa.tGI_us       = sa.tFFT_us/4;
                sa.tGI2_us      = sa.tFFT_us/2;
                sa.tSIGNAL_us   = sa.tGI_us + sa.tFFT_us;
                sa.tSYM_us      = sa.tGI_us + sa.tFFT_us;
                sa.tSHORT_us    = 10 * sa.tFFT_us/4;
                sa.tLONG_us     = sa.tGI2_us + 2*sa.tFFT_us;
                sa.tPREAMBLE_us = (sa.tSHORT_us + sa.tLONG_us);
            end
            phy_tab = struct2table(sa);
            phy_tab.Properties.RowNames = sa_rowName;
        case {"tr", "timing-related", "Table 36-18"}
            clear phy

            guardInterval = varValues_cell{1}; % get guard interval setting
            encoderType   = varValues_cell{2};
            signalExtension = varValues_cell{3};
            guardInterval = 0.8; % get guard interval setting
            encoderType   = 'BCC';
            signalExtension = 0;
            if exist("encoderType",'var')&&~isempty(encoderType)
                encoderType = 'NaN'
            end

            phy_tab.deltaF_PreEHT_kHz   = 312.5;
            phy_tab.deltaF_EHT_kHz      = 78.125;
            phy_tab.tDFT_PreEHT_us      = 3.2;
            phy_tab.tDFT_EHT_us         = 12.8;
            tGI_PreEHT_us       = 0.8;
            tGI_LLFT_us         = 1.6;
            switch guardInterval
                case 0.8
                    tGI_Data_us = 0.8;
                case 1.6
                    tGI_Data_us = 1.6;
                case 3.2
                    tGI_Data_us = 3.2;
            end
            tGI_ehtLTF_us       = tGI_Data_us;
            tSYM_us             = tDFT_EHT_us + tGI_Data_us;
            t_LSTF_us           = 10 * tDFT_PreEHT_us / 4;
            t_LSTF_us           = 2 * tDFT_PreEHT_us + tGI_LLFT_us;
            t_LSIG_us           = 4;
            t_RLSIG_us          = 4;
            t_USIG_us           = 2 * 4;
            t_USIGR_us          = 4 * 4;
            t_ehtSIG_us         = tDFT_PreEHT_us + tGI_PreEHT_us;
            t_ehtSTFT_us        = 5 * 1.6;                                  % EHT-STF field duration for an EHT TB PPDU
            t_ehtSTFNT_us       = 5 * 0.8;                                  % EHT-STF field duration for an EHT MU PPDU
            switch tGI_ehtLTF_us
                case 0.8
                    t_ehtLTF_us = 3.2;                                      % Duration of each OFDM symbol without GI in the EHT-LTF field
                case 1.6
                    tGI_Data_us = 6.4;
                case 3.2
                    tGI_Data_us = 12.8;
            end
            t_ehtLTFSYM_us      = t_ehtLTF_us + tGI_ehtLTF_us;              % Duaration of each OFDM symbol including GI in the EHT-LTF field
            nService            = 16;                                       % Number of bits in the SERVICE field            
            switch encoderType
                case "BCC"
                    NTail    = 6;
                    NTail_u  = 6;
                case "LDPC"
                    NTail    = 0;
                    NTail_u  = 0;
                otherwise
            end
            tSYML_us           = 4;                                        % OFDM symbol duration including GI in the pre-EHT modulated fields
            tPE_us             = [0, 4, 8, 12, 16, 20];                    % Duration of the PE field
                   
            phy = table2struct(phy_tab);
    
    end

    % % struct to talbe
    % spec_tab = struct2table(pec);
    % try
    %     spec_tab.Properties.RowNames = spec_row;
    % end

    switch string(phyCase)
        case 'sensitivity'
        otherwise
            % input conditions
            [phy_tab, phy_val] = tabFilter(phy_tab, rowName_str, varNames_cell, varValues_cell);
    end
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
    if isVarFilter
        for k=1:numel(varNames_cell)
            varNames_cell_str = [varNames_cell_str, string(varNames_cell{k})];
            output_tab_tmp = output_tab(:, varNames_cell_str(k));
            if isVarSelect
                output_tab = output_tab(ismember(output_tab_tmp.Variables, varValues_cell{k}), :);
            end
        end
        output_val = output_tab(:,varNames_cell_str).Variables;     
    else
        output_val = [];
    end

end

