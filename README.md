# ic_WLAN_TX
Generate and transmit WIFI signal.

### Set WIFI line-up
```matlab
Rohm        = 1 % set resistor in ohm
RbwKHz      = 10 % set resolution bw for aclr plot
folderPath  = 'C:\Documents\MATLAB\waveform';
isPaModel   = 0; % turn on/off pa model
```
```matlab
simulation_trx = 'tx'; % set simulation case as 'tx'
paModel_LUT_sheet = '0527SIMCW'; % set LUT sheet
paModel_methology = 'RFSIMPA'; % set pa model function
paModel_type = 'LUT'; % set pa's amam/ampm from LUT
isPaMultiStage = 1; % set pa as multi-stages
isPaFlatness = 1;
isPaFlatnessCompensateLoss = 0;
isMIMO = 0;
```
```matlab
paModel_LUT_sheet   = '0527SIMCW'  % set LUT sheet
paModel_methology   = 'RFSIMPA' % set pa model function
paModel_type        = 'LUT' % set pa's amam/ampm from LUT
isPaMultiStage      = 1 % set pa as multi-stages
isPaFlatness        = 1
isPaFlatnessCompensateLoss = 0
isMIMO = 0
```
### Set WIFI Waveform Configurations
```matlab
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
configParamsWG.GuardInterval = 0.8; % set wv guard interval
if strcmpi(signalFormat, 'HT20')
   configParamsWG.bbFilter_en  = 1;
end
configParamsWG.osr          = 4; % set wv over sampling ratio
configParamsWG.resampleFilter_method    = 'waveformGen';
configParamsWG.resampleFilter_en        = 0;
configParamsWG.NumPackets = 4; % set wv number of packets and
configParamsWG.IdleTime = 10e-6; % set wv idle time
% save waveform
configParamsWG.save_tvs     = 1;
configParamsWG.fileFormat   = "TXT_ADS";
configParamsWG.wlanBand     = "5GHz";
configParamsWG.note         = "";
```
