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
### Generate WIFI Waveform
```matlab
% run
configWG = ovwrConfig(waveformGenWLAN,configParamsWG); % get wv configurations
passParams = waveformGenWLAN.runWG(configWG); % get wv parameters
wv = powerDbm(wv,'set',paOutputPower_target_dBm,'dBm',dlSlots); % set transmitter power
```
|**ACLR**  |**Timing**  |**CCDF**|
|:--:|:--:|:--:|
| ![Image 1](https://github.com/user-attachments/assets/904f4054-8a8c-43ca-927a-a884c6efb7a9) | ![Image 2](https://github.com/user-attachments/assets/c1d62df2-7b45-464a-b53c-3dcbb1cde898) |![Image 3](https://github.com/user-attachments/assets/d649dd1c-a2d0-4fe5-aea9-5bb376667daf) |
### Transmitter -  Apply Noise - Before PA
```matlab
isNoiseFloorBeforePA = 0; % set noise floor to signal before PA
n0_dBmHz_set = -150 % set noise psd in dBm/Hz
% run
[wvNoise, noise0] = rfSim_noise(wv, 'NoisePowerDbmHz', n0_dBmHz_set, fs); % create noise wv
```
### Transmitter - Characterize PA and Signal
```matlab
paLinearGain_dB = 34.8; % set pa linear gain in dB
% run multi-stages pa lut
[lut_mx1, LookupTable1] = paLookupTable_process('paData.xlsx','0614SIMCW1M',[-50,50],[],1); % get pa1 lut from excel
[lut_mx2, LookupTable2] = paLookupTable_process('paData.xlsx','0614SIMCW2M',[-50,50],[],1); % get pa2 lut from excel
[lut_mx3, LookupTable3] = paLookupTable_process('paData.xlsx','0614SIMCW3M',[-50,50],[],1); % get pa3 lut from excel
```
### Transmitter -  Generate PA Signal and Tuning to Target Output Power
```matlab
paOutputPower_tolerance_dB = 0.2; % set pa output power tolerance
% run
[paSignal, paSignal_info_stages, paSignal_stages] = rfSim_pa_stages(paInSignal,'RFSIM',paModel,dlSlots); % generate pa output wv
```
### Transmitter - Flatness - After PA
```matlab
paFlatnessMax_dB    = 2 % set max flatness
paFlatFreqMHz       = cbw/2/1e6 * [0:1:configParamsWG.osr]; % implement flat. freq. vector in MHz
paFlatMagnDb        = -paFlatnessMax_dB * [0:1:configParamsWG.osr] + paFlatnessMax_dB/2; % implement flat. mag. vector in Db
nTapsRange_paFlat   = 10e3; % set range of taps for flat. fir
% run
[paSignalFlat, b_paFlat] = rfSim_fltaness(paSignal, fs, paFlatMagnDb, paFlatFreqMHz, nTapsRange_paFlat); % set wv with flatness
```
### Transmitter - wlanDemodulation and Measurements
```matlab
% run
demod_finalSignal = wlanDemodulation(passParams.finalSignal, passParams, 'tx') % demodulate wv
```
### Transmitter - Summary




