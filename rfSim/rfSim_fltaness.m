% 2024-03-15, Draft rfSim: RF Simulation
% 2024-07-23, Modified from rfSim_fltaness
% 2024-07-30, Change nTaps to nTaps_isFFT
% 2024-07-30, Add interpGainCurve
% 2024-10-17, Change gainDb to gainDbVec, freqMHz to freqMHzVec
% 2024-10-17, Bug! for applyFilter='FFT'
% 2024-10-17, Normalize filter coefficients


function [signalOut, b] = rfSim_fltaness_symmetric(signal, sampleRate, gainDbVec, freqMHzVec, nTaps_isFFT, isGainNorm, fnum)
% Design a filter coefficients based frequency(freqMHzVec) and gain flatness(gainDbVec),
% and apply filter to signal
% The filter is symmetric type between nyquist 1 and nyquist 2
%
% Inputs:
%   signal
%   freqMHzVec         - Frequency response points in MHz where gain values are specified.
%   Suggest to set freqMHzVec(1) to 0.
%   gainDbVec          - Gain response in dB corresponding to the freqMHzVec.
%   Suggest to set gainDbVec(1) to 0.
%   sampleRate      - Sampling rate in Hz.
%   nTaps_isFFT     - Number of filter coefficients or set to 'FFT' to
%   creaste filter by freq. domain.
%   isGainNorm      - Optional. If fure, normalize the gain for max gain is
%   0dB. Default is false.
%   interpMethod    - Optional. Interpolation mehtod. Default is 'linear'.
%   fnum            - Optional.
%
% Outputs:
%   varargout{1}    - signalOut or filter coefficients
%   varargout{2}    - filter coefficients

if 0 % Example
    freqMHzVec = [0 20]
    gainDbVec = [0 -3]
    nTaps_isFFT = 2^10 % applyFilter = 'CONV'
    sampleRate = 430080000;
elseif 1*0
    freqMHzVec = [-20 0 20]
    gainDbVec = [-3 0 -3]
    nTaps_isFFT = 'FFT'
    sampleRate = 430080000;
end

if length(gainDbVec)~=length(freqMHzVec)
    error('check the legnth of gainDbVec and freqMHzVec')
end
isFnum = exist('fnum','var')&&~isempty(fnum);
if ~exist('isGainNorm','var')||isempty(isGainNorm)
    isGainNorm = 0;
end
if ~exist('nTaps_isFFT','var')||isempty(nTaps_isFFT)
    nTaps = 0;
    applyFilter = 'CONV';
elseif isnumeric(nTaps_isFFT)
    nTaps = nTaps_isFFT;
    applyFilter = 'CONV';
elseif strcmpi(nTaps_isFFT, 'FFT')
    nTaps = length(signal);
    applyFilter = 'FFT';
elseif strcmpi(nTaps_isFFT, 'CONV')
    nTaps = 0;
    applyFilter = 'CONV';
else
    error('!')
end
isSignal = exist('signal','var')&&~isempty(signal);
isFreqNegative = any(sign(freqMHzVec)==-1); % check frequency includes negative freq. points
if isFreqNegative * strcmpi(applyFilter,'CONV')
    error('isFreqNegative and CONV cannot coexist for EVM')
end

freqHz = freqMHzVec(:)*1e6;
gainDbVec = gainDbVec(:);

% nyquist frequency
fnz = sampleRate/2;

switch applyFilter
    case {'FFT','Apply filter in freq. domain'}
        error('2024-10-17 bug')
        x = signal;
        nTaps = length(x);

        % create gain curve
        [gain_vec, freqHz_vec] = interpGainCurve(freqHz,gainDbVec,fnz,nTaps,isFreqNegative,isGainNorm);

        if 0
            figure(123), plot(freqHz_vec/1e6, 20*log10(gain_vec))
            plotSpectrum(y, fnz, 1235789)
            plotSpectrum(b, fnz, 1235789)
        end
        
        % create filter coefficients
        if 1
            b = ifft(fftshift(gain_vec));
        else
            b = ifft(fftshift(gain_vec.*exp(1i*linspace(-pi,pi,nTaps)')));
        end

        if 0
            B = abs(fftshift( fft(b,numel(b)) )).^2;
            figure(234), plot(freqHz_vec/1e6, 10*log10(B))
        end

        % apply filter by fft
        if isSignal
            Y = fft(x(:)).*fftshift(gain_vec);
            y = ifft(Y);
            signalOut = y;

            varargout{1} = signalOut;
            varargout{2} = b;
        else
            varargout{1} = b;
        end

    case {'CONV','Apply filter in time domain'}
        % set number of taps ranges
        nTapsRang = [2^(16-1), 2^6];

        % get number of required taps
        if nTaps==0
            df = min(abs(diff(freqMHzVec))) * 1e6;
            nTaps = max([2^nextpow2((fnz) / df), nTapsRang(1)]);
            nTaps = min([nTaps, nTapsRang(end)]);
        else
            nTaps = 2^nextpow2(nTaps);
        end

        % % % update freq step
        % % df = fnz / nTaps;

        % design the arbitrary response magnitude filter using while loop for
        % converging the magnitude er
        errorDb = 1000; % initial errorDb
        toleranceDb = 0.2; % errDb tolearnce
        while (errorDb > toleranceDb && nTaps <= max(nTapsRang))   
            % create gain curve
           [gain_vec, freqHz_vec, taps_vec] = interpGainCurve(freqHz,gainDbVec,fnz,nTaps,isFreqNegative,isGainNorm);

            if 0
                figure(123), plot(taps_vec, gainDb_vec), hold on
                figure(123), plot(taps_vec, gain_vec), hold on
            end

            % Arbitrary response magnitude filter object
            d = fdesign.arbmag('N,F,A',nTaps-1*0,taps_vec,gain_vec);

            % filter design
            if 0
                hd = design(d,'freqsamp','SystemObject',true);
            else
                hd = design(d,'freqsamp','window','hann','SystemObject',true);
            end
            if 0
                fvtool(hd,'magnitude');
                fvtool(hd,'phase');
                figure(123), plot(freqHz_vec, gainDb_vec), hold on
                figure(123), plot(freqHz_hd, gainDb_hd), hold on
            end
            gainDb_hd = 20*log10(hd.measure.Amplitude(:));
            freqHz_hd = hd.measure.Frequencies(:) * fnz;
            
            gainDb_vec = interp1(freqHz_hd, gainDb_hd, freqHz_vec, 'linear');
            % calculate error
            errorDb = max(abs(gainDb_vec - gainDb_hd));

            % update number of taps by x2
            nTaps = nTaps*2;
        end

        % create filter coefficients
        b = hd.Numerator(:);

        % 2024-10-17, Normalize filter coefficients
        b = b/sum(b);

        % apply filter by conv 'same'
        if isSignal
            signalOut = conv(signal(:), b, 'same');

            Nsamps = numel(signal);
            
            varargout{1} = signalOut;
            varargout{2} = b;
        else
            varargout{1} = b;
        end

        if 0 % evl and debug
            [H, w] = freqz(b, 1, 2^14);
            figure(1234), plot(w/pi, angle(H)/pi*180); hold on
            figure(1234), plot(w/pi, abs(w)); hold on
            evm(signal, signalOut)
        end
end

if isFnum % plot
    figure(fnum),
    Nfft = 2^14;
    freqMHz_s = (-Nfft/2:Nfft/2-1)'.*sampleRate/Nfft/1e6;
    if isSignal
        plot(freqMHz_s, 10*log10(1000*abs(fftshift(fft(signal,Nfft)/Nfft).^2)), '-.'),grid on, hold on
        plot(freqMHz_s, 10*log10(1000*abs(fftshift(fft(signalOut,Nfft)/Nfft).^2)), '-.'),grid on, hold on
    else
        plot(freqMHz_s, 10*log10(1000*abs(fftshift(fft(bCoef,Nfft)/Nfft).^2))),grid on, hold on
    end
    xlabel('Frequency [MHz]'), ylabel('Power [dBm]'), title('Spectrum')
end

end

% 2024-07-30, Add interpGainCurve
function [gain_vec, freqHz_vec, taps_vec] = interpGainCurve(freqHz,gainDb,fnz,nTaps,isNagtiveFreq,isGainNorm)
% get magnitude using linear interpolation
if ~isNagtiveFreq
    % taps_vec = (0:nTaps)' / nTaps;
    taps_vec = linspace(0,nTaps,nTaps)/nTaps;
    freqHz_vec = taps_vec' * fnz;
else
    taps_vec = linspace(-1,1,nTaps);
    freqHz_vec = taps_vec' * fnz;
end

if 1
    gainDb_vec = interp1(freqHz, gainDb, freqHz_vec, 'linear', 'extrap');
else
    gainDb_vec = diff(gainDb) / diff(freqHz) * freqHz_vec;
end

scale = 1;
gain_vec = sqrt(1/scale)*10.^(gainDb_vec(:)/20);
gain_vec = gain_vec / (isGainNorm*(max(gain_vec)) + ~isGainNorm); % normalize magnitidue

if 0
    figure(345), plot(freqHz_vec/1e6, 20*log10(gain_vec))
end
end