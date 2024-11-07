% 2024-05-30, Add isLimitSmallSignal
% 2024-06-07, Implement multi-stages by set type CELL to parms_struct_cell
% 2024-06-12, Add dlSlot for power calculation
% 2024-06-19, Add AMAM curve generator
% 2024-07-23, Release v1        
%% 2024-08-08, Apply output saturation

function [y, val_tab, val2] = rfSim_pa_stages(x,modelMethod,parms_struct_cell,dlSlots,fnum)
if iscell(parms_struct_cell)
    nStages = numel(parms_struct_cell); % multi-stages
else
    nStages = 1;
    parms_struct_cell = {parms_struct_cell}; % single-stage
end
if ~exist('dlSlots','var')||isempty(dlSlots)
    dlSlots = true(size(x)); % 2024-06-12, Add dlSlot for power calculation
end
if ~exist('fnum','var')||isempty(fnum)
    isFnum = 0;
else
    isFnum = 1;
end

xin = x(:);
val_tmp = zeros(nStages,6+1);
val2_tmp = zeros(length(x),nStages);
if 0 % debug
    powerDbm2(x)
    xin = powerDbm2(x,'set',-18);
    powerDbm2(xin)
end
for k=1:nStages
    pm = parms_struct_cell{k};
    switch upper(string(modelMethod))
        case 'COMM'
            y = comm_pa_model(xin,pm.Method,pm);
        case 'RFSIM'
            y = rfSim_pa_model(xin,pm.Method,pm);
    end
    if 1
        val_tmp(k,1) = powerDbm(xin(dlSlots),'rms');
        val_tmp(k,2) = powerDbm(xin(dlSlots),'peak');
        val_tmp(k,3) = powerDbm(y(dlSlots),'rms');
        val_tmp(k,4) = powerDbm(y(dlSlots),'peak');
        val_tmp(k,5) = evm(xin(dlSlots),y(dlSlots),'dB');
        val_tmp(k,6) = evm(x(dlSlots),y(dlSlots),'dB');
        val_tmp(k,end) = k;
        val2_tmp(:,k) = y;
    end
    if k==2*0 || 0 % debug
        other to check power ? isTDD ?
        plot_comm([xin, y], [], 'ampm', [], {0611, 'ampm'}); hold on
        powerDbm(y,'peak')
        powerDbm(y,'rms')
    end
    if isFnum % plot
        if isa(fnum,'char')&&strcmpi(fnum,'time')
            isNormalize = 0;
            if k==1
                plot_comm(xin/(isNormalize*(norm(xin)-1)+1), [], 'time', [], {0612, ['pa input']}); hold on
            end
            if 1
                plot_comm(y/(isNormalize*(norm(y)-1)+1), [], 'time', [], {0612, ['pa stage',num2str(k),' output']}); hold on
            end
        else
            % plot_comm([xin, y], [], 'ampm', [], {fnum, ['pa stage',num2str(k)]}); hold on
            plot_comm([x, y], [], 'ampm', [], {fnum, ['pa stage',num2str(k)]}); hold on
        end
    end
    xin = y;
end

if 1 % export - val
    val.Stage = val_tmp(:,end);
    val.PinRmsDbm = val_tmp(:,1);
    val.PinReakDbm = val_tmp(:,2);
    val.PoutRmsDbm = val_tmp(:,3);
    val.PoutPeakDbm = val_tmp(:,4);
    val.EVMDb = val_tmp(:,5);
    val.EVMDb_Cascade = val_tmp(:,6);
    val_tab = struct2table(val);
    val2 = val2_tmp;
end

if isrow(x)
    y = y';
end
end

%% comm_pa_model
function y = comm_pa_model(x,memorylessNonlinearityMethod,paParms_struct)
MLNM = memorylessNonlinearityMethod;
pm = paParms_struct;

if strcmpi(MLNM,"Cubic polynomial")
    switch string(pm.TOISpecification)
        case 'IIP3'
            toiParm = pm.IIP3;
        case 'OIP3'
            toiParm = pm.OIP3;
        case 'IP1dB'
            toiParm = pm.IP1dB;
        case 'OP1dB'
            toiParm = pm.OP1dB;
        case 'IPsat'
            toiParm = pm.IPsat;
        case 'OPsat'
            toiParm = pm.OPsat;
    end
    try
        pm.AMPMConversion;
    catch
        % remove ampm impairments
        pm.AMPMConversion = 0;
        pm.PowerLowerLimit = -1000;
        pm.PowerUpperLimit = 1000;
    end

    % generate MemorylessNonlinearity object 
    paModel = comm.MemorylessNonlinearity(...
        'Method',pm.Method, ...
        'LinearGain',pm.LinearGain,...
        'TOISpecification',pm.TOISpecification, ...
        pm.TOISpecification,toiParm,...
        'AMPMConversion',pm.AMPMConversion, ...
        'PowerLowerLimit',pm.PowerLowerLimit, ...
        'PowerUpperLimit',pm.PowerUpperLimit);

elseif strcmpi(MLNM,"Lookup table")
    if isfield(pm,'LookupTable')
        pm.LookupTable;
    elseif all([isfield(pm,'PinDbm') isfield(pm,'PoutDbm') isfield(pm,'PhaseShiftDeg')])
        pm.LookupTable = [pm.LookupTable.PinDbm, pm.LookupTable.PoutDbm, pm.LookupTable.PhaseShiftDeg];
    end
    try
        lut = pm.Table;
    catch
        lut = pm.LookupTable;
    end
    try
        Rohm = pm.ReferenceImpedance;
    catch
        Rohm = 1;
    end

    % generate MemorylessNonlinearity object
    paModel = comm.MemorylessNonlinearity(...
        'Method',pm.Method, ...
        'Table',lut,...
        'ReferenceImpedance',Rohm);
else
    error("rfSim_pa.comm_pa_model.MLNM")
end
y = paModel(x);

if isreal(x)
    % real export
    y = real(y);
end
end

%% rfSim_pa_model
function y = rfSim_pa_model(x,memorylessNonlinearityMethod,paParms_struct,isLimitSmallSignal)
MLNM = memorylessNonlinearityMethod;
pm = paParms_struct;
pm.ReferenceImpedance = 1;
if ~exist('isAMPMLimit','var')||isempty(isLimitSmallSignal)
    isLimitSmallSignal = 1;
end

if strcmpi(MLNM,"Cubic polynomial")
    % Calculate the linear term c1 and the third order term c3
    c1 = 10^(pm.LinearGain/20);
    c3 = 1;  % init value
    switch pm.TOISpecification
        case 'IIP3'
            c3 = -4*c1/(3*10^((pm.IIP3/10)-3)) / pm.ReferenceImpedance;
        case 'OIP3'
            c3 = -4*c1^3/(3*10^((pm.OIP3/10)-3)) / pm.ReferenceImpedance;
        case 'IP1dB'
            c3 = (2*c1*(10^(19/20)-10)) / ...
                (5*3*10^((pm.IP1dB-30)/10)) / pm.ReferenceImpedance;
        case 'OP1dB'
            c3 = (2*c1*(10^(19/20)-10)) / ...
                (5*3*10^((pm.OP1dB-30-pm.LinearGain+1)/10)) / ...
                pm.ReferenceImpedance;
        case 'IPsat'
            c3 = -(4*c1)/(9*10^((pm.IPsat-30)/10)) / pm.ReferenceImpedance;
        case 'OPsat'
            c3 = -(16*c1^3)/(81*10^((pm.OPsat-30)/10)) / ...
                pm.ReferenceImpedance;
        otherwise
            c3 = 0;
    end

    if 0
        abs_x = (abs(x(:)) + sqrt(realmin)) / pm.ReferenceImpedance; % protect against log10(0)
    else
        pow_x = (real(x.*conj(x)) + 1*sqrt(realmin)) / pm.ReferenceImpedance; % protect against log10(0)
        abs_x = sqrt(pow_x);
    end

    if 0
        sprintf("2024-06-07, Removed")
        if c3 == 0
            saturationSq = Inf;
        else
            saturationSq = -4*c1/(9*c3);
        end
        abs_x(abs_x>=saturationSq) = saturationSq;

    elseif c3 == 0
        % 2024-08-08, Apply output saturation
        saturationSq = 10^((pm.OPsat - pm.LinearGain - 30)/20);
        abs_x(abs_x>=saturationSq) = saturationSq;
        pow_x = abs_x.^2;
    end

    if 0
        yAM = c1*abs_x + c3*(3/4)*abs_x.^3;
    elseif 1*0
        yAM = c1*sqrt(pow_x) + c3*(3/4)*pow_x.*sqrt(pow_x);
    else
        yAM = [sqrt(pow_x), pow_x.*sqrt(pow_x)] * [c1;c3*(3/4)];
    end

    pinDbm = 10*log10(abs_x.^2) + 30;
    if isfield(pm,'AMPMConversion')
        pinPowerLowerLimit = pm.PowerLowerLimit - pm.LinearGain;
        pinPowerUpperLimit = pm.PowerUpperLimit - pm.LinearGain;
        pinDbm(pinDbm<pinPowerLowerLimit) = pinPowerLowerLimit;
        pinDbm(pinDbm>=pinPowerUpperLimit) = pinPowerUpperLimit;
        ampmConversion = pm.AMPMConversion / (pinPowerUpperLimit-pinPowerLowerLimit);
        phsPM = deg2rad(ampmConversion * (pinDbm-pinPowerLowerLimit));
        if 0
            randPM = 2 * randi([1, 2], numel(phsPM), 1) - 3;
        else
            randPM = 1;
        end
    else
        randPM = 1;
        phsPM = 0;
    end

    yPM = randPM.*phsPM;

    isMemoryEffect = 0;
    if isMemoryEffect
        if 0
            x = paInput;
            y = paOutput;
        end
        P = 3;
        M = 2;

        p1 = 0;
        pStep = 1; % 
        pEnd = P-1;
        m1 = 0;
        mEnd = M-1;

        count = 0;
        p_m_vec = [];
        X_pm = zeros(length(x), P*M);
        for p=p1:pStep:pEnd
            for m=m1:mEnd
                count = count + 1;
                x_m = zeros(size(x));
                x_m(m+1:end) = x(1:end-m);
                x_p = sqrt(real(x_m.*conj(x_m))).^p;
                X_pm(:,count) = x_m.*x_p;
                p_m_vec = [p_m_vec; [count p m]];
            end
        end

        if 0 % pa model coefficients
            coef = [c1 0 c3*(3/4) cO1 c11 c21 cO2 c12 c22]';
        else
            coef = zeros((P*M),1);
            coef([1,2,3],:) = [c1 0 c3*(3/4)]';
        end

        y = X_pm * coef;

        % A least-squares solution
        if 0
            coef = X_pm(memLen:xLen,:)\y(memLen:xLen);
        else
            coef = X_pm(1:xLen,:)\y(1:xLen);
        end

        if 0
            plot_comm([x])
        end

        % Output coefficient matrix (coefMat)
        varargout{1} = reshape(coef,memLen,numel(coef)/memLen);
    end

    y = yAM.*exp(1i.*(angle(x) + yPM)) .* 1;

    if 0
        PLOT_COMM.sa(PLOT_COMM, y , 300e6, [], [], {123}, [], 0.1 )

        y_2 = [x, x.^3] * [c1;c3*(3/4)];
        PLOT_COMM.sa(PLOT_COMM, y_2 , 300e6, [], [], {123}, [], 0.1 )

        PLOT_COMM.sa(PLOT_COMM, y , 300e6, [], [], {123}, [], 0.1 )
        y_tmp = c1*x + c3*x.^3;
        PLOT_COMM.sa(PLOT_COMM, y_tmp , 300e6, [], [], {123} )
        y_hd2 = (c1*sqrt(pow_x) + 0*(3/4)*pow_x.*sqrt(pow_x) + 10*abs_x.^2) .* exp(1i.*(angle(x) + yPM));
        PLOT_COMM.sa(PLOT_COMM, y_hd2 , 300e6, [], [], {234} )

        x = 1*x;

        y_hd2_tmp = 5*x + 10*x.^2 + 0*3/4*x.^3;

        c2 = sqrt(c1)

        y_tmp1 = (c1*abs(x) + c2*abs(x).^2 + c3*abs(x).^3) .*  exp(1i.*(angle(x)));
        y_tmp2 = (c1*x + c2*x.^2 + c3*x.^3);

        PLOT_COMM.sa(PLOT_COMM, y_tmp1 , 300e6, [], [], {123} )
        PLOT_COMM.sa(PLOT_COMM, y_tmp2 , 300e6, [], [], {123} )

        y_tmp3 = c1*x + c2*abs(x).*x.* exp(1i.*(angle(x)));
        PLOT_COMM.sa(PLOT_COMM, y_tmp3 , 300e6, [], [], {123} )

        y_tmp4 = c1*x + c2*x.*x + c3*x.*abs(x).^2 + c3*abs(x).^3.*exp(1i.*(angle(x)));
        PLOT_COMM.sa(PLOT_COMM, y_tmp4 , 300e6, [], [], {123} )

        powerDbm(y_tmp4)
        powerDbm(y_tmp1)
    end

elseif strcmpi(MLNM,"Lookup table")
    if iscell(pm.LookupTable)
        isMultiLUT = 1;
    else
        isMultiLUT = 0;
    end

    abs_x = (abs(x(:)) + sqrt(realmin)) / pm.ReferenceImpedance; % protect against log10(0)
    pi_dBw_sig = 10*log10(abs_x.^2);

    % AM/AM, AM/PM interpolation
    sprintf('unit: dBw and rad');
    if isMultiLUT
        lut_amam = pm.LookupTable{1};
        lut_ampm = pm.LookupTable{2};
        pi_dBm_lut = lut_amam.PinDbm;
        po_dBm_lut = lut_amam.PoutDbm;

        pi_dBm_lut_ampm = lut_ampm.PinDbm;
        ps_deg_lut = lut_ampm.PhaseShiftDeg;

        lut_matrix_sig = zeros(length(pi_dBw_sig),2);
        lut_matrix_sig(:,1) = interp1([pi_dBm_lut-30],[po_dBm_lut-30],pi_dBw_sig,'linear','extrap');
        lut_matrix_sig(:,2) = interp1([pi_dBm_lut_ampm-30],ps_deg_lut,pi_dBw_sig,'linear','extrap');
    else
        if istable(pm.LookupTable)
            pi_dBm_lut = pm.LookupTable.PinDbm;
            po_dBm_lut = pm.LookupTable.PoutDbm;
            ps_deg_lut = pm.LookupTable.PhaseShiftDeg;
        elseif ismatrix(pm.LookupTable)
            pi_dBm_lut = pm.LookupTable(:,1);
            po_dBm_lut = pm.LookupTable(:,2);
            ps_deg_lut = pm.LookupTable(:,3);
        else
            error('pm.LookupTable !?')
        end
        lut_matrix = [[pi_dBm_lut, po_dBm_lut]-30, ps_deg_lut*pi/180];
        lut_matrix_sig = interp1(lut_matrix(:,1),lut_matrix(:,2:3),pi_dBw_sig,'linear','extrap');
    end
    voltGain = 10.^((lut_matrix_sig(:,1) - pi_dBw_sig)/20);
    phsShift_rad = lut_matrix_sig(:,2);

    gn_dB_lut = po_dBm_lut - pi_dBm_lut;
    gnln_dB_lut = gn_dB_lut(1); % linear gain
    if isLimitSmallSignal % 2024-05-30, Add isLimitSmallSignal
        pin_dBm_lut_min = min(pi_dBm_lut);
        idx_pin_smsg = find((pi_dBw_sig+30) <= pin_dBm_lut_min); % index of small signal
        phsShift_rad(idx_pin_smsg) = 0; % set phase shift 0 at small signal points
        voltGain(idx_pin_smsg) = 10.^(gnln_dB_lut/20); % set linear voltGain at small signal points
    end

    % combine and scale
    if 0
        randPM = 2 * randi([1, 2], numel(phsShift_rad), 1) - 3;
    else
        randPM = 1;
    end
    y = voltGain .* x .* exp(1i*randPM.*phsShift_rad);

    if 0 % debug
        plot_comm([x, y], [], 'ampm', [], {0530, 'ampm'}); hold on
        figure(0530), subplot(1,3,1), plot(pi_dBm_lut,po_dBm_lut,LineWidth=2,DisplayName='LUT')
        figure(0530), subplot(1,3,2), plot(po_dBm_lut,gn_dB_lut,LineWidth=2,DisplayName='LUT')
        figure(0530), subplot(1,3,3), plot(po_dBm_lut,ps_deg_lut,LineWidth=2,DisplayName='LUT')
    end

    if 0
        pin_dBm = pi_dBw_sig + 30;
        po_dBm = 20*log10(abs(y)) + 30;
        gc_dB = po_dBm - pin_dBm;
        phsShift_deg = phsShift_rad/pi*180;
        figure(0530), subplot(1,3,1), plot(pin_dBm,po_dBm), hold on
        figure(0530), subplot(1,3,1), plot(pi_dBm_lut, po_dBm_lut, 'LineWidth', 2), hold on
        figure(0530), subplot(1,3,2), plot(pi_dBm_lut, gn_dB_lut, 'LineWidth', 2), hold on
        figure(0530), subplot(1,3,2), plot(pin_dBm, gc_dB, 'LineWidth', 2), hold on
        figure(0530), subplot(1,3,3), plot(pin_dBm, phsShift_deg), hold on
        figure(0530), subplot(1,3,3), plot(pi_dBm_lut, ps_deg_lut, 'LineWidth', 2), hold on
    end

    if 0
        plot_comm([paInSignal, paSignal], fs, 'ampm', pltParms, {0530, ["evm:"+string(max(round(output.evmRMS,2)))+"%"]}); hold on
    end

end
end