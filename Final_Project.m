 
 inputfile = 'C:\Users\Norang Tetarwal\Desktop\Input+noise.wav';
 % Input as well as included Noise File
 noiseprintfile = 'C:\Users\Norang Tetarwal\Desktop\whitenoise.wav'; 
 %=inputfile
 outputfile = 'C:\Users\Norang Tetarwal\Desktop\output1.wav';
 %File where our output after analysis is saved
 residuefile = 'C:\Users\Norang Tetarwal\Desktop\residue.wav';
 % File where the remaining residue is stored

 %Algorithms
 filtering_algorithm = 'wiener'; %ssub, wiener, psub, emsr
 % The rules by which we are Supressing the Noise from our file
 
 estimation_algorithm = 'fft'; %fft, welch (mean)
 noise_averaging = 'max'; % max, mean, rms
 
 %STFT
 L = 2*1024; %frame size
 M = L/8; %step size [typ. L/2, L/4, L/8]
 window_function = 'hann'; %rectwin, hamming, hann

 %Noise reduction
 gain = 1.00; %noise reduction gain
 Rmin=0.00;
 AR_model_order = 170; %AR-model order [1 -> (L-1)]

 %Ephiram and Malah Supression Rule
 frame_weight = 0.98; %previous frame weight [0 -> 1]
 Pr_signal_absence = 0.20; %probability of signal absence [0 -> 0.99]

 %Low Pass filter
 lowpass = true; %optional low pass filter
 f3db = 6000; %cutoff frequency [Hz]

 %Save
 save_result = true;
 
 % Load files
 [wav_file, Fs] = audioread(inputfile);
 
 noise_wavform = audioread(noiseprintfile); %extract noise print
 
 win = window_function;
 window_function = eval([window_function, '(L)']); %generate window function

 % Estimate noise print
 
 noise_spectrum = gain * spectrum_estimator(noise_wavform, ...
 window_function, M, estimation_algorithm, noise_averaging);

 % Noise removal
 disp(size(wav_file(:,1)))
 
 denoise_signal = noise_suppressor(wav_file(:,1), noise_spectrum, ...
 window_function, M, filtering_algorithm, frame_weight, Pr_signal_absence);
 
 %denoise_signal(:,2) = noise_suppressor(wav_file(:,2), noise_spectrum, ...
 %window_function, M, filtering_algorithm, frame_weight, Pr_signal_absence);
 
 % Low pass filter signal
 if lowpass
    denoise_signal = lpfilter(denoise_signal, f3db, Fs); %low pass filter
 end

 % Save files, state parameters in filename
 outputfile = [outputfile(1:end-4), '_', filtering_algorithm, ...
     '_', estimation_algorithm, '_', noise_averaging, ...
     '_', win, ...
     '_N', num2str(L), ...
     '_M', num2str(M), ...
     '_g', num2str(gain), ...
     '_Rm', num2str(Rmin), '.wav'];

 if lowpass
    outputfile = [outputfile(1:end-4), ...
    '_f', num2str(f3db), '.wav'];
 end
 
 if max(strcmpi(estimation_algorithm, {'ar', 'ceps'}))
    outputfile = [outputfile(1:end-4), ...
    '_r', num2str(AR_model_order), '.wav'];
 end
 
 if strcmpi(filtering_algorithm, 'emsr')
    outputfile = [outputfile(1:end-4), ...
       '_a', num2str(frame_weight), ...
       '_q', num2str(Pr_signal_absence), '.wav'];
 end
 
 denoise_signal = denoise_signal/max(abs(denoise_signal));
 
 if save_result
    audiowrite(outputfile,denoise_signal,Fs);
 end
 
 
 
 
 % Spectrum estimator of the Given Signals
 function new_spectrum = spectrum_estimator(dn, window_function, M, algorithm,averaging)
 
 L = length(window_function);

 if strcmpi(algorithm, 'welch')
    new_spectrum = pwelch(dn, window_function, (L-M), 'twosided');
    new_spectrum = sum(window_function)^2/norm(window_function)^2 * new_spectrum; %pwelch scale
 else

 N = length(dn);
 dn = [zeros(L,1); dn; zeros((ceil(N/L)*L - N) + L, 1)]; %zero padding
 N = length(dn); %update length

 new_spectrum = zeros(L, 1); % Padding of zeros Initially at all the points

 for n = 1:M:N-(L-M) %n is the signal sample index

    dn_frame = window_function .* dn(n:n+L-1);

    switch algorithm
        case 'fft'
            Dk = fft(dn_frame);

        
        otherwise
            disp('Error: no recognizable choice of algorithm.')
            break;
    end

    switch averaging
        case 'max'

            i = abs(Dk) > abs(new_spectrum);
            new_spectrum(i) = Dk(i);

        case 'mean'
            new_spectrum = (new_spectrum + Dk) / 2;

        case 'rms'
            new_spectrum = sqrt((new_spectrum.^2 + Dk.^2) / 2);

        otherwise
            disp('Error: no recognizable choice of averaging method.')
            break;
    end

 end

 end

 end



% NOISE SUPRESSOR FRO THE FIVEN AUDIO SIGNAL COMPARING WITH ESTMATED
% SPECTRUM OF NOISE
 function xn_hat = noise_suppressor(yn, Dk, window_function, M, rule, frame_weight, Pr_signal_absence)
 
 L = length(window_function);
 N = length(yn);

 yn = [zeros(L,1); yn; zeros((ceil(N/L)*L - N) + L, 1)]; %zero padding
 Npad = length(yn); %update padded length

 xn_hat = zeros(Npad, 1); %preallocate memory
 Gk = zeros(L, 1);
 Xk_hat = zeros(L, 1);
 
 for n = 1:M:Npad-(L-M) %n is the signal sample index
    
    yn_frame = yn(n:n+L-1);
    Yk = fft(window_function .* yn_frame);

    switch rule
        case 'ssub' %Spectral subtraction
            i = (abs(Yk) - abs(Dk)) > 0;
            Gk(i) = (abs(Yk(i)) - abs(Dk(i))) ./ abs(Yk(i));

        case 'wiener' %Wiener filter solution
            i = (abs(Yk) - abs(Dk)) > 0;

            SY = abs(Yk(i)).^2;
            SD = abs(Dk(i)).^2;

            Gk(i) = (SY - SD) ./ SY;

        case 'psub' %Power spectrum subtraction
            i = (abs(Yk) - abs(Dk)) > 0;

            SY = abs(Yk(i)).^2;
            SD = abs(Dk(i)).^2;

            Gk(i) = sqrt( (SY - SD) ./ SY );

        case 'emsr' %Ephraim-Malah suppression rule
            Gk = emsrc(Yk, Dk,Xk_hat, frame_weight, Pr_signal_absence);

        otherwise
            disp('Error: no recognizable suppression rule.')
            break;
    end

    Xk_hat = Gk .* Yk;
    

    xn_frame = ifft(Xk_hat);

    if M > L/2 %gain correction
        xn_frame = xn_frame .* (1./window_function);
    else
        xn_frame = xn_frame .* (M/sum(window_function));
    end

    xn_hat(n:n+L-1) = xn_hat(n:n+L-1) + xn_frame; %overlap add

 end

 xn_hat = xn_hat(L+1:L+N); %zero padding removal
 
 end




% EPHRAIM AND MALAH SUPRESSION RULE:

 function Gk = emsrc(Yk, Dk,Xk_hat, frame_weight, Pr_signal_absence)
 
    L = length(Yk);
    Gk = zeros(L, 1);

    i = (abs(Yk) > 0) & (abs(Dk) > 0); %divide by zero protector

    Rpost = (abs(Yk(i)).^2 ./ abs(Dk(i)).^2) - 1;

    Rprio = (1-frame_weight) * max(Rpost, 0) ... %positives
            + frame_weight * (abs(Xk_hat(i)).^2 ./ abs(Dk(i)).^2);

    Rprio = Rprio ./ (1 - Pr_signal_absence);

    Rpost(Rpost > 700) = 700; %exponential of theta overflow protection

    theta = (1+Rpost) .* (Rprio ./ (1+Rprio));

    I0 = besseli(0, theta/2); %zero order
    I1 = besseli(1, theta/2); %first order

    M_em = exp(-theta/2) .* ((1+theta).*I0 + theta.*I1);

    Gk(i) = (sqrt(pi)/2) ...
            * sqrt( (1./(1+Rpost)) .* (Rprio./(1+Rprio)) ) .* M_em;

    if Pr_signal_absence > 0 %signal presence uncertainty
        mu = (1-Pr_signal_absence) ./ Pr_signal_absence;
        Lambda = mu .* exp(theta) ./ (1+Rprio);
        Gk(i) = Lambda ./ (1+Lambda ) .* Gk(i);
    end

 end


%   LOW PASS BUTTERWORTH FILTER
 function yn = lpfilter(xn, f0, Fs, order)

 if nargin < 4
    order = 4;
 end

 Wn = f0 / (Fs/2);
 [B, A] = butter(order, Wn, 'low'); %butterworth LP coefficients

 yn = filter(B, A, xn); %filter the signal

 end
