 
% SPECTRUM ESTIMATOR OS THE GIVEN NOISE SIGNALS

 function new_spectrum = spectrum_estimator(dn, window_function, M, algorithm,averaging)
 %spectrum_estimator estimates the spectrum of the given Noisy signal using
 %Short time Fourier transform , and various  Algorithms.
 
 % 'window_function' is a window function describing a frame of size 'L'. 'M' is the
 % length of Size between two successive frames. 'algorithm' gives the type of
 % estimation method which we are using like 'fft','AR-Model'
 % 'averaging' gives the averaging process of different frames.
 % 'max' (maximum noise amplitudes)
 % 'mean' (arithmetic mean)
 % 'rms' (root mean square)
 % 'new_spectrum' is the resulting spectrum of size 'L'.

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