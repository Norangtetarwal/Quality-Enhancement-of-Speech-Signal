

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
