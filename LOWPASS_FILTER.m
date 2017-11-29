 %   LOW PASS BUTTERWORTH FILTER
 function yn = lpfilter(xn, f0, Fs, order)

 % 'xn' is the input signal, 'f0' the cut off frequency, 'Fs' the
 % sampling time and the filter 'order'. Then 'yn' is the filtered output
 % signal.

 if nargin < 4
    order = 4;
 end

 Wn = f0 / (Fs/2);
 [B, A] = butter(order, Wn, 'low'); %butterworth LP coefficients

 yn = filter(B, A, xn); %filter the signal

 end