 
 % EPHRAIM AND MALAH SUPRESSION RULE:

 function Gk = emsrc(Yk, Dk,Xk_hat, frame_weight, Pr_signal_absence)
 %emsrc  rule as defined where 'Yk' is the noisy signal spectrum,
 % 'Dk' the noise spectrum and 'Xk_hat' is the previous frame of the estimated noiseless signal.
 % 'Gk'is then the resulting frequency gain vector. 'frame_weight' is the
 % frame weighting, and 'Pr_signal_absence' the probability of signal absence.
  

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