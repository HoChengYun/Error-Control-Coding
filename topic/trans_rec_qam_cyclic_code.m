%para
fc = 3*10^(6);
Ts = 10^(-6);
a = 0.001; %wave amplitude
Es = 127*2*a*a/3;
M = 128;
sample_freq = 10^(7);
sample_time = 1/sample_freq;
No_symbols = 50000; %qam num
bit_frame = 4;
bit_num = bit_frame * No_symbols;
t = [0:sample_time:Ts-sample_time];
phi1 = sqrt(2/Ts)*cos(2*pi*fc*t);
phi2 = sqrt(2/Ts)*sin(2*pi*fc*t);
%set snr
SNR_db = -6:2:30;
SNR = 10.^(SNR_db/10); %Es/N0 = 20 (db)
N0_array = Es./SNR;
%bit error rate
Pe_ecc = [];
%ecc set
g_x = [1 1 0 1];
%bit stream
bin_data = randsample([0 1], bit_num, true);
%S/P
bin_data_P = reshape(bin_data, bit_frame, []).';
dec_data = encoding(bin_data_P, g_x);
m = reshape(dec_data, 1, No_symbols);
m_conj = qammod(m,M);
%generate S vector (basis:phi1, phi2)
S = [];
for i = 0: M-1
    si_wave =  a*real(qammod(i,M))*phi1 + ...
        a*imag(qammod(i,M))*phi2;
    siXphi1 = si_wave .* phi1; siXphi2 = si_wave .* phi2;
    int_siXphi1 = sample_time*sum(siXphi1);
    int_siXphi2 = sample_time*sum(siXphi2);
    Si = [int_siXphi1, int_siXphi2];
    S = [S ; Si];
end
%generate R vactor (basis:phi1, phi2)
for N0 = N0_array
    R = [];
    D = zeros(No_symbols,M);
    D_min = [];
    m_hat = [];
    for i = 1: No_symbols
        si_wave =  a*real(m_conj(i))*phi1 + a*imag(m_conj(i))*phi2;
        siXphi1 = si_wave .* phi1; siXphi2 = si_wave .* phi2;
        int_siXphi1 = sample_time*sum(siXphi1);
        int_siXphi2 = sample_time*sum(siXphi2);
        %add noise (after correlator)
        Ri = [int_siXphi1 + normrnd(0, sqrt(N0/2)), int_siXphi2 + normrnd(0,sqrt(N0/2))];
        R = [R;Ri];
    end
    %generate Di
    for i = 1:M
        distance = -abs((R-S(i,:)).^(2));
        D(:, i) = distance(:, 1) + distance(:, 2);
    end
    D = D.';
    %find m_hat
    [D_max, m_hat] = max(D); m_hat = m_hat-1;
    bin_data_hat_ecc = decoding(m_hat, g_x, 1, M);
    %ber generate
    bit_errorNUM_ecc = sum((bin_data - bin_data_hat_ecc) ~= 0);
    BER_ecc = bit_errorNUM_ecc/bit_num;
    Pe_ecc = [Pe_ecc BER_ecc];
end
%calculate error rate theory
Pe_theory = 4*(1-1/sqrt(M)).*qfunc(sqrt(3*Es./((M-1).*N0_array)));
Pe_theory = M/2*Pe_theory/(M-1);

