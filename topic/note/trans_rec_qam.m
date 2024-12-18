%16QAM
%para
clear
fc = 3*10^(6);
Ts = 10^(-6);
a = 0.001; %wave amplitude
Es = 127*2*a*a/3;
M = 128;
sample_freq = 10^(7);
sample_time = 1/sample_freq;
No_symbols = 50000;
bit_frame = log2(M); m = [];
bit_num = bit_frame * No_symbols;
t = [0:sample_time:Ts-sample_time];
phi1 = sqrt(2/Ts)*cos(2*pi*fc*t);
phi2 = sqrt(2/Ts)*sin(2*pi*fc*t); 
%bit data
bin_data = randsample([0 1], bit_num, true);
%generate S/P
bin_data_P = reshape(bin_data, bit_frame, []).';
%bin to dec
for i = 1:No_symbols
    temp = bin_data_P(i,:);
    temp = num2str(temp);
    temp = temp(~isspace(temp));
    temp = bin2dec(temp);
    m = [m, temp];
end
m_conj = qammod(m,M);
%set snr
SNR_db = [-6:2:30];
SNR = 10.^(SNR_db/10); %Es/N0 = 20 (db)
N0_array = Es./SNR;
%bit error rate
Pe = [];
%generate S vector (basis:phi1, phi2)
S = [];
for i = 0: M-1
    si_wave =  a*real(qammod(i,M))*phi1 + a*imag(qammod(i,M))*phi2;
    siXphi1 = si_wave .* phi1; siXphi2 = si_wave .* phi2;
    int_siXphi1 = sample_time*sum(siXphi1);
    int_siXphi2 = sample_time*sum(siXphi2);
    Si = [int_siXphi1, int_siXphi2];
    S = [S ;Si];
end
%generate R vactor (basis:phi1, phi2)
for snr_i = 1:length(N0_array)
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
        Ri = [int_siXphi1 + normrnd(0, sqrt(N0_array(snr_i)/2)), int_siXphi2 + normrnd(0,sqrt(N0_array(snr_i)/2))];
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
    %dec to bin
    r_bin = dec2bin(m_hat);
    r_bin = num2cell(r_bin);
    r_bin = cellfun(@str2double, r_bin);
    r_bin = r_bin.';r_bin = reshape(r_bin, [1, bit_num]);
    %ber generate
    bit_errorNUM = sum((r_bin - bin_data) ~= 0);
    BER = bit_errorNUM/bit_num;
    Pe = [Pe BER];
    %%look constellation===============
    % if (mod(snr_i, 3) == 0)
    %     figure;
    %     plot(R(:,1), R(:,2), 'x', S(:,1), S(:,2), 'x');
    %     title(sprintf('%d qam  SNR:%.2f dB', M, SNR_db(snr_i)))
    %     xlabel('\phi_{1}(t)') 
    %     ylabel('\phi_{2}(t)') 
    %     grid on
    % end
    %%=================================
end

Pe_no_encoding = Pe;
clearvars -except Pe_no_encoding