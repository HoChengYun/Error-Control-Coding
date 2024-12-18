trans_rec_qam
trans_rec_qam_cyclic_code
figure
semilogy(SNR_db, Pe_no_encoding, '-+', SNR_db, Pe_theory, '-x', SNR_db, Pe_ecc, '-o');
title('bit error rate 128QAM'); 
xlabel('Es/N0 (db)'); 
ylabel('BER'); 
legend('no encoding(exp)', 'no encoding(theory)', 'encoding-cyclic code(exp)', 'Location', 'southwest');

