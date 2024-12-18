function bin_data_hat = decoding(r_dec, g_x, use_ecc, M)
    bin_data_hat = [];temp = []; r = [];
    [rowr, colr] = size(r_dec);
    if (use_ecc == 1)
        for i = 1:colr
            temp = meggitt_decoder(r_dec(i));
            temp = abs([mod(deconv(temp, g_x),2)]);
            bin_data_hat = [bin_data_hat, temp];
        end
    else
        r_bin = dec2bin(r_dec);
        r_bin = num2cell(r_bin);
        r_bin = cellfun(@str2double, r_bin);
        for i = 1:colr
            temp = r_bin(i, :);
            temp = abs([mod(deconv(temp, g_x),2)]);
            bin_data_hat = [bin_data_hat, temp];
        end
    end
end

