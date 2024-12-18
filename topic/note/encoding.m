function s_dec = encoding(s_bin, g_x) 
    %g_x = 1 + x + x^3 -> [1 1 0 1 ...]
    [rows, cols] = size(s_bin); s_dec = []; temp=[];
    for i = 1:rows
        temp = mod(conv(s_bin(i, :), g_x), 2);
        temp = double(temp);
        temp = num2str(temp);
        temp = temp(~isspace(temp));
        temp = bin2dec(temp);
        s_dec = [s_dec; temp];
    end
end