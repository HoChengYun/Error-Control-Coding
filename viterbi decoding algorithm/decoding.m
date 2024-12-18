function [v_hat, u_hat] = decoding(path)
u_hat = []; v_hat = [];
    for i = 1:length(path)-1
        switch path(i)
            case 0
                if path(i+1) == 0
                    v_hat = [v_hat, logical([0 0])];
                    u_hat = [u_hat, logical([0])];
                elseif path(i+1) == 1
                    v_hat = [v_hat, logical([1 1])];
                    u_hat = [u_hat, logical([1])];
                end
            case 1
                if path(i+1) == 2
                    v_hat = [v_hat, logical([1 0])];
                    u_hat = [u_hat, logical([0])];
                elseif path(i+1) == 3
                    v_hat = [v_hat, logical([0 1])];
                    u_hat = [u_hat, logical([1])];
                end
            case 2
                if path(i+1) == 0
                    v_hat = [v_hat, logical([1 1])];
                    u_hat = [u_hat, logical([0])];
                elseif path(i+1) == 1
                    v_hat = [v_hat, logical([0 0])];
                    u_hat = [u_hat, logical([1])];
                end
            case 3
                if path(i+1) == 2
                    v_hat = [v_hat, logical([0 1])];
                    u_hat = [u_hat, logical([0])];
                elseif path(i+1) == 3
                    v_hat = [v_hat, logical([1 0])];
                    u_hat = [u_hat, logical([1])];
                end
        end
end