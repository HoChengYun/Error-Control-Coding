%consider (7,4) cyclic code
r_num = 7;
g_num = 4;
%LSB g_X , MSB r_X
g_X = logical([1 1 0 1]);
r_X = logical([1 1 0 1 1 0 1]);
v_X = r_X;
error_switch = logical([zeros(1, r_num-1), ones(1, r_num+1)]);
error_current = 0;
%MSB error_compensation
error = [];
count = 0;
register = logical([0, 0, 0]);
update_register = logical([0 0 0]);
for r_in = [r_X , logical(zeros(1, r_num))]
    count = count+1;
    update_register(1) = xor(xor(register(3)&g_X(1),r_in), error_current&error_switch(count));
    update_register(2) = xor(register(3)&g_X(2),register(1));
    update_register(3) = xor(register(3)&g_X(3),register(2));
    register = update_register;
    if error_switch(count) == 1
        error_current = register(1) & (~register(2)) & register(3);
        error = [error, error_current];
        fprintf("\n========%d shift==========\n", count-7);
        fprintf("syndrome register:");disp(register);
        fprintf("buffer register:");disp(fliplr(v_X));
        fprintf("correction:%d", error_current)
        v_X(1) = xor(error_current ,v_X(1));
        v_X = circshift(v_X, -1);
    end
end