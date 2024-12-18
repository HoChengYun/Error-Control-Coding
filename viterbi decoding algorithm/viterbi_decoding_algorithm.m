clear
path=[0]; distance = [0];
v_hat = []; u_hat = [];
r = logical([1 1 0 1 1 0 1 1 0 0 0 1 0 1]);
v = logical([1 1 0 1 1 0 0 1 0 0 1 0 1 1]);
u = logical([1 1 1 0 1 0 0]);
for k = 1: length(r)/2
    [path, distance] = add_path(path, distance, r);
end
[v_hat, u_hat] = decoding(path(1,:));
fprintf("mini distance path:\n");
for i = 1:length(path(1,:))
    fprintf("S%d",path(1,i));
    if i ~= length(path(1,:))
        fprintf(" -> ");
    end
end
fprintf("\n r = ");fprintf(mat2str(double(r)));
fprintf("\n v = ");fprintf(mat2str(double(v)));
fprintf("\n v_hat = ");fprintf(mat2str(double(v_hat)));
fprintf("\n u = ");fprintf(mat2str(double(u)));
fprintf("\n u_hat = ");fprintf(mat2str(double(u_hat)));
