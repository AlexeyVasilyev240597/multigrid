% func returns full matrix constructed from 5 vectors
function A_glob = getFullMat(obj)
    [~, N] = size(obj.a);
    N = N + 2;
    
    A = zeros(N, N, N);
    B = zeros(N, N, N);
    B(:, :, 1) = eye(N);
    B(:, :, N) = eye(N);
    C = zeros(N, N, N);
    for i = 2:N-1
        A(:, :, i) = diag([0, obj.d(:, i-1)', 0]);
        B(:, :, i) = diag([1, obj.a(:, i-1)', 1]) +...
            diag([0, obj.e(:, i-1)'], 1) + diag([obj.c(:, i-1)', 0], -1);
        C(:, :, i) = diag([0, obj.b(:, i-1)', 0]);    
    end

    A_glob = zeros(N*N);
    A_glob(1:N, 1:N) = B(:, :, 1);
    for i = 2:N-1
        A_glob(((i-1)*N+1):(i*N), ((i-2)*N+1):((i-1)*N)) = A(:, :, i);
        A_glob(((i-1)*N+1):(i*N), ((i-1)*N+1):(i*N)) = B(:, :, i);
        A_glob(((i-1)*N+1):(i*N), (i*N+1):((i+1)*N)) = C(:, :, i);
    end
    A_glob(((N-1)*N+1):(N*N), ((N-1)*N+1):(N*N)) = B(:, :, N);
end