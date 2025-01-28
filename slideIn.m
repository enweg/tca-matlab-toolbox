function B = slideIn(B, A)
    % Ensure A can slide into B
    if mod(size(B, 1), size(A, 1)) ~= 0
        error('A cannot slide into B because the number of rows of B is not an integer multiple of the number of rows of A.');
    end

    % Ensure B is square
    if size(B, 1) ~= size(B, 2)
        error('B must be square.');
    end

    % Number of horizontal blocks
    nHorizontalBlocks = floor(size(B, 1) / size(A, 1));
    K = size(A, 1);

    for i = 1:nHorizontalBlocks
        % Extract the block from A
        block = A(:, max(1, end - i * K + 1):end);
        
        % Define row and column indices for insertion
        r = ((i - 1) * K + 1):(i * K);
        c = (i * K - size(block, 2) + 1):(i * K);
        
        % Insert block into B
        B(r, c) = block;
    end
end
