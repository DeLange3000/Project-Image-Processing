% Example input matrices
A = [1 2 2; 3 3 1; 1 1 2];
B = [10 20 30; 40 50 60; 70 80 90];

% Initialize an empty array to store the mean values
mean_values = [];

% Iterate through each unique value in the first matrix
for value = unique(A)'
    
    % Find the indices where the value occurs in the first matrix
    indices = find(A == value);
    
    % Extract the values from the second matrix that correspond to the indices
    values = B(indices);
    
    % Compute the mean value of the extracted values
    mean_value = mean(values);
    
    % Append the mean value to the array
    mean_values = [mean_values mean_value];
end

% Display the mean values
disp(mean_values);


mean([1 1 1 1 1 1 1; 2 2 2 2 2 2 2; 3 3 3 3 3 3 3 ])