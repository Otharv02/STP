%EXPERIMENT 1
% Exercise No:1
% Solution
clear
clc

a = 10;
x = 10;
y = 20;

z = exp(-a) * sin(x) + 10 * sqrt(y);

disp('Exp 1')
disp(z)
% Exercise No:2
% Solution


clc; clear; close all;
% Experiment No.: 5

% Exercise.: 1
% 1) Plot the function y = 3x^3 – 26x + 10 , and its first and second derivatives, for –2 ≤ x 4 , all in the same plot.
% Define x range
x = -2:0.1:4;

% Define the function
y = 3*x.^3 - 26*x + 10;

% First derivative: y' = 9x^2 - 26
y1 = 9*x.^2 - 26;

% Second derivative: y'' = 18x
y2 = 18*x;

% Plot all functions
figure;
plot(x, y, 'b', 'LineWidth', 2); hold on;
plot(x, y1, 'r--', 'LineWidth', 2);
plot(x, y2, 'g:', 'LineWidth', 2);

% Add grid and labels
grid on;
xlabel('x');
ylabel('y, y'', y''''');
title('Function y = 3x^3 - 26x + 10 and its derivatives');
legend('y = 3x^3 - 26x + 10', 'y'' = 9x^2 - 26', 'y'''' = 18x');




% Exercise.: 2
% Define the range of x
x = -3:0.1:5;

% Define the function f(x)
f = ((x + 5).^2) ./ (4 + 3*x.^2);

% Plot the function
figure;
plot(x, f, 'b', 'LineWidth', 2);

% Add labels and title
xlabel('x');
ylabel('f(x)');
title('Plot of f(x) = (x + 5)^2 / (4 + 3x^2)');
grid on;

% Exercise.: 3
% Wheatstone Bridge Voltage Difference Plots
% Given: v = 12 V, R3 = R4 = 250 Ω

% Given constants
v = 12;
R3 = 250;
R4 = 250;

R2 = 120;
R1 = linspace(0, 500, 500);

vAB_R1 = v * ((R2 ./ (R1 + R2)) - (R4 / (R3 + R4)));

R1_fixed = 120;
R2 = linspace(0, 500, 500);

vAB_R2 = v * ((R2 ./ (R1_fixed + R2)) - (R4 / (R3 + R4)));

figure;

subplot(2,1,1);
plot(R1, vAB_R1, 'b', 'LineWidth', 1.5);
grid on;
title('v_{AB} vs R_1  (R_2 = 120 Ω)');
xlabel('R_1 (Ω)');
ylabel('v_{AB} (V)');

subplot(2,1,2);
plot(R2, vAB_R2, 'r', 'LineWidth', 1.5);
grid on;
title('v_{AB} vs R_2  (R_1 = 120 Ω)');
xlabel('R_2 (Ω)');
ylabel('v_{AB} (V)');

sgtitle('Wheatstone Bridge Voltage Difference Characteristics');


% Exercise.: 4
% Data
x = 0:0.5:100;
y = cos(x);

% Plot
figure;
plot(x, y, 'r--s', ...
    'MarkerEdgeColor', 'm', ...    % magenta edges
    'MarkerFaceColor', 'c', ...    % cyan face
    'MarkerSize', 5, ...           % marker size 5
    'LineWidth', 1.2);             % optional for better visibility

% Labels and title
title('Cosine Signal');
xlabel('Time (s)');
ylabel('Amplitude (mV)');

% Axes background color
set(gca, 'Color', 'y');  % yellow background

% Grid for clarity
grid on;


% Exercise.: 5
% Audio waveform generation with specified tones and spectrogram

% Parameters
Fs = 8000; % Sampling frequency (8 kHz)
tone_duration = 0.5; % seconds per tone
t = 0:1/Fs:tone_duration-1/Fs;

% Frequencies for the 9 tones (in Hz)
freqs = [659 622 659 622 659 494 587 523 440];

% Generate the tone sequence
y = [];
for f = freqs
    tone = sin(2*pi*f*t);
    y = [y tone];
end

% Normalize amplitude
y = y / max(abs(y));

% Create a copy with every other sample multiplied by -1
y2 = y;
y2(2:2:end) = -y2(2:2:end);

% Append the modified copy to the original
final_waveform = [y y2];

% Play the waveform
sound(final_waveform, Fs);

% Write to a WAV file
audiowrite('nine_tones_sequence.wav', final_waveform, Fs);

% Plot spectrogram
figure;
spectrogram(final_waveform, 256, 128, 256, Fs, 'yaxis');
title('Spectrogram of Generated Nine-Tone Sequence');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

clc; clear; close all;

% Experiment No.: 8

% Exercise 2
% Daily maximum temperatures in °F for April 2002 (Washington, DC)
T = [58 73 73 53 50 48 56 73 73 66 69 63 74 82 84 91 93 89 91 80 ...
     59 69 56 64 63 66 64 74 63 69];

days = 1:length(T); % Days of the month (1–30)

%% (a) Number of days temperature was above 75°F
above75 = T > 75;
num_above75 = sum(above75);

fprintf('EXERCISE: Q2\n');
fprintf('(a) Number of days above 75°F: %d\n', num_above75);

%% (b) Number of days temperature was between 65°F and 80°F
between65_80 = (T >= 65) & (T <= 80);
num_between65_80 = sum(between65_80);

fprintf('(b) Number of days between 65°F and 80°F: %d\n', num_between65_80);

%% (c) Days of the month when temperature was between 50°F and 60°F
between50_60 = (T >= 50) & (T <= 60);
days_between50_60 = days(between50_60);

fprintf('(c) Days between 50°F and 60°F: ');
disp(days_between50_60);


% Exercise 3
fprintf('EXERCISE: Q3\n');
% Input from user
hours = input('Enter the number of hours worked: ');
wage  = input('Enter the hourly wage: ');

% Compute pay
if hours <= 40
    pay = hours * wage;
else
    overtime_hours = hours - 40;
    pay = (40 * wage) + (overtime_hours * wage * 1.5);
end

% Display result
fprintf('Total pay = $%.2f\n', pay);

fprintf('\n');


% Exercise 4
fprintf('EXERCISE: Q4\n')
% Input from user
n = input('Enter number of rows (n): ');
m = input('Enter number of columns (m): ');

% Initialize matrix
A = zeros(n, m);

% Fill the first row (column numbers)
A(1, :) = 1:m;

% Fill the first column (row numbers)
A(:, 1) = (1:n)';

% Fill remaining elements
for i = 2:n
    for j = 2:m
        A(i, j) = A(i-1, j) + A(i, j-1);
    end
end

% Display the matrix
disp('Generated Matrix:');
disp(A);


fprintf('\n');

% Exercise 5
fprintf('EXERCISE: Q5\n')
% Parameters
P0 = 300000;        % initial principal $
r = 0.05;           % annual interest rate (5%)
infl = 0.02;        % annual inflation rate on withdrawals (2%)
w0 = 25000;         % first-year withdrawal (end of year 1)
maxYears = 200;     % safety cap

% Initialize
balance = P0;
year = 0;
withdrawals = [];
balances_after = [];  % balance after each withdrawal

while balance > 0 && year < maxYears
    year = year + 1;
    w = w0 * (1 + infl)^(year-1);       % withdrawal at end of year 'year'
    balance = balance * (1 + r) - w;    % interest then withdrawal at year end
    withdrawals(end+1) = w;
    balances_after(end+1) = balance;
end

% Determine number of full years the planned withdrawal was fully paid
if isempty(balances_after)
    full_years = 0;
else
    % If last balance is negative, the last full year was year-1
    if balances_after(end) < 0
        full_years = length(balances_after) - 1;
    else
        full_years = length(balances_after);
    end
end

% Display results
fprintf('Initial principal: $%0.2f\n', P0);
fprintf('Annual interest: %0.2f%%, Withdrawal inflation: %0.2f%%\n', r*100, infl*100);
fprintf('Initial withdrawal (end of year 1): $%0.2f\n\n', w0);
fprintf('Number of FULL years the account can fully cover the planned withdrawal: %d years\n', full_years);

if balances_after(end) < 0
    fprintf('In year %d the planned withdrawal of $%0.2f cannot be fully paid; remaining balance after previous year was $%0.2f\n', ...
        full_years+1, withdrawals(full_years+1), balances_after(full_years));
else
    fprintf('Account remains positive after %d years (within simulation cap of %d years).\n', full_years, maxYears);
end

% Prepare vectors for plotting (balance just after each withdrawal)
years_vec = 1:length(withdrawals);

figure('Color','w','Position',[200 200 800 600]);

subplot(2,1,1);
plot(years_vec, balances_after, '-o', 'LineWidth', 1.6);
hold on;
yline(0, '--k', 'LineWidth', 1); % zero line
if full_years>0
    xline(full_years, '--r', sprintf('Last full year = %d', full_years), 'LabelHorizontalAlignment','left');
end
grid on;
xlabel('Year');
ylabel('Account balance ($)');
title('Account balance after each year''s withdrawal');

subplot(2,1,2);
bar(years_vec, withdrawals);
hold on;
if full_years>0
    % highlight the first year that cannot be fully paid (if exists)
    if length(withdrawals) > full_years
        bar(full_years+1, withdrawals(full_years+1), 'FaceColor',[0.85 0.33 0.10]); % different color
        legend('Planned withdrawals','First unpaid planned withdrawal','Location','NorthWest');
    else
        legend('Planned withdrawals','Location','NorthWest');
    end
else
    legend('Planned withdrawals','Location','NorthWest');
end
grid on;
xlabel('Year');
ylabel('Withdrawal ($)');
title('Yearly planned withdrawals (inflation-adjusted)');

% Improve layout
sgtitle('Retirement withdrawals and account balance over years');


% Exercise 6
fprintf('EXERCISE: Q6\n')
% Define the matrix size
rows = 5;   % You can change this
cols = 5;   % You can change this

% Initialize the matrix with zeros
A = zeros(rows, cols);

% Loop through each element
for i = 1:rows
    for j = 1:cols
        A(i, j) = (i + j) / (j^2);
    end
end

% Display the result
disp(A);


% Exercise 7
fprintf('EXERCISE: Q7\n')
n = 1;  % start from 1

while true
    if mod(n, 2) == 1 && mod(n, 11) == 0 && sqrt(n) > 132
        fprintf('The required number is: %d\n', n);
        break;
    end
    n = n + 1;
end
fprintf('\n')

% Exercise 8
fprintf('EXERCISE: Q8\n')

x = [-3.5 5 -6.2 11.1 0 7 -9.5 2 15 -1 3 2.5];  % Given vector
n = length(x);

% Sorting using loops and conditional statements (Selection Sort)
for i = 1:n-1
    % Assume the current element is the minimum
    min_index = i;
    
    % Check the rest of the vector for a smaller element
    for j = i+1:n
        if x(j) < x(min_index)
            min_index = j;
        end
    end
    
    % Swap elements if a smaller one was found
    if min_index ~= i
        temp = x(i);
        x(i) = x(min_index);
        x(min_index) = temp;
    end
end

% Display the result
disp('The sorted vector is:')
disp(x)
fprintf('\n')

% Exercise 9
fprintf('EXERCISE: Q9\n')

% Script: average_top8_scores.m
% This program calculates the average of the top 8 exam scores
% without using MATLAB's built-in sort function.

scores = [73 91 37 81 63 66 50 90 75 43 88 80 79 69 26 82 89 99 71 59];
n = length(scores);

% Sort scores in descending order using loops (Selection Sort)
for i = 1:n-1
    max_index = i;
    for j = i+1:n
        if scores(j) > scores(max_index)
            max_index = j;
        end
    end
    % Swap if a larger value is found
    if max_index ~= i
        temp = scores(i);
        scores(i) = scores(max_index);
        scores(max_index) = temp;
    end
end

% Take the top 8 scores
top8 = scores(1:8);

% Calculate the average
avg_top8 = sum(top8) / length(top8);

% Display the results
disp('The top 8 scores are:')
disp(top8)
fprintf('The average of the top 8 scores is: %.2f\n', avg_top8);
fprintf('\n')

% Exercise 10
fprintf('EXERCISE: Q10\n')
% 13,15,4,18,19,3
% Ask the user to enter the coordinates of the three points
x1 = input('Enter x-coordinate of first point: ');
y1 = input('Enter y-coordinate of first point: ');
x2 = input('Enter x-coordinate of second point: ');
y2 = input('Enter y-coordinate of second point: ');
x3 = input('Enter x-coordinate of third point: ');
y3 = input('Enter y-coordinate of third point: ');

% Compute the terms needed for solving the circle equation
A = 2 * (x2 - x1);
B = 2 * (y2 - y1);
C = x2^2 + y2^2 - x1^2 - y1^2;

D = 2 * (x3 - x1);
E = 2 * (y3 - y1);
F = x3^2 + y3^2 - x1^2 - y1^2;

% Solve for circle center (h, k)
denominator = (A * E - B * D);
if denominator == 0
    error('The three points are collinear — no unique circle exists.');
end

h = (C * E - B * F) / denominator;
k = (A * F - C * D) / denominator;

% Compute the radius
r = sqrt((x1 - h)^2 + (y1 - k)^2);

% Display the results
fprintf('\nThe center of the circle is (%.3f, %.3f)\n', h, k);
fprintf('The radius of the circle is %.3f\n', r);

% Plot the circle and points
theta = linspace(0, 2*pi, 300);
xc = h + r * cos(theta);
yc = k + r * sin(theta);

figure;
plot(xc, yc, 'b-', 'LineWidth', 1.5); hold on;
plot([x1 x2 x3], [y1 y2 y3], 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(h, k, 'ko', 'MarkerFaceColor', 'g'); % center point
axis equal;
grid on;
title('Circle Passing Through Three Points');
xlabel('x-axis');
ylabel('y-axis');
legend('Circle','Given Points','Center','Location','best');
fprintf('\n')

% Exercise 11
fprintf('EXERCISE: Q11\n')

n = 1; % start from 1
while true
    % Compute the sum of integers from 1 to n
    s = n * (n + 1) / 2;
    
    % Check if the sum is within the desired range
    if s >= 100 && s <= 1000
        % Convert to string to check if all digits are identical
        s_str = num2str(s);
        if length(s_str) == 3 && all(s_str == s_str(1))
            fprintf('The required integer n is: %d\n', n);
            fprintf('The corresponding sum is: %d\n', s);
            break;
        end
    end
    
    % Stop the loop if the sum exceeds 1000 (no need to go further)
    if s > 1000
        disp('No such integer found.');
        break;
    end
    
    n = n + 1; % increment n
end
fprintf('\n')

% Exercise 12
fprintf('EXERCISE: Q12\n')
% Initialize an empty matrix to store twin primes
twinPrimes = [];

% Function to check if a number is prime
isPrimeFunc = @(n) all(mod(n, 2:n-1) ~= 0) && n > 1;

% Loop through numbers from 10 to 498 (last candidate for twin prime pair)
for n = 10:498
    if isPrimeFunc(n) && isPrimeFunc(n + 2)
        % Append the twin prime pair to the matrix
        twinPrimes = [twinPrimes; n, n+2];
    end
end

% Display the twin primes
disp('Twin primes between 10 and 500 are:');
disp(twinPrimes);

%EXPERIMENT 2
% Exercise No:1
fprintf('EXERCISE Q1\n')
a = 10;
b = 20;

% Overwrite variable a
a = 15;

% Display the results
disp(['The value of a is: ', num2str(a)]);
disp(['The value of b is: ', num2str(b)]);
fprintf('\n')

% Exercise No:2
fprintf('EXERCISE Q2\n')
a = 10;
b = 25;

%c = a +   
fprintf('c = a + \n')
fprintf(['Error: File: untitled2.m Line: 23 Column: 64\n' ...
    'Invalid expression. Check for missing or extra characters.\n'])
c = a + b;  
disp(['The value of c is: ', num2str(c)]);
fprintf('\n')

% Exercise No:3
fprintf('EXERCISE Q3\n')
x = 2 + 3 * 4;
y = (2 + 3) * 4;

disp(['x = ', num2str(x)]);
disp(['y = ', num2str(y)]);
fprintf('\n')

% Exercise No:4
fprintf('EXERCISE Q4\n')

piVal = pi;

format short
disp('format short:');
disp(piVal)

format long
disp('format long:');
disp(piVal)

format bank
disp('format bank:');
disp(piVal)
fprintf('\n')

% Exercise No:5
fprintf('EXERCISE Q5\n')

p = 100;
q = 200;

% List variables
disp('Using who:');
who
disp('Using whos:');
whos

% Clear variable p
clear p
disp('After clearing p:');
who

% Clear all variables
%clear all
disp('After clearing all:');
who
fprintf('\n')

% Exercise No:6
fprintf('EXERCISE Q6\n')

a = 10;
b = 20;
c = a + b;
d = a * b;
disp(c);
disp(d);
fprintf('\n')

% Exercise No:7
fprintf('EXERCISE Q7\n')
a = 5; b = 10; c = a + b; d = a * b; disp(c); disp(d);

%EXPERIMENT 3
% Exercise No:1
fprintf('EXERCISE Q1\n')
% Define matrices
M = [1, 3, -1, 6; 2, 4, 0, -1; 0, -2, 3, -1; -1, 2, -5, 1];
N = [-1, -3, 3; 2, 1, 6; 1, 4, -1; 2, -1, 2];

% Multiply M * N
result1 = M * N;
disp('M * N =');
disp(result1);

% Multiply N * M (check if possible)
[mN, nN] = size(N);
[mM, nM] = size(M);

if nN == mM
    result2 = N * M;
    disp('N * M =');
    disp(result2);
else
    disp('N * M cannot be computed because the inner dimensions do not match.');
end
fprintf('\n')

% Exercise No:2
fprintf('EXERCISE Q2\n')

% Define matrices
M = [1, -2, 8, 0];
N = [1 5 6 8; 2 5 6 9];

% Attempt to add M and N
try
    result = M + N;
    disp('M + N =');
    disp(result);
catch ME
    disp('Error:');
    disp(ME.message);
end
fprintf('\n')

% Define matrices
A = [1 6 9 8 5; 9 3 5 8 4; 5 6 3 5 7];
B = [6 5 9 3 5; 6 5 4 8 5; 6 3 5 7 9];
C = [2 5 9 3 4; 5 6 3 7 8; 9 8 6 5 4];

% Compute (A + B + C) transpose
result = (A + B + C)';

% Display the result
disp('[(A + B) + C]^T =');
disp(result);

%EXPERIMENT 4
%  SCALAR FUNCTIONS
disp('SCALAR FUNCTIONS');

A = [-pi/4, pi/2, -3; 2, -5, 4];  % Example matrix

disp('Original Matrix A:');
disp(A);

disp('sin(A):'); disp(sin(A));
disp('cos(A):'); disp(cos(A));
disp('tan(A):'); disp(tan(A));
disp('asin(sin(A)):'); disp(asin(sin(A)));
disp('acos(cos(A)):'); disp(acos(cos(A)));
disp('atan(tan(A)):'); disp(atan(tan(A)));
disp('exp(A):'); disp(exp(A));
disp('log(abs(A)):'); disp(log(abs(A)));   % log requires positive numbers
disp('abs(A):'); disp(abs(A));
disp('sqrt(abs(A)):'); disp(sqrt(abs(A))); % sqrt requires non-negative numbers
disp('rem(A,3):'); disp(rem(A,3));
disp('round(A):'); disp(round(A));
disp('floor(A):'); disp(floor(A));
disp('ceil(A):'); disp(ceil(A));

pause; % Press Enter to continue to next section

%  VECTOR FUNCTIONS
disp('VECTOR FUNCTIONS');
v = [3, -6, 9, 2, 7, -4];

disp('Vector v:');
disp(v);

[max_val, max_idx] = max(v);
[min_val, min_idx] = min(v);

disp(['Maximum value: ', num2str(max_val), ' (at index ', num2str(max_idx), ')']);
disp(['Minimum value: ', num2str(min_val), ' (at index ', num2str(min_idx), ')']);
disp(['Length of v: ', num2str(length(v))]);
disp('Sorted v:'); disp(sort(v));
disp(['Sum of elements: ', num2str(sum(v))]);
disp(['Product of elements: ', num2str(prod(v))]);
disp(['Median value: ', num2str(median(v))]);
disp(['Mean value: ', num2str(mean(v))]);
disp(['Standard deviation: ', num2str(std(v))]);

pause; % Press Enter to continue to next section

%  MATRIX BUILDING FUNCTIONS
disp('MATRIX BUILDING FUNCTIONS');

disp('Identity matrix (eye):');
disp(eye(3));

disp('Matrix of zeros (zeros):');
disp(zeros(3,4));

disp('Matrix of ones (ones):');
disp(ones(2,3));

disp('Create diagonal matrix:');
d = [1 5 9];
disp(diag(d));

M = [4 7 2; 3 8 5; 6 1 9];
disp('Original Matrix M:');
disp(M);

disp('Upper triangular (triu):');
disp(triu(M));

disp('Lower triangular (tril):');
disp(tril(M));

disp('Random matrix (rand):');
disp(rand(3,3));

pause; % Press Enter to continue to next section

%  MATRIX OPERATIONS
disp('MATRIX OPERATIONS');

B = [2 3 1; 4 1 5; 7 2 6];

disp('Matrix B:');
disp(B);

disp(['Size of B: ', num2str(size(B,1)), ' x ', num2str(size(B,2))]);
disp(['Determinant of B: ', num2str(det(B))]);
disp('Inverse of B (inv):');
disp(inv(B));
disp(['Rank of B: ', num2str(rank(B))]);
disp('Reduced Row Echelon Form (rref):');
disp(rref(B));
[V, D] = eig(B);
disp('Eigenvalues and Eigenvectors (eig):');
disp('Eigenvalues (D):'); disp(D);
disp('Eigenvectors (V):'); disp(V);
disp('Characteristic polynomial coefficients (poly):');
disp(poly(B));
fprintf('\n')

% Exercise No:1
fprintf('EXERCISE Q1\n')
A = [4 -2 6;
     2 8 2;
     6 10 3];

% Constant matrix B
B = [8;
     4;
     0];
X = A \ B;
disp('Solution (x, y, z):');
disp(X);
fprintf('\n')

% Exercise No:2
fprintf('EXERCISE Q2\n')
% Given vectors
u = [-3 8 -2];
v = [6.5 -5 -4];

%% (a) Using element-by-element multiplication and sum
dot_a = sum(u .* v);
disp('(a) Dot product using element-by-element and sum:');
disp(dot_a);

%% (b) Using row and column vectors with matrix multiplication
u_row = [-3 8 -2];   % row vector
v_col = [6.5; -5; -4];  % column vector
dot_b = u_row * v_col;
disp('(b) Dot product using matrix multiplication:');
disp(dot_b);

%% (c) Using built-in MATLAB function dot
dot_c = dot(u, v);
disp('(c) Dot product using built-in function:');
disp(dot_c);
fprintf('\n')

%EXPERIMENT 6
% Exercise No:1
fprintf('EXERCISE Q1\n')
% Define range of x
x = -10:0.1:10;

% Define signum function manually
y = zeros(size(x));  % initialize
for i = 1:length(x)
    if x(i) > 0
        y(i) = 1;
    elseif x(i) < 0
        y(i) = -1;
    else
        y(i) = 0;
    end
end

% Plot
figure;
plot(x, y, 'LineWidth', 2);
grid on;
xlabel('x');
ylabel('sign(x)');
title('Signum Function');
ylim([-1.5 1.5]);
xline(0, '--r');  % reference line at x = 0
fprintf('\n')

% Exercise No:2
fprintf('EXERCISE Q2\n')

% Time vector
t = 0:0.1:10;   % from 0 to 10 seconds, step size 0.1

% Parameters
A = 1;          % Initial amplitude
k = 0.5;        % Growth rate (positive for growing signal)

% Exponentially growing signal
x = A * exp(k * t);

% Plot the signal
figure;
plot(t, x, 'LineWidth', 2);
grid on;
xlabel('Time (t)');
ylabel('Amplitude');
title('Exponentially Growing Signal:  x(t) = A * e^{k t}');
legend('x(t) = e^{0.5t}', 'Location', 'northwest');
fprintf('\n')

% Exercise No:3
fprintf('EXERCISE Q3\n')
% Time vector
t = 0:0.1:10;   % from 0 to 10 seconds

% Parameters
A = 1;          % Initial amplitude
k = 0.5;        % Decay rate (positive for decaying signal)

% Exponentially decaying signal
x = A * exp(-k * t);

% Plot the signal
figure;
plot(t, x, 'LineWidth', 2);
grid on;
xlabel('Time (t)');
ylabel('Amplitude');
title('Exponentially Decaying Signal:  x(t) = A * e^{-k t}');
legend('x(t) = e^{-0.5t}', 'Location', 'northeast');
fprintf('\n')

% Exercise No:4
fprintf('EXERCISE Q4\n')

% Time vector
t = 0:0.001:1;   % 1 second duration with fine resolution

% Signal parameters
A1 = 2;          % Amplitude of first signal
A2 = 1.5;        % Amplitude of second signal
f1 = 5;          % Frequency of first signal (Hz)
f2 = 3;          % Frequency of second signal (Hz)
phi1 = 0;        % Phase of first signal
phi2 = pi/4;     % Phase of second signal

% Define the two sinusoidal signals
x1 = A1 * sin(2*pi*f1*t + phi1);
x2 = A2 * sin(2*pi*f2*t + phi2);

% Subtract and multiply
x_sub = x1 - x2;
x_mul = x1 .* x2;   % element-by-element multiplication

% --- Plotting ---
figure;

subplot(3,1,1);
plot(t, x1, 'b', t, x2, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Sinusoidal Signals');
legend('x1(t)', 'x2(t)');

subplot(3,1,2);
plot(t, x_sub, 'm', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Subtracted Signal: x1(t) - x2(t)');

subplot(3,1,3);
plot(t, x_mul, 'g', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Multiplied Signal: x1(t) * x2(t)');
fprintf('\n')

% Exercise No:5
fprintf('EXERCISE Q5\n')

% Given data
alpha = 23e-6;     
L1 = 4.5;          
W1 = 2.25;         
T1_F = 400;        
T2_F = 920;        

% Temperature Conversion (F → C)
T1_C = FtoC(T1_F);
T2_C = FtoC(T2_F);
delta_T = T2_C - T1_C;   

% Thermal Expansion Calculations
delta_L = alpha * L1 * delta_T;   
delta_W = alpha * W1 * delta_T;   
L2 = L1 + delta_L;
W2 = W1 + delta_W;
A1 = L1 * W1;                     
A2 = L2 * W2;                     
delta_A = A2 - A1;                

% Display Results
fprintf('--- THERMAL EXPANSION RESULTS ---\n');
fprintf('Initial Temperature (°C): %.2f\n', T1_C);
fprintf('Final Temperature (°C): %.2f\n', T2_C);
fprintf('Temperature Change (°C): %.2f\n\n', delta_T);

fprintf('Initial Length (m): %.4f\n', L1);
fprintf('Final Length (m): %.4f\n', L2);
fprintf('Change in Length (m): %.6f\n\n', delta_L);

fprintf('Initial Width (m): %.4f\n', W1);
fprintf('Final Width (m): %.4f\n', W2);
fprintf('Change in Width (m): %.6f\n\n', delta_W);

fprintf('Initial Area (m^2): %.4f\n', A1);
fprintf('Final Area (m^2): %.4f\n', A2);
fprintf('Change in Area (ΔA): %.6f m^2\n', delta_A);

% Local Function
function Tc = FtoC(Tf)
    Tc = (Tf - 32) * 5/9;
end

% EXPERIMENT 7 - Corrected Exercises with printed explanations

%% Exercise No:1 - Sum of first 10 integers
fprintf('EXERCISE Q1\n')
disp('Explanation: Avoid using ''sum'' as a variable name to prevent conflict with MATLAB sum() function.')
sumVal = 0;  
for i = 1:10
    sumVal = sumVal + i;
end
disp(['Sum = ', num2str(sumVal)]);
fprintf('\n')

%% Exercise No:2 - Area of a rectangle
fprintf('EXERCISE Q2\n')
disp('Explanation: Width was missing, now defined as 3.')
lengthVal = 5;
width = 3;  
area = lengthVal * width;
disp(['Area = ', num2str(area)]);
fprintf('\n')

%% Exercise No:3 - Average of numbers
fprintf('EXERCISE Q3\n')
disp('Explanation: Removed square (^2) in denominator.')
x = [10 20 30 40];
avg = sum(x)/length(x);  
disp(['Average = ', num2str(avg)]);
fprintf('\n')

%% Exercise No:4 - Division operation
fprintf('EXERCISE Q4\n')
disp('Explanation: Avoid division by zero, changed b=1.')
a = 10; b = 1; 
c = a / b;
disp(['Result = ', num2str(c)]);
fprintf('\n')

%% Exercise No:5 - Factorial
fprintf('EXERCISE Q5\n')
disp('Explanation: Corrected function name from factorials() to factorial().')
n = 5;
f = factorial(n);  
disp(['Factorial = ', num2str(f)]);
fprintf('\n')

%% Exercise No:6 - Access array element
fprintf('EXERCISE Q6\n')
disp('Explanation: Index 5 was out of bounds, corrected to 4.')
A = [3 5 7 9];
val = A(4);  
disp(['Value = ', num2str(val)]);
fprintf('\n')

%% Exercise No:7 - While loop
fprintf('EXERCISE Q7\n')
disp('Explanation: Added increment i=i+1 to prevent infinite loop, corrected ''End'' to ''end''.')
i = 1;
while i <= 5
    disp(i);
    i = i + 1; 
end
fprintf('\n')

%% Exercise No:8 - Temperature conversion
fprintf('EXERCISE Q8\n')
disp('Explanation: Corrected variable name to C (was c).')
C = 25;
F = (9/5) * C + 32;  
disp(['Fahrenheit = ', num2str(F)]);
fprintf('\n')

%% Exercise No:9 - Function to square a number
fprintf('EXERCISE Q9\n')
disp('Explanation: Removed syntax error in function definition.')
sqVal = squareValue(5);  
disp(['Square = ', num2str(sqVal)]);
fprintf('\n')

%% Exercise No:10 - Rectangle area and perimeter function
fprintf('EXERCISE Q10\n')
disp('Explanation: Added semicolons and function call corrected.')
[lArea, lPeri] = rectangleCalc(4, 6); 
disp(['Area = ', num2str(lArea) ', Perimeter = ' num2str(lPeri)]);
fprintf('\n')

%% Exercise No:11 - Sum of first n integers
fprintf('EXERCISE Q11\n')
disp('Explanation: Loop now includes n to sum correctly.')
n = 10;
s = 0;
for i = 1:n  
    s = s + i;
end
disp(['Sum = ', num2str(s)]);
fprintf('\n')

%% Exercise No:12 - Quadratic formula
fprintf('EXERCISE Q12\n')
disp('Explanation: Added parentheses around 2*a in denominator.')
a = 1; b = 5; c = 6;
root1 = (-b + sqrt(b^2 - 4*a*c)) / (2*a); 
root2 = (-b - sqrt(b^2 - 4*a*c)) / (2*a);
disp(['Roots = ' num2str(root1) ', ' num2str(root2)]);
fprintf('\n')

%% Exercise No:13 - Access element in matrix
fprintf('EXERCISE Q13\n')
disp('Explanation: Corrected matrix index to be within bounds.')
M = magic(4);
valM = M(4,2);  
disp(['Value = ', num2str(valM)]);
fprintf('\n')

%% Exercise No:14 - Average of vector function
fprintf('EXERCISE Q14\n')
disp('Explanation: Corrected total = total + x(i), not x.')
vec = [2 4 6 8];
avgVec = myAverage(vec);
disp(['Average = ', num2str(avgVec)]);
fprintf('\n')

%% Exercise No:15 - Recursive Fibonacci
fprintf('EXERCISE Q15\n')
disp('Explanation: Corrected recursion from fib(n-3) to fib(n-2).')
fib5 = fib(5);  
disp(['5th Fibonacci = ', num2str(fib5)]);
fprintf('\n')

%% Exercise No:16 - String concatenation
fprintf('EXERCISE Q16\n')
disp('Explanation: Converted number to string using num2str().')
age = 21;
name = 'Alice';
disp(['Name: ' name ', Age: ' num2str(age)]); 
fprintf('\n')

%% Exercise No:17 - Smallest k such that sum(1:k) > 50
fprintf('EXERCISE Q17\n')
disp('Explanation: Adjusted k-1 after last increment.')
k = 1; sSum = 0;
while sSum <= 50
    sSum = sSum + k;
    k = k + 1;
end
disp(['k = ', num2str(k-1)]); 
fprintf('\n')

%% Exercise No:18 - Multiplication table
fprintf('EXERCISE Q18\n')
disp('Explanation: Nested loops used correctly to fill table.')
for i = 1:5
    for j = 1:5
        T(i,j) = i*j;
    end
end
disp(T);

%% ===== Functions =====
function y = squareValue(x)
y = x*x;
end

function [area, peri] = rectangleCalc(l, w)
area = l*w;
peri = 2*(l+w);
end
function avg = myAverage(x)
n = length(x);
total = 0;
for i = 1:n
    total = total + x(i); 
end
avg = total / n;
end
function f = fib(n)
if n == 1 || n == 2
    f = 1;
else
    f = fib(n-1) + fib(n-2); 
end
end