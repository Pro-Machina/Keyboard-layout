%% Assignment 2 | Kshitij Pisal | 16QM30001
% (date of submission - 9th April 2019)

%% Improting files
bigram = importdata('bigram-probabilities.txt'); % In the data given in the question, this is the first array of 27*27 numbers 
moving_time = importdata('inter-key-intervals.txt'); % In the data given in the question, this is the second array of 27*27 numbers 
alphabet = ['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'i'; 'j'; 'k'; 'l'; 'm'; 'n'; 'o'; 'p'; 'q'; 'r'; 's'; 't'; 'u'; 'v'; 'w'; 'x'; 'y'; 'z'; '.'];
% The 27th letter is denoted by '.'

%% Simulated Annealing

temperature = 100; %Initially set to 100 (please change to required value)
rho = 0.9; %Factor by which temperature is reduced every iteration (please change to required value)
markov_len = 50; %Markov length is set to 50 (please change to required value)
      
while (temperature > 0.1) % Limit to the minimum temperature reached is kept 0.1
    mar = 0; % Markov length is again initiated to zero after one iteration
    
    while (mar < markov_len)
        
    % Initial estimate of MT(i,j) is given as e1
        for i = 1:27
            for j = 1:27
                for k = 1:27
                    for l = 1:27
                        e1 = moving_time(i, j)*bigram(k, l); % The nested for loop is for the esimated value Sum(Sum(Sum(Sum(MT[i,j] * P[k,l]))))
                        % i,j is the slot time in miliseconds, k,l is the
                        % bigram probability between two alphbet, say k and
                        % l
                    end
                end
            end
        end
  
        % Applying swapping and moving operations
        
        % Cumiliative probs:
        % 1 swap = 0.5
        % 2 swap = 0.3
        % 1 move = 0.15
        % 2 move = 0.05
        % Functions of the operations are defined at the end
        check_opr = 0.01*randi([0 100]); % Random number from 0.00 to 1.00
        if(check_opr < 0.5)
            [bi_temp, mo_temp, al_temp] = swap_once(moving_time, bigram, alphabet);
        elseif(check_opr >= 0.5 && check_opr < 0.8)
            [bi_temp, mo_temp, al_temp] = swap_twice(moving_time, bigram, alphabet);
        elseif(check_opr >= 0.8 && check_opr < 0.95)
            [bi_temp, mo_temp, al_temp] = move_once(moving_time, bigram, alphabet);
        elseif(check_opr >= 0.95 && check_opr < 1)
            [bi_temp, mo_temp, al_temp] = move_twice(moving_time, bigram, alphabet);
        end
        
        % For calculating the estimate of the swaped/moved array of
        % alphabet
        for i = 1:27
            for j = 1:27
                for k = 1:27
                    for l = 1:27
                        e2 = mo_temp(i, j)*bi_temp(k, l);
                    end
                end
            end
        end
        
        % Probability calculated from Arrhenius equation
        p_accept = anne_prob(e1, e2, temperature); % Function defined at the end
        check_prob = 0.01*randi([0 100]); 
        if(p_accept >= check_prob) % If the solution is accepted, the node of the markov chain is moved to the new solution
            moving_time = mo_temp;
            bigram = bi_temp;
            alphabet = al_temp;
            mar = mar + 1;
        end
    end
    temperature = temperature*rho; % At the end of the markov chain, the temperature is reduced by a factor of rho
end
        

%% Functions used (Swaping, Moving, Probability calculation)
function [ m1, m2, a1 ] = swap_once(matrix1, matrix2, alpha) % return array of three arrays corresponding to new MT[i,j], P[k,l] and their corresponding alphabet
    
    ind = randi([1 26]); % RAndom number between 1 and 26
    m1 = matrix1;
    row1 = m1(ind, :);
    col1 = m1(:, ind);
    m1(ind, :) = m1(ind+1, :);
    m1(:, ind) = m1(:, ind+1);
    m1(ind+1, :) = row1;
    m1(:, ind+1) = col1;
    
    m2 = matrix2;
    row2 = m2(ind, :);
    col2 = m2(:, ind);
    m2(ind, :) = m2(ind+1, :);
    m2(:, ind) = m2(:, 9);
    m2(ind+1, :) = row2;
    m2(:, ind+1) = col2;
    
    a1 = alpha;
    row = a1(ind);
    a1(ind) = a1(ind+1);
    a1(ind+1) = row;
    
end

function [ m1, m2, a1 ] = swap_twice(matrix1, matrix2, alpha) % Same as above function, just that instead of one there are two swaps (at random positions).
    ind1 = randi([1 26]);
    ind2 = randi([1 26]);
    while (ind1 == ind2)
        ind2 = randi([1 26]);
    end
 
% swap 1    
    m1 = matrix1;
    row1 = m1(ind1, :);
    col1 = m1(:, ind1);
    m1(ind1, :) = m1(ind1+1, :);
    m1(:, ind1) = m1(:, ind1+1);
    m1(ind1+1, :) = row1;
    m1(:, ind1+1) = col1;
    
    m2 = matrix2;
    row2 = m2(ind1, :);
    col2 = m2(:, ind1);
    m2(ind1, :) = m2(ind1+1, :);
    m2(:, ind1) = m2(:, 9);
    m2(ind1+1, :) = row2;
    m2(:, ind1+1) = col2;
 
    a1 = alpha;
    row = a1(ind1);
    a1(ind1) = a1(ind1+1);
    a1(ind1+1) = row;

% swap 2    
    row1 = m1(ind2, :);
    col1 = m1(:, ind2);
    m1(ind2, :) = m1(ind2+1, :);
    m1(:, ind2) = m1(:, ind2+1);
    m1(ind2+1, :) = row1;
    m1(:, ind2+1) = col1;
    
    row2 = m2(ind2, :);
    col2 = m2(:, ind2);
    m2(ind2, :) = m2(ind2+1, :);
    m2(:, ind2) = m2(:, 9);
    m2(ind2+1, :) = row2;
    m2(:, ind2+1) = col2;
    
    row = a1(ind2);
    a1(ind2) = a1(ind2+1);
    a1(ind2+1) = row;
end

function [m1, m2, a1] = move_once(matrix1, matrix2, alpha) % Algorithm similar to swap, insted of consicutive, directly opposite blocks are choosen
    
    ind = randi([1 26]);
    m1 = matrix1;
    row1 = m1(ind, :);
    col1 = m1(:, ind);
    m1(ind, :) = m1(27-ind, :);
    m1(:, ind) = m1(:, 27-ind);
    m1(27-ind, :) = row1;
    m1(:, 27-ind) = col1;
    
    m2 = matrix2;
    row2 = m2(ind, :);
    col2 = m2(:, ind);
    m2(ind, :) = m2(ind+1, :);
    m2(:, ind) = m2(:, 9);
    m2(ind+1, :) = row2;
    m2(:, ind+1) = col2;
    
    a1 = alpha;
    row = a1(ind);
    a1(ind) = a1(ind+1);
    a1(ind+1) = row;
    
end  

function [m1, m2, a1] = move_twice(matrix1, matrix2, alpha)
    ind1 = randi([1 26]);
    ind2 = randi([1 26]);
    while (ind1 == ind2)
        ind2 = randi([1 26]);
    end
 
% move 1    
    m1 = matrix1;
    row1 = m1(ind1, :);
    col1 = m1(:, ind1);
    m1(ind1, :) = m1(27-ind1, :);
    m1(:, ind1) = m1(:, 27-ind1);
    m1(27-ind1, :) = row1;
    m1(:, 27-ind1) = col1;
    
    m2 = matrix2;
    row2 = m2(ind1, :);
    col2 = m2(:, ind1);
    m2(ind1, :) = m2(27-ind1, :);
    m2(:, ind1) = m2(:, 9);
    m2(27-ind1, :) = row2;
    m2(:, 27-ind1) = col2;
 
    a1 = alpha;
    row = a1(ind1);
    a1(ind1) = a1(27-ind1);
    a1(27-ind1) = row;

% move 2    
    row1 = m1(ind2, :);
    col1 = m1(:, ind2);
    m1(ind2, :) = m1(27-ind2, :);
    m1(:, ind2) = m1(:, 27-ind2);
    m1(27-ind2, :) = row1;
    m1(:, 27-ind2) = col1;
    
    row2 = m2(ind2, :);
    col2 = m2(:, ind2);
    m2(ind2, :) = m2(27-ind2, :);
    m2(:, ind2) = m2(:, 9);
    m2(27-ind2, :) = row2;
    m2(:, 27-ind2) = col2;
    
    row = a1(ind2);
    a1(ind2) = a1(27-ind2);
    a1(27-ind2) = row;
end

function p_accept = anne_prob(e1, e2, t)
    del_e = (e2 - e1);
    if(del_e > 0)
        p_accept = exp(-(del_e)/t); % Probability of acceptance according to Arrhenius equation 
    else
        p_accept = 1; % if e2 < e1 then the solution is more favourable and the probability of acceptance is automatically 1
    end
end
