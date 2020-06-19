function [a,b]=levy_estimation(eps_data,s,p,q)

% Levy Estimation of poles for permittivity data fitting
%
% In general, the output of this function will then be fed in to a more
% sophisticated optimization process (e.g. Vector Fitting or Mathematica)
%
% This is done to have better initial conditions for the iterative solver
% to work with.
%
% eps_data is a column vector of complex data representing e(w) = e_re + e_im
%   size: 1xNs
% s is vector of frequency points
%   size: 1xNs
% p is number of zeros to use
% q is number of poles to use (these are nominally guesses)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code only works accurately for p = q currently!!!!!%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

re_data = real(eps_data);
im_data = imag(eps_data);

[~, Ns] = size(eps_data);

lambda = zeros(1,2*p);
S = zeros(1,p+q);
T = zeros(1,p+q);
U = zeros(1,2*q);

for h = 1:2:2*p+1
    lambda(h) = sum(s.^(h-1));
%     lambda(h) = h-1;
end

for h = 1:2:2*(p+q)+1
    S(h) = sum( (s.^(h-1)).*re_data );
%     S(h) = h-1;
end

for h = 1:2:2*(p+q)+1
    T(h) = sum( (s.^(h)).*im_data );
%     T(h) = h;
end

for h = 3:2:2*q+1
    U(h) = sum( (s.^(h-1)) .* (re_data.^2 + im_data.^2) );
%     U(h) = h-1;
end

size_pq = 2*floor((p+q-1)/2) + 1; % as we want nearest lower odd integer
min_size = (size_pq + 1)/2;



M_ul = zeros(p+1);
M_ur = zeros(p+1,q);
M_ll = zeros(q,p+1);
M_lr = zeros(q);

% Iterate over diagonals:
%% First Upper Left Matrix (M_ul):
M_size = size(M_ul,1); % For computing linear index for sub/super diagonals
for d = 1:2:p+1
    sup_diagonal = zeros(p-d+2,1);
    sub_diagonal = zeros(p-d+2,1);
    index = 1;
    for h = d:2:2*p+2-d
        tri = (-1)^( ((index)^2 + index + 2)/2 ); % Triangular sequence to get 1,1,-1,-1,1,1...
        sup_diagonal(index) = (-1)^((d-1)/2)*tri*lambda(h);
        sub_diagonal(index) = tri*lambda(h);
        index = index + 1;
    end
    if d == 1
        M_ul = diag(sup_diagonal);
    else
        M_ul((d-1)*M_size+1:M_size+1:end) = sup_diagonal; % super-diagonal elements
        M_ul(d:M_size+1:1+M_size*(M_size-(d-1))) = sub_diagonal; % sub-diagonal
    end
end

%% Upper Right Matrix (M_ur):

% Do Super-Diagonal Terms first:
M_size = size(M_ur,1); % For computing linear index for sub/super diagonals
for d = 1:q
    sup_diagonal = zeros(q-d+1,1);
    index = 1;
    for h = d:2:2*(q)-d
        tri = (-1)^( ((index)^2 + index + 2)/2 ); % Triangular sequence to get 1,1,-1,-1,1,1...
        if mod(d,2) ~= 0
            sup_diagonal(index) = (-1)^((d-1)/2)*tri*T(h);
        else
            sup_diagonal(index) = (-1)^((d/2)+1)*tri*S(h+1);
        end
        index = index + 1;
    end
    if d == 1
        M_ur(1:M_size+1:end) = sup_diagonal; % Actually the 'central' diagonal
    else
        M_ur((d-1)*M_size+1:M_size+1:end) = sup_diagonal; % super-diagonal elements
    end
end

% Then sub-diagonal:
for d = 2:p+1
    sub_diagonal = zeros(p+1-d+1,1);
    index = 1;
    for h = d:2:2*(p+1)-d
        tri = (-1)^( ((index)^2 + index + 2)/2 ); % Triangular sequence to get 1,1,-1,-1,1,1...
        tri_sub = (-1)^( ((index+3)^2 + index+3 + 2)/2 );
        if mod(d,2) ~= 0
            sub_diagonal(index) = tri*T(h);
        else
            sub_diagonal(index) = tri_sub*S(h+1);
        end
        index = index + 1;
    end
    M_ur(d:M_size+1:1+M_size*(min(M_size-(d-1),size(M_ur,2))) ) = sub_diagonal; % sub-diagonal
end

%% Lower Left Matrix (M_ll):

% Do Super-Diagonal Terms first:
M_size = size(M_ll,1); % For computing linear index for sub/super diagonals
for d = 2:p+1
    sup_diagonal = zeros(p+1-d+1,1);
    index = 1;
    for h = d:2:2*(p+1)-d
        tri = (-1)^( ((index)^2 + index + 2)/2 ); % Triangular sequence to get 1,1,-1,-1,1,1...
        if mod(d,2) ~= 0
            sup_diagonal(index) = (-1)^((d-1)/2)*tri*T(h);
        else
            sup_diagonal(index) = (-1)^(d/2)*tri*S(h+1);
        end
        index = index + 1;
    end
    M_ll((d-1)*M_size+1:M_size+1:end) = sup_diagonal; % super-diagonal elements
end

% Then sub-diagonal:
for d = 1:q
    sub_diagonal = zeros(q-d+1,1);
    index = 1;
    for h = d:2:2*(q)-d
        tri = (-1)^( ((index)^2 + index + 2)/2 ); % Triangular sequence to get 1,1,-1,-1,1,1...
        tri_sub = (-1)^( ((index+1)^2 + index+1 + 2)/2 );
        if mod(d,2) ~= 0
            sub_diagonal(index) = tri*T(h);
        else
            sub_diagonal(index) = tri_sub*S(h+1);
        end
        index = index + 1;
    end
    if d == 1
        M_ll(1:M_size+1:end) = sub_diagonal; % Actually the 'central' diagonal
    else
        M_ll(d:M_size+1:1+M_size*(min(M_size-(d-1),size(M_ll,2))) ) = sub_diagonal; % sub-diagonal
    end
end

%% Lower Right Matrix (M_lr):
M_size = size(M_lr,1); % For computing linear index for sub/super diagonals
for d = 3:2:2*q-1
    sup_diagonal = zeros(q-d+3,1);
    sub_diagonal = zeros(q-d+3,1);
    index = 1;
    for h = d:2:2*q+4-d
        tri = (-1)^( ((index)^2 + index + 2)/2 ); % Triangular sequence to get 1,1,-1,-1,1,1...
        sup_diagonal(index) = (-1)^((d+1)/2)*tri*U(h);
        sub_diagonal(index) = tri*U(h);
        index = index + 1;
    end
    if d == 3
        M_lr(1:M_size+1:end) = sup_diagonal; % Actually the 'central' diagonal
    else
        M_lr((d-3)*M_size+1:M_size+1:end) = sup_diagonal; % super-diagonal elements
        M_lr(d-2:M_size+1:1+M_size*(M_size-(d-3))) = sub_diagonal; % sub-diagonal
    end
end


%% Put it all together and Solve
% size(M_ul)
% size(M_ur)
% size(M_ll)
% size(M_lr)

% M_ul
% M_ur
% M_ll
% M_lr

% The resizing is needed here to make sure we can actually assemble this
% matrix. It might be better(?) to pad with 0's to expand the smaller
% matrices rather than clip the larger one...
total_M = [M_ul M_ur; M_ll M_lr];
cu_vec = zeros(max(p,q)+1,1);

for h = 1:max(p,q)+1
    if mod(h,2) ~= 0
        cu_vec(h) = S(h);
    else
        cu_vec(h) = T(h-1);
    end
end


c1_vec = U';
c1_vec = c1_vec(2:q+1);

c_vec = [cu_vec; c1_vec];

% size(total_M)
% size(c_vec)

n = total_M\c_vec;

a = n(1:p+1);
b = n(p+2:end);




