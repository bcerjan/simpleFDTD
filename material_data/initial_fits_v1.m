% Script to perfrom e(w) fitting for given input data
%
% Assumes the input spectra is in the format: lambda (um), n, k
%

% Set some inital parameters:
% Max / Min wavelength (um)
l_min = 0.2; 
l_max = 1.0;
speedC = 2.99792458e14; % c in um/sec
h = 4.135667e-15; % Planck constant in eV/sec

% Number of poles / zeros:
num_poles = 9; % Needs to be two times + 1 what is used in Mathematica as
               % complex conjugate-pairs are individual poles in this 
               % description
p = num_poles;
q = num_poles;

% Read in spectral data:
fname = 'Ag-Wu.txt'; % whitespace delimited

data = importdata(fname, ' ', 0); % 0 Header rows

% Trim to our range of interest:
data = data(data(:,1) >= l_min & data(:,1) <= l_max,:);

f = zeros(size(data,1),2);

% Convert to eV and condense to a single complex number:
f(:,1) = h * speedC ./ data(:,1);
f(:,2) = data(:,2).^2 - data(:,3).^2 + 2*i*data(:,2).*data(:,3);

% Now we're ready for the fitting
% First, Levy method to find good starting points:
[a, b] = levy_estimation(f(:,2),f(:,1),p,q);

% Add b0 term (always = 1):
b = [1; b];

% Now we need to invert the order of a and b to use Matlab's residue
% function
b = flip(b);
a = flip(a);


% Now turn these into the poles we feed in to the VF step:
[res, poles, k] = residue(a,b)

% Plot the results at this point -- they shouldn't be terrible

s = 1i.*f(:,1)';
test_fit = zeros(1,length(s));
test_fit = k; % Set constant term
for n = 1:length(poles)
    test_fit = test_fit + res(n) ./ (s - poles(n));
end

figure(3)
plot(imag(s),imag(test_fit),imag(s),-imag(f(:,2)'))
figure(4)
plot(imag(s),real(test_fit),imag(s),real(f(:,2)'))
