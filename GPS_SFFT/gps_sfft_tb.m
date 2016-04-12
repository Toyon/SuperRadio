%--------------------------------------------------------------------------
% sparse FFT GPS function test-bench for HDL coder
% Toyon Research 2016
% Author: B.Weals
%--------------------------------------------------------------------------
clear all;

% Parameters: start
Ncode = 1024;
SNR_dB = -15;
offset = 20000; % code offset
n = 2^16; % number of samples to process
% Parameters: end

p = log2(n);
B = n/p;
Nsamp_per_chip = n/Ncode;

code_orig = 2*((rand(1,Ncode)>0.5)-0.5);
%code_orig(1) = 1;
%code_orig(2) = -1;
code = ones(Nsamp_per_chip,1)*code_orig;
code = code(:)';
fft_code = fft(code);
sub_sampled_fcode = fft_code(1:p:end);
conj_sub_sampled_fcode = conj(sub_sampled_fcode);

signal = ones(Nsamp_per_chip,1)*code_orig;
signal = signal(:)';
signal = circshift(signal,offset,2);
uncorrupted_signal = signal;
noise = (10^(-SNR_dB/20))*randn(size(signal))/sqrt(2);
noise = noise + 1i*(10^(-SNR_dB/20))*randn(size(signal))/sqrt(2);
signal = signal + noise;

x = reshape(signal',B,p)';
x1 = sum(x);
fx1 = fft(x1);

cfx1 = conj_sub_sampled_fcode.*fx1;
fft_code = [real(conj_sub_sampled_fcode)' imag(conj_sub_sampled_fcode)'];
save('fft_code.txt','fft_code','-ascii','-double');
fft_code = conj_sub_sampled_fcode;
save('fft_code.mat','fft_code');
x2 = ifft(cfx1);
mx2 = abs(x2);
[mm,k] = max(mx2);

figure(1),plot(mx2);
cor_results = zeros(1,p);
code_to_save = code(1:p:end)'; %
save('code.txt','code_to_save','-ascii');
code0 = code;
code = code(1:p:end);
save('code.mat','code');
code = code0;
code0 = [];

for m=1:p
    shift_code = circshift(code,k+(m-1)*B-1,2);
    sub_shift_code = shift_code(1:p:end);
    cor_results(m) = abs(sum(sub_shift_code.*signal(1:p:end)));
end

figure(2),plot(cor_results);
[mm,max_cor_idx] = max(cor_results);

snr = 10*log10(var(uncorrupted_signal)/var(noise))
offset
offset_estimate = (max_cor_idx-1)*B+k-1
err = offset - offset_estimate
chip_error = err/Nsamp_per_chip


% Now we'll do the hardware:
validOut = 0;
validIn = 1;
m = 1;
while ~validOut
    [timeOffset,validOut] = gps_sfft(real(signal(m)), imag(signal(m)), validIn);
    m = m+1;
    if m > n
        m = 1;
        validIn = 0;
    end
end

display(['offset: estimated: ' num2str(timeOffset) ' actual: ' num2str(20000)]);







