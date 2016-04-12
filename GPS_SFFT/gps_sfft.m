%--------------------------------------------------------------------------
% sparse FFT GPS function for HDL coder
% Toyon Research 2016
% Author: B.Weals
%--------------------------------------------------------------------------
function [timeOffset,validOut] = gps_sfft(yIn_Re, yIn_Im, validIn)

% fft objects & related
persistent fft4096 ifft4096
persistent alias_buf fft_buf mult_buf fft_code code corr_buf samp_buf
persistent ifft_buf
persistent ka kf km kc kp ks ks0 ksb ksc
persistent max_bin imax_bin max_corr imax_corr
% states:
persistent state

% States:
do_alias = 0;
do_fft = 1;
do_mult = 2;
do_ifft = 3;
do_bin = 4;
do_corr = 5;
do_offset = 6;

% constants:
n = 2^16; % 65536
p = 2^4;  % 16
B = 2^12; % 4096

% Block to set persistent variables:
if isempty(fft4096)

    samp_buf = complex(zeros(1,65536));
    alias_buf = complex(zeros(1,4096));
    fft_buf = complex(zeros(1,4096));
    ifft_buf = complex(zeros(1,4096));
    mult_buf = complex(zeros(1,4096));
    %tmp1 = coder.load('fft_code.txt','-ascii');
    %fft_code = [1 1i]*tmp1';
    fft_code = coder.load('fft_code.mat','fft_code');
    %coder.load('fft_code.mat');
    %tmp = coder.load('code.txt','-ascii');
    %code = tmp';
    code = coder.load('code.mat','code');
    %coder.load('code.mat','code');
    corr_buf = complex(zeros(1,p));
    
    ka = 1;
    kf = 1;
    km = 1;
    kc = 1;
    kp = 1;
    ks = 1;
    ks0 = 1;
    ksb = 1;
    ksc = 0;
    
    max_bin = 0;
    imax_bin = 0;
    max_corr = 0;
    imax_corr = 0;
    
    state = 0;
    fft4096 = dsp.HDLFFT('FFTLength',4096,'OverflowAction','Saturate',...
        'RoundingMethod', 'Nearest','BitReversedOutput',false);
    ifft4096 = dsp.HDLIFFT('FFTLength',4096,'OverflowAction','Saturate',...
        'RoundingMethod', 'Nearest','BitReversedOutput',false);
end

% Initialize output variables:
timeOffset = 0;
validOut = 0;

% Load sample buffer and alias samples in time:
if state == do_alias
    % Alias in time:
    if kp <= p && validIn
        samp = complex(yIn_Re,yIn_Im);
        samp_buf(ks) = samp;
        ks = mod(ks,n)+1;
        if kp == 1
            alias_buf(ka) = samp;
        else
            alias_buf(ka) = alias_buf(ka) + samp;
        end
        if ka == B
            ka = 1; 
            kp = kp+1;
        else
            ka = ka+1;
        end
    else 
        state = do_fft;
        ka = 1;
        kp = 1;
        %save('alias_buf.mat','alias_buf');
    end
end

% Take FFT of aliased samples:
fidx_a = 1;
fvalid = false;
if state == do_fft 
    fvalid = ka <= B;
    if fvalid
        fidx_a = ka;
    else
        fidx_a = 1;
    end
end
% Note that the streaming FFT call cannot be inside a conditional block or
% a loop:
[yOut,fft_valid] = step(fft4096,alias_buf(fidx_a),fvalid);
if state == do_fft
    %[state ka kf fft_valid]
    if fvalid
        ka = ka+1;
    end

    if kf <= B && fft_valid
        fft_buf(kf) = yOut;
        kf = kf+1;
    elseif kf > B
        kf = 1;  
        ka = 1;
        state = do_mult;
        %save('fft_buf.mat','fft_buf');
    end            
end

% Multiply FFT of aliased samples by down-sampled FFT of code:
if state == do_mult
    if km <= B
        mult_buf(km) = fft_code.fft_code(km)*fft_buf(km);
        km = km+1;
    else
        km = 1;
        state = do_ifft;
        %save('mult_buf.mat','mult_buf');
    end
end

% Take the IFFT of the result of the multiplication:
iidx_a = 1;
ivalid = false;
if state == do_ifft    
    ivalid = ka <= B;
    if ivalid
        iidx_a = ka;
    else
        iidx_a = 1;
    end
end
[iyOut,ifft_valid] = step(ifft4096,mult_buf(iidx_a),ivalid);       
if state == do_ifft
    if ivalid
        ka = ka+1;
    end      
    if kf <= B && ifft_valid
        ifft_buf(kf) = iyOut;
        kf = kf+1;
    elseif kf > B
        kf = 1;
        ka = 1;
        state = do_bin;
        max_bin = 0;
        imax_bin = 0;            
        %save('ifft_buf.mat','ifft_buf');
    end
end
   
% Find the maximum magnitude bin (IFFT output):
if state == do_bin 
    if kf <= B
        rb = abs(real(ifft_buf(kf)));
        ib = abs(imag(ifft_buf(kf)));
        bin_mag = max(rb,ib) + 0.375*min(rb,ib); % abs(fft_buf(kf));
        if bin_mag > max_bin
            max_bin = bin_mag;
            imax_bin = kf;
        end
        kf = kf+1;
    else
        kf = 1;
        state = do_corr;
        ks0 = imax_bin;
        ksc = imax_bin-1;
        %figure(3),plot(abs(fft_buf));
    end
end

% Do p time correlations to associated with the maximum bin:
if state == do_corr
    if kp <= p
        if kc == 1
            corr_buf(kp) = samp_buf(ksc+1)*code.code(kc);
        else
            corr_buf(kp) = corr_buf(kp) + samp_buf(ksc+1)*code.code(kc);
        end
        ksc = mod(ksc+p,n);
        ksb = ksb+p;
        if kc < B
            kc = kc+1;
        else
            kc = 1;
            kp = kp+1;
            ks0 = ks0+B;
            ksc = ks0-1;
            ksb = 1;
        end
    else
        kp = 1;
        ksc = 0;
        ksb = 1;
        kc = 1;
        state = do_offset;
        %figure(4),plot(abs(corr_buf));
    end
end

% Find the maximum magnitude time correlation and report result:
if state == do_offset
    if kp <= p
        rc = abs(real(corr_buf(kp)));
        ic = abs(imag(corr_buf(kp)));
        mag_corr = max(rc,ic) + 0.375*min(rc,ic); % abs(corr_buf(kp));                
        if mag_corr > max_corr
            max_corr = mag_corr;
            imax_corr = kp;
        end
        kp = kp+1;
    else
        kp = 1;      
        timeOffset = (imax_corr-1)*B + imax_bin - 1;
        max_corr = 0;
        imax_corr = 0;          
        validOut = 1;
    end
end     
%[state validOut]        
        
        
    
    
    
        
    


        


