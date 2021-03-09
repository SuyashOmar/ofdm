clear
N = 8192; % number of bits or symbols
Eb_N0_dB = [0:10]; 
fftsize = 32;
cp = 2;
pkg load communications

for ii = 1:length(Eb_N0_dB)
    % Transmitter
    ip1 = rand(1,N)>0.5; % generating 0,1 with equal probability
    ip2 = rand(1,N)>0.5;
    ip = ip1+j*ip2;
    s1 = 2*ip1-1; % QPSK modulation 0
    s2 = 2*ip2-1;
    s = s1+j*s2;
    a = cp*length(s)/fftsize;
    xti = ofdm(s, fftsize, cp);
    xti1 = transpose(xti);
    xtil = reshape(xti1,1,[]);
    n = 1/sqrt(2)*[randn(1,(N+a)) + j*randn(1,(N+a))];
    y = xtil + 10^(-Eb_N0_dB(ii)/20)*n;
    %y = awgn(xtil, Eb_N0_dB(ii));
    yti = dofdm(y, fftsize, cp);
    yti1 = transpose(yti);
    ytil = reshape(yti1,1,[]);
    iprHat = real(ytil)>0;
    ipiHat = imag(ytil)>0;
    ipHat = iprHat+j*ipiHat;
    nErr(ii) = size(find([ip- ipHat]),2);
end
function ytil = ofdm(x, fftsize, cp)
  ytil = zeros(length(x)/fftsize, fftsize+cp,1);
  j = 1;
  for i = 1:fftsize:length(x)
    jp = x(i:i+fftsize-1);
    jp1 = ifft(jp,fftsize);  
    ytil(j,:) = [jp1(length(jp1)-cp+1:length(jp1)) jp1]; % cyclic prefix size.
    j++;
  endfor
endfunction
function ybar = dofdm(x, fftsize, cp)
  a = fftsize+cp;
  j = 1;
  ybar = zeros(length(x)/a, fftsize, 1);
  for i = 1:a:length(x)
    b = x(i+cp:a+i-1);
    ybar(j,:) = fft(b,fftsize);
    j++;
  endfor
endfunction  
simBer = nErr/N; % simulated ber
EbN0Lin = 10.^(Eb_N0_dB/10);
close all
figure
semilogy(Eb_N0_dB,simBer,'mo-');
grid on
legend('Ofdm with 32-FFT and cp=2');
xlabel('Average Eb/No,dB');
ylabel('Bit Error Rate');

  