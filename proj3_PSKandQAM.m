% clear all;
% close all;

num_bits=12000;
Eb=1;
SNR_vector= -4:14;
%Generate random bits 
Data = randi([0 1] , 1 , num_bits);
%Define positions on signal space
BPSK_positions=[1 -1];
QPSK_positions=[-1-1i -1+1i 1-1i 1+1i];
QPSK_positions_NotGray=[-1-1i -1+1i 1+1i 1-1i];
PSK8_positions=[1+0i sqrt(1/2)+sqrt(1/2)*1i -sqrt(1/2)+sqrt(1/2)*1i 0+1i sqrt(1/2)-sqrt(1/2)*1i 0-1i -1+0i -sqrt(1/2)-sqrt(1/2)*1i];
QAM16_positions=[-3-3i -3-1i -3+3i -3+1i -1-3i -1-1i -1+3i -1+1i 3-3i 3-1i 3+3i 3+1i 1-3i 1-1i 1+3i 1+1i];

BPSK(Data,num_bits,Eb,SNR_vector,BPSK_positions);
QPSK(Data,num_bits,Eb,SNR_vector,QPSK_positions,QPSK_positions_NotGray);
QAM(Data,num_bits,Eb,SNR_vector,QAM16_positions);
PSK8(Data,num_bits,Eb,SNR_vector,PSK8_positions);

        %---------------BPSK-----------------%
function BPSK(Data,num_bits,Eb,SNR_vector,BPSK_positions)
bits_per_symbol=1;
%%%the mapper%%%%
%mapping to BPSK (convert logic 0 to -1)
mapper_out = mapp_function(Data,num_bits,'BPSK',BPSK_positions,bits_per_symbol);

%%%The Channel%%%%
%Generate a unity variance, zero mean additive white Gaussian noise signal
%with the same size as transmitted signal.
noise=randn(size(mapper_out));
%Get the Noisy signal for each SNR value
for i=1:length(SNR_vector)
     %calculate variance from SNR (SNR=Eb/No)
     No_vector(i)=Eb / 10^(SNR_vector(i)/10);
     variance_vector(i)=No_vector(i)/2; %variance=No/2
     %Scale the noise sequence to have variance = N0/2 by multiplying the sequence
     %by sqrt(N0/2). 
     scaled_noise=sqrt(variance_vector(i)) * noise; 
     %Add the noise to the transmitted sequence
     noisy_signal_BPSK(i,:)=mapper_out + scaled_noise;
end
%%%%The Demapper and BER calculation%%%
BER=get_BER(SNR_vector,noisy_signal_BPSK,num_bits,Data,'BPSK',bits_per_symbol,BPSK_positions); 
%calculate the theoritical BER
BER_theoretical=0.5 * erfc(sqrt(Eb ./ No_vector));
%plot Demapper output BER vs theoretical
plot_BER(SNR_vector,BER,BER_theoretical,'BPSK - practical BER Vs theoritical BER');
end
        %---------------QPSK-----------------%
function QPSK(Data,num_bits,Eb,SNR_vector,QPSK_positions,QPSK_positions_NotGray)
bits_per_symbol=2;
%mapping to QPSK
mapper_out = mapp_function(Data,num_bits,'QPSK',QPSK_positions,bits_per_symbol);
mapper_out_NotGray = mapp_function(Data,num_bits,'QPSK',QPSK_positions_NotGray,bits_per_symbol);

%%%%The Channel%%%%
%Generate a unity variance, zero mean additive white Gaussian noise signal
%with the same size as transmitted signal.
I_noise=randn(size(mapper_out));
Q_noise=randn(size(mapper_out));
%Get the Noisy signal for each SNR value
Eavg=sum(abs(mapper_out).^2)/(num_bits/bits_per_symbol);
for j=1:length(SNR_vector)
     %calculate variance from SNR (SNR=Es/Ns)
     No_vector(j)=Eb / 10^(SNR_vector(j)/10);
     variance_vector(j)=(No_vector(j)*(Eavg/bits_per_symbol))/2; %variance=Ns/2
     %Scale the noise sequence to have variance = N0/2 by multiplying the sequence
     %by sqrt(N0/2). 
     I_scaled_noise=sqrt(variance_vector(j)) * I_noise; 
     Q_scaled_noise=sqrt(variance_vector(j)) * Q_noise; 
     %Add the noise to the transmitted sequence
     noisy_signal_QPSK(j,:)=mapper_out + I_scaled_noise + Q_scaled_noise*1i;
     noisy_signal_QPSK_NotGray(j,:)=mapper_out_NotGray + I_scaled_noise + Q_scaled_noise*1i;
end

%%%%The Demapper and BER calculation%%%
BER=get_BER(SNR_vector,noisy_signal_QPSK,num_bits,Data,'QPSK',bits_per_symbol,QPSK_positions); 
BER_NotGray=get_BER(SNR_vector,noisy_signal_QPSK_NotGray,num_bits,Data,'QPSK',bits_per_symbol,QPSK_positions_NotGray); 
%calculate the theoritical BER
BER_theoretical=0.5 * erfc(sqrt(Eb ./ (No_vector*(Eavg/bits_per_symbol))));
%plot Demapper output BER vs theoretical
plot_BER(SNR_vector,BER,BER_theoretical,'QPSK - practical BER Vs theoritical BER');
plot_BER(SNR_vector,BER_NotGray,BER_theoretical,'QPSK Not Gray - practical BER Vs theoritical BER');

figure;
semilogy(SNR_vector,BER);
 hold on;
semilogy(SNR_vector,BER_NotGray);
hold on;
semilogy(SNR_vector,BER_theoretical);
title('QPSK - practical BER with and without Gray coding Vs theoritical BER'); 
ylim([1e-4 1]);
xlabel('Eb/No');
ylabel('BER');
legend("practical BER","practical BER without Gray","Theoretical BER");
hold off;
end
        %---------------8-PSK----------------%
function PSK8(Data,num_bits,Eb,SNR_vector,PSK8_positions)
bits_per_symbol=3;
%mapping to 16_QAM
mapper_out = mapp_function(Data,num_bits,'8-PSK',PSK8_positions,bits_per_symbol);
%%%%The Channel%%%%
%Generate a unity variance, zero mean additive white Gaussian noise signal
%with the same size as transmitted signal.
I_noise=randn(size(mapper_out));
Q_noise=randn(size(mapper_out));
Eavg=sum(abs(mapper_out).^2)/(num_bits/bits_per_symbol);
%Get the Noisy signal for each SNR value
for j=1:length(SNR_vector)
     %calculate variance from SNR (SNR=Es/Ns)
     No_vector(j)=Eb / 10^(SNR_vector(j)/10);
     variance_vector(j)=(No_vector(j)*(Eavg/bits_per_symbol))/2; %variance=No/2
     %Scale the noise sequence to have variance = N0/2 by multiplying the sequence
     %by sqrt(N0/2). 
     I_scaled_noise=sqrt(variance_vector(j)) * I_noise; 
     Q_scaled_noise=sqrt(variance_vector(j)) * Q_noise; 
     %Add the noise to the transmitted sequence
     noisy_signal_8PSK(j,:)=mapper_out + I_scaled_noise + Q_scaled_noise*1i;
end

%%%The Demapper and BER calculation%%%
BER=get_BER(SNR_vector,noisy_signal_8PSK,num_bits,Data,'8-PSK',bits_per_symbol,PSK8_positions); 

%calculate the theoritical BER
BER_theoretical=(1/3) * erfc(sin(pi/8)*sqrt(Eb./(No_vector*(Eavg/bits_per_symbol))));
%plot Demapper output BER vs theoretical
plot_BER(SNR_vector,BER,BER_theoretical,'8-PSK - practical BER Vs theoritical BER');

end
        %---------------16_QAM---------------%
function QAM(Data,num_bits,Eb,SNR_vector,QAM16_positions)
bits_per_symbol=4;
%mapping to 16_QAM
mapper_out = mapp_function(Data,num_bits,'16_QAM',QAM16_positions,bits_per_symbol);
%%%%The Channel%%%%
%Generate a unity variance, zero mean additive white Gaussian noise signal
%with the same size as transmitted signal.
I_noise=randn(size(mapper_out));
Q_noise=randn(size(mapper_out));
Eavg=sum(abs(mapper_out).^2)/(num_bits/bits_per_symbol);
%Get the Noisy signal for each SNR value
for j=1:length(SNR_vector)
     %calculate variance from SNR (SNR=Es/Ns)
     No_vector(j)=Eb / 10^(SNR_vector(j)/10);
     variance_vector(j)=(No_vector(j)*(Eavg/bits_per_symbol))/2; %variance=No/2
     %Scale the noise sequence to have variance = N0/2 by multiplying the sequence
     %by sqrt(N0/2). 
     I_scaled_noise=sqrt(variance_vector(j)) * I_noise; 
     Q_scaled_noise=sqrt(variance_vector(j)) * Q_noise; 
     %Add the noise to the transmitted sequence
     noisy_signal_16QAM(j,:)=mapper_out + I_scaled_noise + Q_scaled_noise*1i;
end

%%%The Demapper and BER calculation%%%
BER=get_BER(SNR_vector,noisy_signal_16QAM,num_bits,Data,'16_QAM',bits_per_symbol,QAM16_positions); 

%calculate the theoritical BER
BER_theoretical=(3/8) * erfc(sqrt(Eb ./ (No_vector*(Eavg/bits_per_symbol))));
%plot Demapper output BER vs theoretical
plot_BER(SNR_vector,BER,BER_theoretical,'16-QAM - practical BER Vs theoritical BER');
end
%___________________________________________________________________________________%
function mapper_out=mapp_function(Data,num_bits,mod_type,positions,bits_per_symbol)
 if strcmp(mod_type,'BPSK')
     mapper_out = (2*Data) - 1; 
 elseif strcmp(mod_type,'QPSK')|| strcmp(mod_type,'16_QAM') || strcmp(mod_type,'8-PSK')
    num_symbols = num_bits/bits_per_symbol;
    mapper_vector_bin=reshape(Data,bits_per_symbol,num_bits/bits_per_symbol)';
    for s=1:num_symbols
        mapper_vector_Dec(s) = bin2dec(num2str(mapper_vector_bin(s,:)));
        mapper_out(s)=positions(mapper_vector_Dec(s)+1);
    end
  
 end
end
function BER=get_BER(SNR_vector,recieved_signal,num_bits,Data,mod_type,bits_per_symbol,positions)
estimate_signal=complex(zeros(1,(num_bits/bits_per_symbol)));
phases=rad2deg(phase(recieved_signal(1,:)));
if strcmp(mod_type,'BPSK')
 for i=1:length(SNR_vector)
     %Estimation of symbols and Demapping to binary bits
     for c=1:num_bits
         if recieved_signal(i,c)>=0
            recieved_signal(i,c)=1;
         else
            recieved_signal(i,c)=0;
         end
     end
     %Calculate the bit error rate for each SNR value
     error_counter=0; 
     for c=1:num_bits
         if recieved_signal(i,c)~=Data(c)
            error_counter=error_counter+1;
         end
     end
     BER(i)=error_counter/num_bits;
 end
elseif strcmp(mod_type,'QPSK') || strcmp(mod_type,'16_QAM') || strcmp(mod_type,'8-PSK')
     for i=1:length(SNR_vector)
         %Estimation of symbols
         if strcmp(mod_type,'QPSK')
             for s=1:num_bits/bits_per_symbol
                 if real(recieved_signal(i,s))>=0 && imag(recieved_signal(i,s))>=0
                     estimate_signal(s)=1+1i;
                 elseif real(recieved_signal(i,s))>=0 && imag(recieved_signal(i,s))<0
                     estimate_signal(s)=1-1i;
                 elseif real(recieved_signal(i,s))<0 && imag(recieved_signal(i,s))>=0
                     estimate_signal(s)=-1+1i;
                 else
                     estimate_signal(s)=-1-1i;
                 end
             end
         elseif strcmp(mod_type,'16_QAM')
             for s=1:num_bits/bits_per_symbol
                  if real(recieved_signal(i,s))>0 
                      if real(recieved_signal(i,s))>2
                          estimate_signal(s)=3+imag(recieved_signal(i,s))*1i;
                      else
                          estimate_signal(s)=1+imag(recieved_signal(i,s))*1i;
                      end
                  else
                      if real(recieved_signal(i,s))<-2
                          estimate_signal(s)=-3+imag(recieved_signal(i,s))*1i;
                      else
                          estimate_signal(s)=-1+imag(recieved_signal(i,s))*1i;
                      end
                  end
                  if imag(estimate_signal(s))>0 
                      if imag(estimate_signal(s))>2
                          estimate_signal(s)=real(estimate_signal(s))+3i;
                      else
                          estimate_signal(s)=real(estimate_signal(s))+1i;
                      end
                  else
                      if imag(estimate_signal(s))<-2
                          estimate_signal(s)=real(estimate_signal(s))-3i;
                      else
                          estimate_signal(s)=real(estimate_signal(s))-1i;
                      end
                  end
             end
         elseif strcmp(mod_type,'8-PSK')
             for s=1:num_bits/bits_per_symbol
                      if rad2deg(phase(recieved_signal(i,s)))+360>337.5 && rad2deg(phase(recieved_signal(i,s)))+360<360
                          estimate_signal(s)=1+0i;
                      elseif rad2deg(phase(recieved_signal(i,s)))>0 && rad2deg(phase(recieved_signal(i,s)))<22.5
                          estimate_signal(s)=1+0i;
                      elseif rad2deg(phase(recieved_signal(i,s)))>22.5 && rad2deg(phase(recieved_signal(i,s)))<67.5
                          estimate_signal(s)=sqrt(1/2)+sqrt(1/2)*1i;
                      elseif rad2deg(phase(recieved_signal(i,s)))>67.5 && rad2deg(phase(recieved_signal(i,s)))<112.5
                          estimate_signal(s)=0+1i;
                      elseif rad2deg(phase(recieved_signal(i,s)))>112.5 && rad2deg(phase(recieved_signal(i,s)))<157.5
                          estimate_signal(s)=-sqrt(1/2)+sqrt(1/2)*1i;
                      elseif rad2deg(phase(recieved_signal(i,s)))>157.5 && rad2deg(phase(recieved_signal(i,s)))<180
                          estimate_signal(s)=-1+0i;
                      elseif rad2deg(phase(recieved_signal(i,s)))+360>180 && rad2deg(phase(recieved_signal(i,s)))+360<202.5
                          estimate_signal(s)=-1+0i;
                      elseif rad2deg(phase(recieved_signal(i,s)))+360>202.5 && rad2deg(phase(recieved_signal(i,s)))+360<247.5
                          estimate_signal(s)=-sqrt(1/2)-sqrt(1/2)*1i;
                      elseif rad2deg(phase(recieved_signal(i,s)))+360>247.5 && rad2deg(phase(recieved_signal(i,s)))+360<292.5
                          estimate_signal(s)=0-1i;
                      elseif rad2deg(phase(recieved_signal(i,s)))+360>292.5 && rad2deg(phase(recieved_signal(i,s)))+360<337.5
                          estimate_signal(s)=sqrt(1/2)-sqrt(1/2)*1i;
                      end
             end
         end
         %Demapping estimated symbols to the corresponding decimal values then to
         %binary bits
         for s=1:num_bits/bits_per_symbol
             for p=1:length(positions)
                 if estimate_signal(s) == positions(p)
                     Demapper_vector_Dec(s)=p-1;
                 end
             end
         end
         Demapper_vector_str=dec2bin(Demapper_vector_Dec);
         for b=1:bits_per_symbol
             Demapper_vector_bin(:,b)=str2num(Demapper_vector_str(:,b));
         end
         Demapper_out=reshape(Demapper_vector_bin',1,num_bits);
         %Calculate the bit error rate for each SNR value
         error_counter=0; 
         for c=1:num_bits
             if Demapper_out(c)~=Data(c)
               error_counter=error_counter+1;
             end
         end
         BER(i)=error_counter/num_bits;
         
     end
end
end
function plot_BER(SNR_vector,BER,BER_theoretical,str)
    figure;
    semilogy(SNR_vector,BER);
    hold on;
    semilogy(SNR_vector,BER_theoretical);
    ylim([1e-4 1]);
    title(str); 
    xlabel('Eb/No');
    ylabel('BER');
    legend("practical BER","Theoretical BER");
    hold off;
end
