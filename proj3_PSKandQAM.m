clear all;
close all;

num_bits=12000;
Eb=1;
SNR_vector= -4:14;
%Generate random bits 
Data = randi([0 1] , 1 , num_bits);
%Define positions on signal space
BPSK_positions=[1 -1];
QPSK_positions=[-1-1i -1+1i 1-1i 1+1i];
QAM16_positions=[-3-3i -3-1i 3+3i -3+1i -1-3i -1-1i -1+3i -1+1i 3-3i 3-1i 3+3i 3+1i 1-3i 1-1i 1+3i 1+1i];
%---------------BPSK----------------%
bits_per_symbol=1;
%%%%the mapper%%%%
%mapping to BPSK (convert logic 0 to -1)
mapper_out = mapp_function(Data,num_bits,'BPSK',BPSK_positions,bits_per_symbol);

%%%%The Channel%%%%
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
 
%---------------QPSK----------------%
bits_per_symbol=2;
%mapping to QPSK
mapper_out = mapp_function(Data,num_bits,'QPSK',QPSK_positions,bits_per_symbol);
%%%%The Channel%%%%
%Generate a unity variance, zero mean additive white Gaussian noise signal
%with the same size as transmitted signal.
I_noise=randn(size(mapper_out));
Q_noise=randn(size(mapper_out));
%Get the Noisy signal for each SNR value
for j=1:length(SNR_vector)
     %calculate variance from SNR (SNR=Es/Ns)
     Eavg=sum(abs(mapper_out).^2)/(num_bits/bits_per_symbol);
     No_vector(j)=Eb / 10^(SNR_vector(j)/10);
     variance_vector(j)=(No_vector(j)*(Eavg/bits_per_symbol))/2; %variance=Ns/2
     %Scale the noise sequence to have variance = N0/2 by multiplying the sequence
     %by sqrt(N0/2). 
     I_scaled_noise=sqrt(variance_vector(j)) * I_noise; 
     Q_scaled_noise=sqrt(variance_vector(j)) * Q_noise; 
     %Add the noise to the transmitted sequence
     noisy_signal_QPSK(j,:)=mapper_out + I_scaled_noise + Q_scaled_noise*1i;
end

%%%%The Demapper and BER calculation%%%
BER=get_BER(SNR_vector,noisy_signal_QPSK,num_bits,Data,'QPSK',bits_per_symbol,QPSK_positions); 

%calculate the theoritical BER
BER_theoretical=0.5 * erfc(sqrt(Eb ./ No_vector));
%plot Demapper output BER vs theoretical
plot_BER(SNR_vector,BER,BER_theoretical,'QPSK - practical BER Vs theoritical BER');

%---------------16_QAM----------------%

bits_per_symbol=4;
%mapping to QPSK
mapper_out = mapp_function(Data,num_bits,'16_QAM',QAM16_positions,bits_per_symbol);
%%%%The Channel%%%%
%Generate a unity variance, zero mean additive white Gaussian noise signal
%with the same size as transmitted signal.
I_noise=randn(size(mapper_out));
Q_noise=randn(size(mapper_out));
%Get the Noisy signal for each SNR value
for j=1:length(SNR_vector)
     %calculate variance from SNR (SNR=Es/Ns)
     Eavg=sum(abs(mapper_out).^2)/(num_bits/bits_per_symbol);
     No_vector(j)=Eb / 10^(SNR_vector(j)/10);
     variance_vector(j)=(No_vector(j)*(Eavg/bits_per_symbol))/2; %variance=No/2
     %Scale the noise sequence to have variance = N0/2 by multiplying the sequence
     %by sqrt(N0/2). 
     I_scaled_noise=sqrt(variance_vector(j)) * I_noise; 
     Q_scaled_noise=sqrt(variance_vector(j)) * Q_noise; 
     %Add the noise to the transmitted sequence
     noisy_signal_16QAM(j,:)=mapper_out + I_scaled_noise + Q_scaled_noise*1i;
end

%%%%The Demapper and BER calculation%%%
BER=get_BER(SNR_vector,noisy_signal_16QAM,num_bits,Data,'16_QAM',bits_per_symbol,QAM16_positions); 

%calculate the theoritical BER
BER_theoretical=0.5 * erfc(sqrt(Eb ./ No_vector));
%plot Demapper output BER vs theoretical
plot_BER(SNR_vector,BER,BER_theoretical,'QPSK - practical BER Vs theoritical BER');
%_____________________________________________________%
function mapper_out=mapp_function(Data,num_bits,mod_type,positions,bits_per_symbol)
 if strcmp(mod_type,'BPSK')
     mapper_out = (2*Data) - 1; 
 elseif strcmp(mod_type,'QPSK')|| strcmp(mod_type,'16_QAM')
    num_symbols = num_bits/bits_per_symbol;
    mapper_vector_bin=reshape(Data,bits_per_symbol,num_bits/bits_per_symbol)';
    for s=1:num_symbols
        mapper_vector_Dec(s) = bin2dec(num2str(mapper_vector_bin(s,:)));
        mapper_out(s)=positions(mapper_vector_Dec(s)+1);
    end
  
 end

end
function BER=get_BER(SNR_vector,noisy_signal,num_bits,Data,mod_type,bits_per_symbol,positions)
if strcmp(mod_type,'BPSK')
 for i=1:length(SNR_vector)
     %Demapper output
     for c=1:num_bits
         if noisy_signal(i,c)>=0
            noisy_signal(i,c)=1;
         else
            noisy_signal(i,c)=0;
         end
     end
     %Calculate the bit error rate for each SNR value
     error_counter=0; 
     for c=1:num_bits
         if noisy_signal(i,c)~=Data(c)
            error_counter=error_counter+1;
         end
     end
     BER(i)=error_counter/num_bits;
 end
elseif strcmp(mod_type,'QPSK')
     for i=1:length(SNR_vector)
         for s=1:num_bits/bits_per_symbol
             if real(noisy_signal(i,s))>=0 && imag(noisy_signal(i,s))>=0
                 estimate_signal(s)=1+1i;
             elseif real(noisy_signal(i,s))>=0 && imag(noisy_signal(i,s))<0
                 estimate_signal(s)=1-1i;
             elseif real(noisy_signal(i,s))<0 && imag(noisy_signal(i,s))>=0
                 estimate_signal(s)=-1+1i;
             else
                 estimate_signal(s)=-1-1i;
             end
         end
         %Demapper output
         for s=1:num_bits/bits_per_symbol
             for p=1:length(positions)
                 if estimate_signal(s) == positions(p)
                     symbols_vector_Dec(s)=p-1;
                 end
             end
         end
         symbols_vector_str=dec2bin(symbols_vector_Dec);
         for b=1:bits_per_symbol
             symbols_vector_bin(:,b)=str2num(symbols_vector_str(:,b));
         end
         Demapped_out=reshape(symbols_vector_bin',1,num_bits);
         %Calculate the bit error rate for each SNR value
         error_counter=0; 
         for c=1:num_bits
             if Demapped_out(c)~=Data(c)
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
    title(str); 
    xlabel('Eb/No');
    ylabel('BER');
    legend("practical BER","Theoretical BER");
    hold off;
end