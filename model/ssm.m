clear all
close all
clc

[s, Fs, Bt, opts] = wavread('speech.wav');
M = length(s);
d = wavread('noise.wav');
Z = length(d);

if Z<M then
    s=wavread('speech.wav',Z);
    M=Z;
else 
    d=wavread('noise.wav',M);
end

x=s+d;
wavwrite(x, Fs, Bt, 'inputMAT.wav');

N =256;        
h =hann(N); 
N2=N/2;       
L =floor(M/N2)-1;                  
C =1.0;        
P =4;          

Npow=zeros(N,1);    
for k=1:P
  n =N2*(k-1);
  x1=h.*x( n+1 : n+N );  
  X1=fft(x1);      
  Npow=Npow+abs(X1)/P; 
end
y=x(1:P*N2);    

past_tail=zeros(N2,1);  
for k=P+1:L
  n =N2*(k-1);
  x1=h.*x( n+1 : n+N );  
  X1=fft(x1);               
  G =ones(N,1)-C*Npow./abs(X1);  
  G =( G + abs(G) )/2;             
  Y1=G.*X1;                     
  y1=ifft(Y1);                
  y =[y; past_tail + y1(1:N2) ]; 
                                
  past_tail=y1(N2+1:N);         
end
wavwrite(y,Fs,Bt,'outputMAT.wav');