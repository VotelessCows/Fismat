
 
clear all
close all
%-----------------------------
% Generación y guardado de muestras de audio
%-----------------------------
Fs = 44100; % Frecuencia de sampleo típica para una señal de audio.
nBits = 16; % Se trabaja con una resolución de 16 bits.
nChannel = 1; % Se grabará solo un canal.
duracion = 2; % Tiempo total de grabación.
load flag.mat flag % Se carga la variable flag de base de datos. Si es 1 entonces el sistema se encuentra encendido. Si es 0 se encuentra apagado.
error = zeros(1,2); % Vector que contendrá los errores provenientes de la diferencia entre las DEE.
%-----------------------------
% Grabado de la muestra y presentación de instrucciones. 
%-----------------------------
sample = audiorecorder(Fs,nBits,nChannel);
disp('Instrucciones: se procederá a grabar una muestra de la canción. Los comandos disponibles para interactuar con el sistema son "prender" y "apagar".')
%-----------------------------  
if(flag)
    disp('El sistema se encuentra: ENCENDIDO.');
else                                                    % Se le informa al usuario el estado del sistema que fue cargado desde la base de datos.
    disp('El sistema se encuentra: APAGADO.');
end
%------------------------------
disp('A continuación se grabará el comando que desee ejecutar.')
disp('Presione una tecla para continuar.')
pause;
disp('Grabando...')
recordblocking(sample, duracion);
disp('Grabación finalizada.')
disp('Procesando...')
recorded_sample = getaudiodata(sample);
audiowrite('muestra.wav', recorded_sample, Fs); % Se guarda la muestra en el directorio de trabajo.

[y,Fs] = audioread('muestra.wav');
y = y/max(abs(y)); % Se normaliza la amplitud (volumen) de la señal en el tiempo.
T = 1/Fs; % Tiempo representativo de cada muestra.
N = length(y); % Total de muestras.
t = (0:N-1)*T; % Vector de tiempo.
NFFT = 2^nextpow2(N); % Optimización del algoritmo FFT a través de un vector de muestras que sea una potencia de 2.
Y = fft(y, NFFT); % Se calcula la transformada de Fourier.
Y = Y(1:NFFT/2); % Como se va a trabajar con señales reales, se recorta el espectro y se trabaja con el espectro unilateral.
Yabs = abs(Y); % Se calcula el módulo de la transformada.
DEEs = Yabs.^2; % Densidad espectral de energía de la muestra.
f = (0:NFFT/2-1)*Fs/NFFT; % Vector de frecuencia en Hz.
%------------------------------
[v,Fs] = audioread('BeatIt.wav');
v = v/max(abs(v)); % Normalización volumen
T1 = 1/Fs;
N1 = length(v);
t1 = (0:N1-1)*T1;
NFFT1 = 2^nextpow2(N1);
V = fft(v, NFFT);
V = V(1:NFFT/2);
Vabs = abs(V);
DEEa = Vabs.^2; 
f = (0:NFFT1/2-1)*Fs/NFFT1;
%------------------------------
[x,Fs] = audioread('Bohemian.wav');
x = x/max(abs(x));
X = fft(x, NFFT);
X = X(1:NFFT/2);
Xabs = abs(X);
DEEc = Xabs.^2;
%------------------------------
[b,Fs] = audioread('InDaClub.wav');
b = b/max(abs(b)); % Normalización volumen
B = fft(b, NFFT);
B = B(1:NFFT/2);
Babs = abs(B);
DEEb = Babs.^2; 
%------------------------------
[c,Fs] = audioread('WithoutMe.wav');
c = c/max(abs(c)); % Normalización volumen
C = fft(c, NFFT);
C = C(1:NFFT/2);
Cabs = abs(C);
DEEd = Cabs.^2; 
%------------------------------
[d,Fs] = audioread('DarkHorse.wav');
d = d/max(abs(d)); % Normalización volumen
D = fft(d, NFFT);
D = D(1:NFFT/2);
Dabs = abs(D);
DEEe = Dabs.^2; 
%------------------------------
%Condidión de decisión: mínimo valor medio entre la resta de Densidades.
%------------------------------
error(1) = mean(abs(DEEs(1:Fs/2)-DEEa(1:Fs/2))); % Se calculan los errores en base al módulo de la diferencia entre las densidades espectrales de energía.
error(2) = mean(abs(DEEs(1:Fs/2)-DEEc(1:Fs/2))); % Se tienen en cuenta los valores que llegan hasta la frecuencia de Nyquist (fmax = fs/2).
%------------------------------
%------------------------------
%Condición de decisión: máximo valor de la correlación cruzada de las FFT unilaterales
%------------------------------
[c1,lag1] = xcorr(Yabs,Vabs);
[c2,lag2] = xcorr(Yabs,Xabs);
c1_max = max(c1);
c2_max = max(c2);
%------------------------------
disp('Utilizando correlación cruzada del módulo del espectro de la FFT:')
switch (c1_max > c2_max)
    case 1
        disp('El comando mencionado fue: Beat It de Michael Jackson.')
    case 0
        disp('El comando mencionado fue: Bohemian Rhapsody de Queen.')
end
%------------------------------
disp('Utilizando el valor promedio del módulo de la resta de las DEE:')
switch min(error)
   case error(1)
       disp('El comando mencionado fue: Beat It de Michael Jackson.')
   case error(2)
       disp('El comando mencionado fue: Bohemian Rhapsody de Queen.')
end
%------------------------------
coincidencia = (min(error)== error(1)) && (c1_max > c2_max) || (min(error)== error(2)) && (c1_max < c2_max); % Si la misma canción es detectada por ambos métodos, entonces coincidencia = 1.

if (coincidencia)
    if(flag)
        if((min(error)== error(1)))
            disp('El sistema ya se encuentra en funcionamiento')
        else
            disp('Bohemian Rhapsody de Queen')
            flag = 0;
        end
    else
        if((min(error)== error(1)))
            disp('Beat It de Michael Jackson')
            flag = 1;
        else
            disp('Bohemian Rhapsody de Queen')
        end
    end
else
    disp('No se tomará acción por falta de coincidencia')
end            
%------------------------------
save flag.mat flag % Se guarda la variable flag modificada en base de datos.
%------------------------------
%Gráficos
%------------------------------
subplot(3,2,3);
plot(f(1:Fs/4),DEEa(1:Fs/4));
xlabel('Frecuencia [Hz]');
ylabel('|V(f)|²');
title('Densidad espectral de energía: Beat It');
grid on;
axis tight;
%------------------------------
subplot(3,2,4);
plot(f(1:Fs/4),DEEc(1:Fs/4));
xlabel('Frecuencia [Hz]');
ylabel('|X(f)|²');
title('Densidad espectral de energía: Bohemian Rhapsody');
grid on;
axis tight;
%------------------------------
subplot(3,2,5);
plot(lag1/Fs,c1);
xlabel('Tiempo [s]');
ylabel('Amplitud');
title('Correlación cruzada muestra-Beat It');
grid on;
axis tight;
%------------------------------
subplot(3,2,6);
plot(lag2/Fs,c2);
xlabel('Tiempo [s]');
ylabel('Amplitud');
title('Correlación cruzada muestra-Bohemian Rhapsody');
grid on;
axis tight;
%------------------------------
%------------------------------
figure;
plot(f(1:Fs/4),DEEb(1:Fs/4));
xlabel('Frecuencia [Hz]');
ylabel('|V(f)|²');
title('Densidad espectral de energía: In Da Club');
grid on;
axis tight;
%------------------------------
figure;
plot(f(1:Fs/4),DEEd(1:Fs/4));
xlabel('Frecuencia [Hz]');
ylabel('|V(f)|²');
title('Densidad espectral de energía: Without Me');
grid on;
axis tight;
%------------------------------
figure;
plot(f(1:Fs/4),DEEe(1:Fs/4));
xlabel('Frecuencia [Hz]');
ylabel('|V(f)|²');
title('Densidad espectral de energía: Dark Horse');
grid on;
axis tight;