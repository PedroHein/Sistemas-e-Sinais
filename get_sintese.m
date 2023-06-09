function [impulso_repetido] = get_sintese(x, y, ajuste)
    % Calculo a frequencia e as formantes da vogal            
    % Vetor tempo
    x = x(:, 1);
    N = length(x);          
    T = inv(y);                   
    D = (N - 1) * T;            
    tempo = linspace(0, D, N);     
    
    % Gráfico sinal 'puro'
    figure()
    plot(tempo, x)
    title("Sinal da vogal 'pura'")
    xlabel('Frequência')
    ylabel('Amplitude')

    % Janelamento
    janela = hamming(N);
    vogal = x .* janela;
    
    % Frequencia
    fft_vog = abs(fftshift(fft(vogal)));
    freq = linspace(-y/2, y/2, length(vogal));
    [~, index] = max(fft_vog);
    freq_fund = abs(freq(index));       
    
    % Fator de segurança
    fator_fft = freq_fund*0.9;
    
    % Base logaritmica
    fft_log = log10(fft_vog);     
    
    % Gráfico
    figure()
    plot(abs(freq), fft_log)

    % Valores positivos
    indice_positivo = freq > 0;
    freq_positivo = freq(indice_positivo);
    fft_log_pos = fft_log(indice_positivo);

    % Macropicos
    [picos_x,picos_y] = findpeaks(fft_log_pos, freq_positivo,'MinPeakDistance',fator_fft);
    title('Espectro de amplitude da vogal ');
    xlabel('Frequência em Hz');
    ylabel('Magnitude em dB');

    % Gráfico do sinal e a derivada sobre os picos
    figure()
    plot(picos_y,picos_x)
    hold on
    findpeaks(fft_log_pos, freq_positivo,'MinPeakDistance',fator_fft);
    hold off
    title("Função descontinua da vogal")

    [~, LOCSdentro] = findpeaks(abs(picos_x), picos_y, 'MinPeakDistance', 10);

    n = find(LOCSdentro>60);
    
    FREQUENCIASfn = LOCSdentro(n);
    
    % Resultados
    f0 = FREQUENCIASfn(1)       % Frequencia fundamental
    f1 = FREQUENCIASfn(2)       % Formante 1
    f2 = FREQUENCIASfn(3)       % Formante 2

    % Sintese da vogal para a reprodução
    % Picos em um determinado trecho
    [~,ip] = findpeaks(diff(x),tempo(1:end-1),'MinPeakDistance',0.001,'MinPeakDistance',0.005);
    inicio = find(tempo==ip(31));
    fim = find(tempo==ip(32));

    X = tempo(inicio:fim) - tempo(inicio);
    Y = x(inicio:fim);
    
    [vexp,texp] = findpeaks(Y,X,'MinPeakDistance',0.001);
    f = fit(texp',vexp,'exp1')
    
    % Gráfico da curva
    figure()
    plot(f,texp',vexp)
    
    % Transformando a função temporal em uma função de Laplace
    l =[0 f.a f.b ajuste]';
    Y_new = l(1) + l(2).*exp(l(3).*X).*sin(2*pi*l(4).*X);
    
    % Gráfico sinal original e sintetizado
    figure()
    plot(X,Y,'+r',X,Y_new,'b')
    title('Sinal de original X Sinal ajustado');
    xlabel('Tempo em s');
    ylabel('Tensão em Volts');  
    legend('Original','Ajustado')
    
    %Continuação de Laplace
    k = l(2)*2*pi*l(4);
    p1 = l(3) + j*2*pi*l(4);
    p2 = l(3) - j*2*pi*l(4);
    z = [];
    Gs = zpk(z,[p1 p2],k)
    
    % Gráfico do impulso
    figure()
    impulse(Gs)
    
    % Variaveis para reprodução da vogal
    repeticoes =    round(freq_fund*D);
    [impulso, t] = impulse(Gs);
    
    impulso_repetido = repmat(impulso, 1, repeticoes);
    impulso_repetido = impulso_repetido(:)';
end




