function [] = get_freq_e_form(VG, fs)
    %Calcula a frequencia e as formantes da vogal
    % Vetor tempo
    VG = VG(:, 1);
    Na = length(VG);          
    Ta = inv(fs);                   % taxa de amostragem 
    Da = (Na - 1) * Ta;             % Duração do sinal
    tempo = linspace(0, Da, Na);     
    figure(1)
    plot(tempo, VG)
    title("Sinal da vogal 'pura'")
    xlabel('Frequência')
    ylabel('Amplitude')

    % Janelamento do sinal
    janela = hamming(Na);
    vogal = VG .* janela;
    
    fft_vog = abs(fftshift(fft(vogal)));
    freq = linspace(-fs/2, fs/2, length(vogal));
    [~, index] = max(fft_vog);
    freq_fund = abs(freq(index));       % Frequencia fundamental
    
    fator_fft = freq_fund*0.9;          
    disp(freq_fund);
    
    fft_log = log10(fft_vog);     
    figure(2)
    plot(freq, fft_log)

    % Isolando apenas os valores positivos de frequência e seus respectivos
    % valores de energia

    indice_positivo = freq > 0;
    freq_positivo = freq(indice_positivo);
    fft_log_pos = fft_log(indice_positivo);

    % obtendo os macropicos
    [picos_x,picos_y] = findpeaks(fft_log_pos, freq_positivo,'MinPeakDistance',fator_fft);
    figure(3)
    plot(picos_y,picos_x)
    title("Função descontinua da vogal")

    % Utilização do comando interp1

    resolucao_atual = length(picos_x);
    resolucao_desejada = 5 * resolucao_atual;

    indices_atual = 1:resolucao_atual;
    indices_desejados = linspace(1, resolucao_atual, resolucao_desejada);

    peaks_y_resampled = interp1(indices_atual, picos_y, indices_desejados, 'spline');
    peaks_x_resampled = interp1(indices_atual, picos_x, indices_desejados, 'spline');

    % Ajustando uma curva suave aos pontos originais
    ajuste_curva = fit(peaks_y_resampled', peaks_x_resampled', 'smoothingspline');

    % Definindo a resolução desejada para a função contínua
    resolucao_desejada = 600;

    % Criando um vetor de valores suavemente espaçados para a variável independente
    x_continuo = linspace(min(peaks_y_resampled), max(peaks_y_resampled), resolucao_desejada);

    % Avaliando a curva ajustada nos valores contínuos de x
    y_continuo = feval(ajuste_curva, x_continuo);

    % Plot da função contínua
    figure(4)
    plot(x_continuo, y_continuo)
    title("Função continua da vogal")

    derivada = diff(y_continuo);
    segunda_derivada = diff(derivada);
    
    x_continuo2 = x_continuo(1:end-2);

    % Indices de maximo e minimo
    indices_mudanca_sinal = (find(diff(sign(derivada)) ~= 0)+1);

    % Verificando se a segunda derivada é negativa nos pontos de mudança de sinal
    pontos_segunda_derivada_negativa = x_continuo2(indices_mudanca_sinal(2:end)); % Ignorando o primeiro ponto (frequencia fundamental)

    % Formantes (em HZ)
    formants = pontos_segunda_derivada_negativa(segunda_derivada(indices_mudanca_sinal(2:end)) < 0)
end