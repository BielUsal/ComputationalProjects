# Evaluación Contínua FísComp 1 - Teorema de Wiener-Khinchin

#### Por Gabriel D'Andrade Furlanetto XDD204950 y Álvaro Gamarra Ralero

## Análisis de la función autocorrelación de $t=0ps$ hasta $t=5ps$

Tenemos qu en este intervalo de tiempo, las velocidades están distribuidas segundo:

![P0](velol.png)

Y que la autocorrelación es:

![P1](autocor.png)

## Cálculo de la PSD usando el teorema de Wiener-Kinchin
Tomamos la transformada de Fourier de la función de autocorrelación y obtenemos que la PSD es:

![P2](PSD.png)

*En virtud del teorema de muestreo de Nyquist‐Shannon ¿hasta qué frecuencia, Bmax, en GHz podría calcular la densidad espectral SXX(f)?*

El teorema de mostreo de Shannon nos dice que la frecuencia máxima que podemos obtener con una transformada de Fourier es:

$$B_{max} = \frac{f_s}{2} = \frac{1}{2t_s}$$

De manera que, como tenemos un mostreo de $t_s = 10 fs$, se puede obtener una frecuencia máxima de 

$$  B_{max} = \frac{1}{20 fs} = 0.05 10^{15} Hz= 5 10^{4} GHz = 50 THz $$

## Cálculo de la PSD directamente de la Transformada de Fourier

### Transformada de Fourier "directa"