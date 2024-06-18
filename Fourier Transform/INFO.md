# Fourier Transform and Quantum Computing
## What is Fourier Transform?
### Fourier Series
Functions can be described as a summation of sinusoidal functions. More specifically, different frequencies of sines and cosines all with different amplitudes can estimate a given function over a given period, T. 
$$
s_n(x) = A_0 + \sum_{n=1}^N(A_n cos(2\pi\frac{n}{T}x) + B_n sin(2\pi\frac{n}{T}x))
$$
Cosine and sine in their polar form are $\frac{e^{i\theta}+e^{-i\theta}}{2}$ and $\frac{e^{i\theta}-e^{-i\theta}}{2i}$ respectively. Replacing the summation above with these new definition, we get a new summation where we will be using the rest of the derivation.
$$
s_n(x) = \sum_{n=-N}^NC_ne^{i2\pi\frac{n}{T}x}
$$
### Sampling
Most audio is captured at around 44.1kHz or 48kHz. This means every second, we measure the amplitude of audio 44100 times. The wav file contains 44100 points of data for one second, which we will analyze in this notebook. The csv file contains the time the amplitude it was recorded and the amplitude measured in the format: "time", "amplitude".
```code
0.00000,-0.46933
0.00002,-0.46011
0.00005,-0.44931
0.00007,-0.41455
0.00009,-0.38632
0.00011,-0.34164
...
0.99995,0.12454
0.99998,0.13571
```
To find the frequencies that make up the audio, we'll use the *discrete Fourier Transform* where $a_0=-0.46933$, $a_1=-0.46011$, ... and so on. Since there are 44100 points of data, we'll set $N=44100$. Each frequency, $\phi$ will be defined as:
$$
\phi_{k} = \frac{1}{\sqrt{N}}\sum_{j=0}^{N-1}a_je^{2\pi ijk/N}
$$
For example:
$$
\phi_{0} = \frac{1}{44100}(-0.46933e^{2\pi i(0)(0)/44100}-0.46011e^{2\pi i(1)(0)/44100}+...+0.13571e^(44099)(0)(44100))=-0.0973861\\
\phi_{1} = \frac{1}{44100}(-0.46933e^{2\pi i(0)(1)/44100}-0.46011e^{2\pi i(1)(1)/44100}+...+0.13571e^(44099)(1)(44100))=−0.118737+0.136405i\\
. \\
. \\
. \\
\phi_{44099} = \frac{1}{44100}(-0.46933e^{2\pi i(0)(44099)/44100}-0.46011e^{2\pi i(1)(44099)/44100}+...+0.13571e^(44099)(44099)(44100))=−0.118737−0.136405i
$$
One thing to notice is that $\phi_{1}=\phi_{44099}^*$. This remain true such that $\phi_{k}=\phi_{N-k}^*$ for k = $1, 2, ... , N/2 - 1$ so $\phi_{0}$ and $\phi{N/2}$ are unique. Taking the norm and plotting the frequencies against amplitude results in a frequency spectrum. 
### Classical Computing
Calculating one $\phi_{k}$ requires summing N terms with N different $\phi_{k}$, this computation becomes $O(N^2)$. Although this is polynomial time, which is better than some other problems out there such as Traveling Salseman, when analyzing large amounts of data this algorithm can become very slow. Another way to approach this problem is to think of it as a matrix-vector multiplication where $\omega = e^{2\pi i/N}$. The fourier transform becomes:
$$
\phi_{k} = \frac{1}{\sqrt{N}}\sum_{j=0}^{N-1}a_j\omega^{jk}
$$
For example:
$$
\phi_{0} = \frac{1}{\sqrt{N}}(a_0 + a_1 + a_2 + ... + a_{N-1}),\\
\phi_{1} = \frac{1}{\sqrt{N}}(a_0 + a_1\omega +  a_2\omega^2 + ... + a_{N-1}\omega^{N-1}),\\
\phi_{2} = \frac{1}{\sqrt{N}}(a_0 + a_1\omega^2 +  a_2\omega^4 + ... + a_{N-1}\omega^{2(N-1)}),\\
. \\
. \\
. \\
\phi{k} = \frac{1}{\sqrt{N}}(a_0 + a_1\omega^{N-1} +  a_2\omega^{2(N-1)} + ... + a_{N-1}\omega^{(N-1)^2})
$$
Writing the equations as a matrix results in:
$$
\begin{pmatrix}
    \phi_{0}\\
    \phi_{1}\\
    \phi_{2}\\
    \vdots\\
    \phi_{N-1}\\
\end{pmatrix} = 
\begin{pmatrix}
    1 & 1 & 1 & \cdot & 1 \\
    1 & \omega & \omega^2 & \cdots & \omega^{N-1} \\
    1 & \omega^2 & \omega^4 & \cdots & \omega^{2(N-1)} \\
    \vdots & \vdots & \vdots & \ddots & \vdots \\
    1 & \omega^{N-1} & \omega^{2(N-1)} & \cdots & \omega{(N-1)^2}
\end{pmatrix}
\begin{pmatrix}
    a_0 \\
    a_1 \\
    a_2 \\
    \vdots\\
    a_{N-1}
\end{pmatrix}
$$
Although this is useful conceptually, the amount of operations are still the same, $O(N^2)$. There are some algorithms out there that have runtimes of $O(N\log N)$ but biggest thing to note is that the $N$ x $N$ matrix is unitary.

### Quantum Computing
