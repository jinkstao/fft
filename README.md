# Fast Fourier Transform

---

This is a fast fourier transform algorithm.

- ## Usage

1. Call **fft\_init(int sampling\_point\_count, void \*\*fft\_id\_ptr)** to initialize the FFT. The parameter **sampling\_point\_count** is the length of the time-domain sequence involved in the FFT, whose value must be an integer power of 2. And the parameter **fft\_id\_ptr** is used to obtain the identity of the FFT, which is a two-level pointer. You should check the return value of the function, which equals to 0 means the initialization failed.

2. When you need to do FFT or IFFT, you can call **fft\_d2fft\_real(const void \*fft\_id, const float \*input\_seq, float \*output\_seq, int mode)** or **fft\_d2fft\_complex(const void \*fft\_id, const fft\_complex \*input\_seq, fft\_complex \*output\_seq, int mode)**. The former is a real number domain FFT, while the latter is a complex number domain FFT. The parameter **fft\_id** is the first-level pointer to the second-level pointer **fft\_id\_ptr** you obtained from **fft\_init**. The parameter **input\_seq/output\_seq** means input/output sequence, whose length equals to **int sampling_point_count** in **fft\_init**. And parameter **mode(=FFT or IFFT)** determines whether the function performs FFT or IFFT. You can call this function any number of times after **fft\_init**. This function has no return value.

3. When you finish your work, remember to call **fft\_close(void \*fft\_id)** to free memory. The parameter **fft\_id** is same to the one in **fft\_d2fft\_real** and **fft\_d2fft\_complex**. This function has no return value.