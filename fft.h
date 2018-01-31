#ifndef _FFT_H_
#define _FFT_H_

#ifdef WIN32
#ifdef DLL_EXPORT
#define API_DECLSPEC  __declspec(dllexport)
#else
#define API_DECLSPEC  __declspec(dllimport)
#endif  // DLL_EXPORT
#endif  // WIN32

#ifdef unix
#define API_DECLSPEC
#endif  // UNIX

#ifdef _cplusplus
extern "C"
{
#endif

static const int FFT = 0;
static const int IFFT = 1;

typedef struct _fft_complex
{
    float real;
    float imag;
}fft_complex;

typedef struct _fft_info
{
    int sampling_point_count;
    int sampling_point_level;
    fft_complex *polar_sequence;
    fft_complex *tmp_sequence;
    fft_complex *tmp_output_sequence;
}fft_info;

API_DECLSPEC int fft_init(int sampling_point_count, void **fft_id_ptr);

API_DECLSPEC void fft_close(void *fft_id);

API_DECLSPEC void fft_d2fft_real(const void *fft_id,
                                 const float *input_seq,
                                 float *output_seq,
                                 int mode);

API_DECLSPEC void fft_d2fft_complex(const void *fft_id,
                                    const fft_complex *input_seq,
                                    fft_complex *output_seq,
                                    int mode);

#ifdef _cplusplus
}
#endif

#endif  // _FFT_H_