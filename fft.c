#include <stdio.h>
#include <malloc.h>
#include <math.h>
#define DLL_EXPORT
#include "fft.h"

/******************************************MACRO DEFINE******************************************/

#define	    PI      3.14159265	    // PI	    

/************************************************************************************************/

/****************************************GLOBAL VARIABLES****************************************/

// static int fft_sampling_point_count;            // count of sampling points(2^n)
// static int fft_sampling_point_level;            // log2(fft_sampling_point_count)
// static fft_complex *fft_polar_sequence;         // this array store temp values of polar
// static fft_complex *fft_tmp_sequence;           // this array store temp complex values of input
// static fft_complex *fft_tmp_output_sequence;    // this array store temp complex values of output

/************************************************************************************************/


/**************************************Function Declaration**************************************/

static void _fft_complex_plus(const fft_complex *complex1,
                       const fft_complex *complex2,
                       fft_complex *result);

static void _fft_complex_minus(const fft_complex *complex1,
                        const fft_complex *complex2,
                        fft_complex *result);

static void _fft_complex_multiple(const fft_complex *complex1,
                           const fft_complex *complex2,
                           fft_complex *result,
                           int mode);

static float _fft_complex_module(const fft_complex *complex);

/************************************************************************************************/


/*************************************Function Implementation*************************************/

static inline void _fft_complex_plus(const fft_complex *complex1,
                                     const fft_complex *complex2,
                                     fft_complex *result)
{
    result->real = complex1->real + complex2->real;
	result->imag = complex1->imag + complex2->imag;
}

static inline void _fft_complex_minus(const fft_complex *complex1,
                                      const fft_complex *complex2,
                                      fft_complex *result)
{
    result->real = complex1->real - complex2->real;
	result->imag = complex1->imag - complex2->imag;
}

static inline void _fft_complex_multiple(const fft_complex *complex1,
                                         const fft_complex *complex2,
                                         fft_complex *result,
                                         int mode)
{
    float complex2imag = mode == IFFT ? -complex2->imag : complex2->imag;
    result->real = complex1->real * complex2->real - complex1->imag * complex2imag;
	result->imag = complex1->real * complex2imag + complex1->imag * complex2->real;
}

static inline float _fft_complex_module(const fft_complex *complex)
{
    return sqrt(complex->real * complex->real + complex->imag * complex->imag);
}

int fft_init(int sampling_point_count, void **fft_id_ptr)
{
    int i;
    float angle;
	float fValue;
    *fft_id_ptr = (fft_info*)malloc(sizeof(fft_info));
    fft_info *info = (fft_info*)(*fft_id_ptr);
    // fft_sampling_point_count = sampling_point_count;
    info->sampling_point_count = sampling_point_count;
    int tmp = 1;
    info->sampling_point_level = 0;
    while(tmp < info->sampling_point_count)
    {
        tmp *= 2;
        ++info->sampling_point_level;
    }
    if(tmp != info->sampling_point_count)
    {
        // sampling point count must be 2^n!
        printf("FATAL ERROR: SAMPLING POINT COUNT MUST BE 2^N.\n");
        return 0;
    }
    // Initialize polar sequence and temp sequence
    // fft_polar_sequence = (fft_complex*)malloc(fft_sampling_point_count / 2 * sizeof(fft_complex));
    // fft_tmp_sequence = (fft_complex*)malloc(fft_sampling_point_count * sizeof(fft_complex));
    // fft_tmp_output_sequence = (fft_complex*)malloc(fft_sampling_point_count * sizeof(fft_complex));
    info->polar_sequence = (fft_complex*)malloc(info->sampling_point_count / 2 * sizeof(fft_complex));
    info->tmp_sequence = (fft_complex*)malloc(info->sampling_point_count * sizeof(fft_complex));
    info->tmp_output_sequence = (fft_complex*)malloc(info->sampling_point_count * sizeof(fft_complex));
    // Euler Formula
    for (i = 0; i < info->sampling_point_count / 2; i++)
    {
        angle = -2 * PI / info->sampling_point_count * i;
        // (fft_polar_sequence + i)->real = cos(angle);
        // (fft_polar_sequence + i)->imag = sin(angle);
        (info->polar_sequence + i)->real = cos(angle);
        (info->polar_sequence + i)->imag = sin(angle);
    }
    return 1;
}

void fft_close(void *fft_id)
{
    fft_info *info = (fft_info*)fft_id;
    free(info->polar_sequence);
    free(info->tmp_sequence);
    free(info->tmp_output_sequence);
    info->polar_sequence = 0;
    info->tmp_sequence = 0;
    info->tmp_output_sequence = 0;
    free(info);
    info = 0;
}

void fft_d2fft_real(const void *fft_id,
                    const float *input_seq,
                    float *output_seq,
                    int mode)
{
    fft_info *info = (fft_info*)fft_id;
    int i, j, k, interval, half_n, gap;
    unsigned int rev, num;
    fft_complex tmp_complex;

    // Copy real input array to temp complex array
    for (i = 0; i < info->sampling_point_count; ++i)
    {
        (info->tmp_sequence + i)->real = *(input_seq + i);
        (info->tmp_sequence + i)->imag = 0.0f;
    }

    // Butterfly Transform
    for (i = 0; i < info->sampling_point_level; i++)
    {
        interval = 1 << i;	// store count of groups
        half_n = 1 << (info->sampling_point_level - i);	// store length of every group
        for (j = 0; j < interval; j++)
        {
            // j means group j
            gap = j * half_n;
            for (k = 0; k < half_n / 2; k++)
            {
                _fft_complex_plus(info->tmp_sequence + k + gap,
                                  info->tmp_sequence + k + gap + half_n / 2,
                                  info->tmp_output_sequence + k + gap);
                _fft_complex_minus(info->tmp_sequence + k + gap,
                                   info->tmp_sequence + k + gap + half_n / 2,
                                   &tmp_complex);
                _fft_complex_multiple(&tmp_complex,
                                      info->polar_sequence + k * interval,
                                      info->tmp_output_sequence + k + gap + half_n / 2,
                                      mode);
            }
        }
        // copy result to input
        for (j = 0; j < info->sampling_point_count; j++)
        {
            *(info->tmp_sequence + j) = *(info->tmp_output_sequence + j);
        }
    }

    // invert the output code bit
    for (i = 0; i < info->sampling_point_count; i++)
    {
        rev = 0;
        num = (unsigned int)i;
        for (j = 0; j < info->sampling_point_level; j++)
        {
            rev <<= 1;
            rev |= num & 1;
            num >>= 1;
        }
        if(mode == IFFT)
        {
            *(output_seq + rev) = (info->tmp_sequence + i)->real / info->sampling_point_count;
        }
        else
        {
            *(output_seq + rev) = _fft_complex_module(info->tmp_sequence + i);
        }
    }
}

void fft_d2fft_complex(const void *fft_id,
                       const fft_complex *input_seq,
                       fft_complex *output_seq,
                       int mode)
{
    fft_info *info = (fft_info*)fft_id;
    int i, j, k, interval, half_n, gap;
    unsigned int rev, num;
    fft_complex tmp_complex;

    // Copy complex input array to temp complex array
    for (i = 0; i < info->sampling_point_count; ++i)
    {
        *(info->tmp_sequence + i) = *(input_seq + i);
    }

    // Butterfly Transform
    for (i = 0; i < info->sampling_point_level; i++)
    {
        interval = 1 << i;	// store count of groups
        half_n = 1 << (info->sampling_point_level - i);	// store length of every group
        for (j = 0; j < interval; j++)
        {
            // j means group j
            gap = j * half_n;
            for (k = 0; k < half_n / 2; k++)
            {
                _fft_complex_plus(info->tmp_sequence + k + gap,
                                  info->tmp_sequence + k + gap + half_n / 2,
                                  info->tmp_output_sequence + k + gap);
                _fft_complex_minus(info->tmp_sequence + k + gap,
                                   info->tmp_sequence + k + gap + half_n / 2,
                                   &tmp_complex);
                _fft_complex_multiple(&tmp_complex,
                                      info->polar_sequence + k * interval,
                                      info->tmp_output_sequence + k + gap + half_n / 2,
                                      mode);
            }
        }
        // copy result to input
        for (j = 0; j < info->sampling_point_count; j++)
        {
            *(info->tmp_sequence + j) = *(info->tmp_output_sequence + j);
        }
    }

    // invert the output code bit
    for (i = 0; i < info->sampling_point_count; i++)
    {
        rev = 0;
        num = (unsigned int)i;
        for (j = 0; j < info->sampling_point_level; j++)
        {
            rev <<= 1;
            rev |= num & 1;
            num >>= 1;
        }
        if(mode == IFFT)
        {
            (output_seq + rev)->real = (info->tmp_sequence + i)->real / info->sampling_point_count;
            (output_seq + rev)->imag = (info->tmp_sequence + i)->imag / info->sampling_point_count;
        }
        else
        {
            *(output_seq + rev) = *(info->tmp_sequence + i);
        }
    }
}
/*************************************************************************************************/
