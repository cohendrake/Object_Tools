/**************************************************************************//*
* @file  hisf_common_fft_rel.c
*
* @date  2017-12-29
*
* @brief  Copyright (c) 2017 - 2022. Horizon Robotics. All rights reserved.	
*	
*****************************************************************************/

#include <string.h>
#include <assert.h>
#include <memory.h> 
#include "../Include/hisf_common.h"
#include "../Include/hisf_common_types.h"

/*---------------------------------------------------------------------*
 * Local constants
 *---------------------------------------------------------------------*/

 /*---------------------------------------------------------------------*
 * Local constants
 *---------------------------------------------------------------------*/
#define N_MAX_FFT     4096
#define N_MAX_DIV2  (N_MAX_FFT>>1)
#define N_MAX_DIV4  (N_MAX_DIV2>>1)
#define INV_SQR2  0.70710676908493f

int16 Base2N(const uint16 a)
{
	int16 zeros = 0;

	if (a == 0) return 0;

	if (!(0xFF00 & (a << zeros))) zeros += 8;
	if (!(0xF000 & (a << zeros))) zeros += 4;
	if (!(0xC000 & (a << zeros))) zeros += 2;
	if (!(0x8000 & (a << zeros))) zeros += 1;

	zeros = 15 - zeros;
	return zeros;
	//return (a - ((unsigned int16)1 << zeros));
}

/***********************************************************************************//**
 * @fn void HisfCommonFFTRealValue( float32 x[], const int16 n, const int16 m )
 *
 * @brief *  Computes the split-radix FFT in place for the real-valued
 *  signal x of length n.  The algorithm has been ported from
 *  the Fortran code of [1].
 *
 *  The function  needs sine and cosine tables t_sin and t_cos,
 *  and the constant N_MAX_FFT.  The table  entries  are defined as
 *  sin(2*pi*i) and cos(2*pi*i) for i = 0, 1, ..., N_MAX_FFT-1. The
 *  implementation  assumes  that any entry  will not be needed
 *  outside the tables. Therefore, N_MAX_FFT and n must be properly
 *  set.  The function has been tested  with the values n = 16,
 *  32, 64, 128, 256, and N_MAX_FFT = 1280.
 *
 *  References
 *  [1] H.V. Sorensen,  D.L. Jones, M.T. Heideman, C.S. Burrus,
 *      "Real-valued fast  Fourier transform  algorithm,"  IEEE
 *      Trans. on Signal Processing,  Vol.35, No.6, pp 849-863,
 *      1987.
 *
 *  OUTPUT
 *      x[0:n-1]  Transform coeffients in the order re[0], re[1],
 *                ..., re[n/2], im[n/2-1], ..., im[1].
 *
 * @author zhuchangbao.
 * @date 2016/12/20.
 *
 * @param x i/o: input/output vector.
 * @param n vector length .
 * @param m log2 of vector length.
 *
 * @seealso .
 **************************************************************************************/
void HisfCommonFFTRealValue(float32 x[], 	const int16 n,	const int16 m )
{
	int16 i, j, k, n1, n2, n4;
	int16 step;
	float32 xt, t1, t2;
	float32 *x0, *x1, *x2;
	float32 *xi2, *xi3, *xi4, *xi1;
	const float32 *s, *c;
	const int16 *idx;
	uint16 base = (uint16)(N_MAX_FFT / n);
	int16 right_shift = Base2N(base);

	/* !!!! NMAX = 4096 is hard coded here  !!! */

	float32 temp[2 * N_MAX_FFT];

	if (n > 7 && 0 == N_MAX_FFT % n && (1 == base || 0 == base % 2))
	{
		idx = kFftReadIndexes;

		/* Combined Digit reverse counter & Length two butterflies */

		if (n == N_MAX_FFT)
		{
			x2 = temp;
			for (i = 0; i < N_MAX_FFT / 2; i++)
			{
				j = *idx++;
				k = *idx++;

				*x2++ = x[j] + x[k];
				*x2++ = x[j] - x[k];
			}
		}
		else
		{
			x2 = temp;
			for (i = 0; i < n / 2; i++)
			{
				j = *idx++;
				k = *idx++;

				*x2++ = x[j >> right_shift] + x[k >> right_shift];
				*x2++ = x[j >> right_shift] - x[k >> right_shift];
			}
		}

		/*-----------------------------------------------------------------*
		* 1st Stage Loop has been Unrolled because n4 is '1' and that
		* allows the elimination of the 'for_ (j = 1; j < n4; j++)' loop
		* and the associated pointers initialization.
		* Also, it allows to Put the Data from 'temp' back into 'x' due
		* to the previous Combined Digit Reverse and Length two butterflies
		*-----------------------------------------------------------------*/

		/*for_ (k = 2; k < 3; k++)*/
		{
			x0 = temp;
			x1 = x0 + 2;
			x2 = x;

			for (i = 0; i < n; i += 4)
			{
				*x2++ = *x0++ + *x1;    /* x[i] = xt + x[i+n2];    */
				*x2++ = *x0;
				*x2++ = *--x0 - *x1++;  /* x[i+n2] = xt - x[i+n2];      */
				*x2++ = -*x1;         /* x[i+n2+n4] = -x[i+n2+n4];     */

				x0 += 4;
				x1 += 3; /* x1 has already advanced */
			}
		}
	}
	else
	{
		/*-----------------------------------------------------------------*
		* Digit reverse counter
		*-----------------------------------------------------------------*/
		return;
		j = 0;
		x0 = &x[0];
		for (i = 0; i < n - 1; i++)
		{
			if (i < j)
			{
				xt = x[j];
				x[j] = *x0;
				*x0 = xt;
			}
			x0++;
			k = n / 2;
			while (k <= j)
			{
				j -= k;
				k = k >> 1;
			}
			j += k;
		}

		/*-----------------------------------------------------------------*
		* Length two butterflies
		*-----------------------------------------------------------------*/

		x0 = &x[0];
		x1 = &x[1];
		for (i = 0; i < n / 2; i++)
		{
			*x1 = *x0 - *x1;
			*x0 = *x0 * 2 - *x1;

			x0++;
			x0++;
			x1++;
			x1++;
		}

		/*-----------------------------------------------------------------*
		* 1st Stage Loop has been Unrolled because n4 is '1' and that
		* allows the elimination of the 'for_ (j = 1; j < n4; j++)' loop
		* and the associated pointers initialization.
		*-----------------------------------------------------------------*/

		/* for_ (k = 2; k < 3; k++) */
		{
			x0 = x;
			x1 = x0 + 2;

			for (i = 0; i < n; i += 4)
			{
				*x1 = *x0 - *x1;      /* x[i+n2] = xt - x[i+n2];      */
				*x0 = *x0 * 2 - *x1++;    /* x[i] = xt + x[i+n2];    */
				*x1 = -*x1;             /* x[i+n2+n4] = -x[i+n2+n4];     */

				x0 += 4;
				x1 += 3; /* x1 has already advanced */
			}
		}
	}

	/*-----------------------------------------------------------------*
	* Other butterflies
	*
	* The implementation described in [1] has been changed by using
	* table lookup for evaluating sine and cosine functions.  The
	* variable ind and its increment step are needed to access table
	* entries.  Note that this implementation assumes n4 to be so
	* small that ind will never exceed the table.  Thus the input
	* argument n and the constant N_MAX_FFT must be set properly.
	*-----------------------------------------------------------------*/

	n4 = 1;
	n2 = 2;
	n1 = 4;

	step = N_MAX_DIV4;

	for (k = 3; k <= m; k++)
	{
		step >>= (1);
		n4 <<= 1;
		n2 <<= 1;
		n1 <<= 1;

		x0 = x;
		x1 = x0 + n2;
		x2 = x1 + n4;

		for (i = 0; i < n; i += n1)
		{
			*x1 = *x0 - *x1;      /* x[i+n2] = xt - x[i+n2];      */
			*x0 = *x0 * 2 - *x1;      /* x[i] = xt + x[i+n2];    */
			*x2 = -*x2;             /* x[i+n2+n4] = -x[i+n2+n4];     */

			s = kSinCosExt;
			c = s + N_MAX_FFT / 4;          /* 1024/4 = 256, 256/4=64 */
			xi1 = x0;
			xi3 = xi1 + n2;
			xi2 = xi3;
			x0 += n1;
			xi4 = x0;

			for (j = 1; j < n4; j++)
			{
				xi3++;
				xi1++;
				xi4--;
				xi2--;
				c += step;
				s += step;          /* auto increment by ar0 */

				t1 = *xi3* *c + *xi4 * *s;    /* t1 = *xi3**(pt_c+ind) + *xi4**(pt_s+ind);   */
				t2 = *xi3 * *s - *xi4 * *c;    /* t2 = *xi3**(pt_s+ind) - *xi4**(pt_c+ind);     */

				*xi4 = *xi2 - t2;
				*xi2 = *xi1 - t1;
				*xi1 = *xi1 * 2 - *xi2;
				*xi3 = -2 * t2 - *xi4;
			}

			x1 += n1;
			x2 += n1;
		}
	}

	return;
}

/***********************************************************************************//**
 * @fn void HisfCommonIFFTRealValue( float32 io[], const int16 n, const int16 m )
 *
 * @brief  Calculate the inverse FFT of a real signal
 *
 * Based on the FORTRAN code from the article "Real-valued Fast Fourier Transform Algorithms"
 * by Sorensen, ... in IEEE Trans. on ASSP, Vol. ASSP-35, No. June 6th 1987.
 *
 * Input: the io[] signal containing the spectrum in the following order :
 *
 * Re[0], Re[1], ..  Re[n/2], Im[n/2-1], .. Im[1]
 *
 * @author zhuchangbao.
 * @date 2012/12/20.
 *
 * @param io input/output vector.
 * @param n  vector length.
 * @param m  log2 of vector length.
 *
 * @seealso .
 **************************************************************************************/
void HisfCommonIFFTRealValue(float32 io[], const int16 n, const int16 m)
{
	int16 i, j, k;
	int16 step;
	int16 n2, n4, n8, i0;
	int16 is, id;
	float32 *x, *xi0, *xi1, *xi2, *xi3, *xi4, *xup1, *xdn6, *xup3, *xdn8;
	float32 r1;
	float32 t1, t2, t3, t4, t5;
	float32 cc1, cc3, ss1, ss3;
	const float32 *s, *s3, *c, *c3;
	const int16 *idx;
	float32 temp[2 * N_MAX_FFT];
	float32 n_inv;
	uint16 base = (uint16)(N_MAX_FFT / n);
	int16 right_shift = Base2N(base);

	n_inv = 1.0f / n;

	/*-----------------------------------------------------------------*
	* IFFT
	*-----------------------------------------------------------------*/

	x = &io[-1];
	n2 = 2 * n;

	for (k = 1; k < m; k++)
	{
		is = 0;
		id = n2;
		n2 = n2 >> 1;
		n4 = n2 >> 2;
		n8 = n4 >> 1;

		while (is < n - 1)
		{
			xi1 = x + is + 1;
			xi2 = xi1 + n4;
			xi3 = xi2 + n4;
			xi4 = xi3 + n4;

			for (i = is; i < n; i += id)
			{
				t1 = *xi1 - *xi3;
				*xi1 += *xi3;
				*xi2 = 2.0f * *xi2;
				*xi3 = t1 - 2.0f * *xi4;
				*xi4 = t1 + 2.0f * *xi4;

				if (1 != n4)
				{
					t1 = (*(xi2 + n8) - *(xi1 + n8)) * INV_SQR2;
					t2 = (*(xi4 + n8) + *(xi3 + n8)) * INV_SQR2;

					*(xi1 + n8) += *(xi2 + n8);
					*(xi2 + n8) = *(xi4 + n8) - *(xi3 + n8);
					*(xi3 + n8) = (float32)(2.0f * (-t2 - t1));
					*(xi4 + n8) = (float32)(2.0f * (-t2 + t1));
				}
				xi1 += id;
				xi2 += id;
				xi3 += id;
				xi4 += id;
			}
			is = 2 * id - n2;
			id = 4 * id;
		}
		step = N_MAX_FFT / n2;

		s = kSinCosExt + step;
		c = s + N_MAX_FFT / 4;
		s3 = kSinCosExt + 3 * step;
		c3 = s3 + N_MAX_FFT / 4;

		for (j = 2; j <= n8; j++)
		{
			cc1 = *c;
			ss1 = *s;
			cc3 = *c3;
			ss3 = *s3;

			is = 0;
			id = 2 * n2;

			c += step;
			s += step;

			c3 += 3 * step;
			s3 += 3 * step;
			while (is < n - 1)
			{
				xup1 = x + j + is;
				xup3 = xup1 + 2 * n4;
				xdn6 = xup3 - 2 * j + 2;
				xdn8 = xdn6 + 2 * n4;

				for (i = is; i < n; i += id)
				{
					t1 = *xup1 - *xdn6;
					*xup1 = *xup1 + *xdn6;
					xup1 += n4;
					xdn6 -= n4;

					t2 = *xdn6 - *xup1;
					*xdn6 = *xup1 + *xdn6;

					xdn6 += n4;
					t3 = *xdn8 + *xup3;
					*xdn6 = *xdn8 - *xup3;

					xup3 += n4;
					xdn8 -= n4;

					t4 = *xup3 + *xdn8;
					*xup1 = *xup3 - *xdn8;

					t5 = t1 - t4;
					t1 = t1 + t4;
					t4 = t2 - t3;
					t2 = t2 + t3;
					*xup3 = t1 * cc3 - t2 * ss3;
					xup3 -= n4;
					*xup3 = t5 * cc1 + t4 * ss1;
					*xdn8 = -t4 * cc1 + t5 * ss1;

					xdn8 += n4;
					*xdn8 = t2*cc3 + t1*ss3;

					xup1 -= n4;
					xup1 += id;
					xup3 += id;
					xdn6 += id;
					xdn8 += id;
				}
				is = 2 * id - n2;
				id = 4 * id;
			}
		}
	}

	/*-----------------------------------------------------------------*
	* Length two butterflies
	*-----------------------------------------------------------------*/

	is = 1;
	id = 4;
	while (is < n)
	{
		xi0 = x + is;
		xi1 = xi0 + 1;

		for (i0 = is; i0 <= n; i0 += id)
		{
			r1 = *xi0;
			*xi0 = r1 + *xi1;
			*xi1 = r1 - *xi1;
			xi0 += id;
			xi1 += id;
		}
		is = 2 * id - 1;
		id = 4 * id;
	}

	/*-----------------------------------------------------------------*
	* Digit reverse counter
	*-----------------------------------------------------------------*/

	idx = kFftReadIndexes;
	xi0 = temp - 1;

	if (n > 7 && 0 == N_MAX_FFT % n && (1 == base || 0 == base % 2))
	{
		if (n == N_MAX_FFT)
		{
			for (i = 0; i < n; i++)
			{
				j = *idx++;
				temp[i] = x[1 + j];
			}
		}
		else
		{
			for (i = 0; i < n; i++)
			{
				j = *idx++;
				temp[i] = x[1 + (j >> right_shift)];
			}
		}
	}

	/*-----------------------------------------------------------------*
	* Normalization
	*-----------------------------------------------------------------*/

	for (i = 1; i <= n; i++)
	{
		x[i] = xi0[i] * n_inv;
	}

	return;
}

/**********************************************************************************************//**
 * @fn	void ComIfft(float32 *x, int32 N, int32 log2n, int32 data_type, float32 *out)
 *
 * @brief	公共的IFFT操作
 *
 * @author	zhuchangbao.
 * @date	2017/7/21.
 *
 * @param [int32]	x  	输入信号
 * @param	N		   	N,
 * @param	log2n	   	log2 N
 * @param	data_type  	1为紧密排列，即实部，虚部紧密排列
 * @param [out]	out	输出信号
 *
 * @seealso	.
 **************************************************************************************************/
void ComIfft(float32 *x, int32 N, int32 log2n, int32 data_type, float32 *out)
{
	int32 i;
	assert(x != out);
	memcpy(out, x, N*sizeof(out[0]));

	if (1 == data_type)
	{
		out[0] = 0;
		out[N / 2] = x[N];
		for (i = 1; i < N / 2; i++)
		{
			out[i] = x[2 * i];
			out[N - i] = x[2 * i + 1];
		}
	}

	HisfCommonIFFTRealValue(out, N, log2n);		
}

void ComFft(float32 *x, int32 N, int32 log2n, int32 data_type, float32 *out)
{
	int32 i;
	assert(x != out);


	HisfCommonFFTRealValue(x, N, log2n);

	

	memcpy(out, x, N*sizeof(out[0]));

	if (1 == data_type)
	{
		out[0] = x[0];
		out[1] = 0;
		out[N - 2] = x[N/2];
		out[N - 1] = 0;
		for (i = 1; i < N / 2; i++)
		{
			out[2*i] = x[i];
			out[2 * i + 1] = x[N - i];
		}
	}
}


/***********************************************************************************//**
 * @fn void HisfCommonAddWin(float32 *au, const float32 *win, int32 len, float32 *win_out)
 *
 * @brief Hisf common add window.
 *
 * @author zhuchangbao.
 * @date 2016/12/20.
 *
 * @param [in,out] au	   If non-null, the au.
 * @param win			   The window.
 * @param len			   The length.
 * @param [in,out] win_out If non-null, the window out.
 *
 * @seealso .
 **************************************************************************************/
void HisfCommonAddWin(float32 *au, const float32 *win, int32 len, float32 *win_out)
{
	int32 i;
	for (i = 0; i < len; i++)
	{
		win_out[i] = win[i] * au[i];
	}
}


/***********************************************************************************//**
 * @fn void HisfCommonGetPsd(float32 *f, int32 fft_len, float32 *psd)
 *
 * @brief Hisf common get psd.
 *
 * @author zhuchangbao.
 * @date 2016/12/20.
 *
 * @param [in,out] f   If non-null, a float32 to process.
 * @param fft_len	   Length of the FFT.
 * @param [in,out] psd If non-null, the psd.
 *
 * @seealso .
 **************************************************************************************/
void HisfCommonGetPsd(float32 *f, int32 fft_len, float32 *psd)
{
	int32 i;

	psd[0] = f[0] * f[0];
	psd[fft_len/2] = f[fft_len/2] * f[fft_len/2];

	for ( i = 1; i < fft_len/2; i++)
	{
		psd[i] = f[i] * f[i] + f[fft_len - i] * f[fft_len - i];
	}
}
