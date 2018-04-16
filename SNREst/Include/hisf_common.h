/**************************************************************************//*
* @file  hisf_common.h
*
* @date  2017-12-29
*
* @brief  Copyright (c) 2017 - 2022. Horizon Robotics. All rights reserved.	
*	
*****************************************************************************/
#pragma once
#ifndef __HISF_COMMON_H_
#define __HISF_COMMON_H_
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "./hisf_common_types.h"
/* macro define*/
#ifndef NULL
#define NULL 0
#endif

/* 宏函数 */
#define  MULT_C(ar,ai, br, bi, cr, ci) {cr = (ar)*(br) - (ai)*(bi); ci = (ai)*(br) + (ar)*(bi); }
#define  ADD_C(ar,ai, br, bi, cr, ci)  {cr = ar + br; ci = ai + bi;}
 
/***********************************************************************************//**
 * @def TRIMFILEPATH(x) strrchr((x), ' ') ? strrchr((x), ' ') + 1 : (x)
 *
 * @brief A macro that defines trimfilepath.
 *
 * @author JiaJie.
 * @date 2018/1/17.
 *
 * @param x A void to process.
 *
 * @seealso .
 **************************************************************************************/
#ifdef linux
#define TRIMFILEPATH(x) strrchr((x), '/') ? strrchr((x), '/') + 1 : (x) 
#else
#define TRIMFILEPATH(x) strrchr((x), '\\') ? strrchr((x), '\\') + 1 : (x) 
#endif

/***********************************************************************************//**
 * @def CHECK_POINTER(pointer,err_type)
 *
 * @brief A macro that defines check pointer.
 *
 * @author JiaJie.
 * @date 2018/1/17.
 *
 * @param pointer  The pointer.
 * @param err_type Type of the error.
 *
 * @seealso .
 **************************************************************************************/
#define  CHECK_POINTER(pointer,err_type) {if((pointer) == NULL){HisfLogFatalLine(); return (err_type);}}
/***********************************************************************************//**
 * @def PRINT_ERROR_LOG(format,...)
 *
 * @brief A macro that defines print error log.
 *
 * @author JiaJie.
 * @date 2018/1/17.
 *
 * @param format Describes the format to use.
 * @param ...    Variable arguments providing additional information.
 *
 * @seealso .
 **************************************************************************************/
#define FPRINTF_STDERR_LOG(format,...) { fprintf(stderr, format, ##__VA_ARGS__);HisfLogFatal(format,__VA_ARGS__);}
/***********************************************************************************//**
 * @def PRINT_stdout_LOG(format,...)
 *
 * @brief A macro that defines print stdout log.
 *
 * @author JiaJie.
 * @date 2018/1/17.
 *
 * @param format Describes the format to use.
 * @param ...    Variable arguments providing additional information.
 *
 * @seealso .
 **************************************************************************************/
#define FPRINTF_STDOUT_LOG(format,...) { fprintf(stdout, format, ##__VA_ARGS__);HisfLogInfo(format,__VA_ARGS__);}
/* a normalization scale of audio data */
/**********************************************************************************************//**
 * @def AUDIO_NORM_SCALE.
 *
 * @brief   A macro that defines audio normalise scale.
 *
 * @author  zhuchangbao.
 * @date    2017/12/13.
 *
 * @seealso .
 **************************************************************************************************/                     
#define  MAX_WORD16 32767
#define  MIN_WORD16 -32768
#define  PI 3.1415926f
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define HISF_ABS(a) ((a) > 0 ? (a) : (-(a)))
#define MIN_DIV0       1e-8f 
#define MIN_DB_DELTA   1e-5f
/**********************************************************************************************//**
 * @typedef struct _HisfMem
 *
 * @brief   Defines an alias representing the hisf memory state.
 *
 * @seealso .
 **************************************************************************************************/
typedef struct _HisfMem
{
	/** @brief   Size of the memory total. */
	int32 mem_total_size;                                 
	/** @brief   Size of the available memory. */
	int32 available_mem_size;                                   
	/** @brief   The start point of memory */
	char *p_init;
	/** @brief   The current point of memory. */
	char *p_curr;
}THisfMem;


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
 * @param x  i/o: input/output vector 
 * @param n  vector length.
 * @param m  log2 of vector length.
 *
 * @seealso .
 **************************************************************************************/
void HisfCommonFFTRealValue(float32 x[],  const int16 n, const int16 m);
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
void HisfCommonIFFTRealValue( 	float32 io[],	const int16 n,	const int16 m );
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
void ComIfft(float32 *x, int32 N, int32 log2n, int32 data_type, float32 *out);
void ComFft(float32 *x, int32 N, int32 log2n, int32 data_type, float32 *out);
void HisfCommonAddWin(float32 *au, const float32 *win, int32 len, float32 *win_out);

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
void   HisfCommonGetPsd(float32 *f, int32 fft_len, float32 *psd);

/***************************************************************************************************//**
 table data
***************************************************************************************************/
extern const float32 kHammingWin512[512];
extern const float32 kHammingWin1024[1024];
extern const float32 kHammingWin2048[2048];
extern const int16 kFftReadIndexes[];               /* FFT */
extern const float32 kSinCosExt[];

#endif
