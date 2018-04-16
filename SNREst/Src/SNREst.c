#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "wav_head.h"
#include "hisf_common_types.h"
#include "hisf_common.h"
#include "SNREst.h"


/*global virables*/
FILE *f_inwav = NULL;
FILE *f_vad = NULL;


int main(int argc, char *argv[])
{
	/*local virables*/
	hdr_t hd_in;
	int hd_len = 0;
	int flen = 0;
	int vad_len = 0;
	int n_samples, i;
	char inwav_name[MAX_STR_SIZE] = { 0 };
	char vad_name[MAX_STR_SIZE] = { 0 };
	char* ptr;
	short *data;
	float fft[FRAME_LEN] = { 0 };
	float fft_win[FRAME_LEN] = { 0 };
	float psd[FRAME_HALF] = { 0 };
	char str[10];
	int vadflag;
	float signal = 0.0f;
	float noise = 0.0f;
	float signal_f = 0.0f;
	float noise_f = 0.0f;
	float snr;
	float snr_f;
	unsigned int frame_s = 0;
	unsigned int frame_n = 0;
	
	// arguments check
	if (argc == 3)
	{
		strcpy(inwav_name, argv[1]);
		strcpy(vad_name, argv[2]);
	}
	else
	{
		fprintf(stderr, "Usage: SNREst inwav invad\n");
		FILECLOSE(f_inwav);
		FILECLOSE(f_vad);
		exit(1);
	}
	
	// wav suffix check
	if (NULL == (ptr = strrchr(inwav_name, '.')))
	{
		fprintf(stderr,"Can't get the suffix of input file\n");
		FILECLOSE(f_inwav);
		FILECLOSE(f_vad);
		exit(1);
	}
	else if(strcmp(ptr, ".wav"))
	{
		fprintf(stderr,"Please input a wav file\n");
		FILECLOSE(f_inwav);
		FILECLOSE(f_vad);
		exit(1);
	}
	else
	{
		if (-1 == (hd_len = sig_wave_stream_init(&hd_in, inwav_name)) )
		{
			fprintf(stderr,"Input wav head is unrecognized:%s\n", inwav_name);
			FILECLOSE(f_inwav);
			FILECLOSE(f_vad);
			exit(1);
		}
		else
		{
#ifdef INFO_OUT
			sig_wave_head_disp(&hd_in);
#endif
			
			if (NULL == (f_inwav = fopen(inwav_name, "rb")))
			{
				fprintf(stderr,"Open input wav file %s failed!\n", inwav_name);
				FILECLOSE(f_inwav);
				FILECLOSE(f_vad);
				exit(1);
			}
			else
			{
				fseek(f_inwav, 0L, SEEK_END);
				flen = ftell(f_inwav);
				flen -= hd_len;
				if (flen != hd_in.datsize) 
				{
					fprintf(stderr,"Data Size [%u] in wav head is inconsistent with real size [%d]\n", hd_in.datsize,flen);
				}
			}
		}
	}
	
	//vad check
	if (NULL == (f_vad = fopen(vad_name, "rb")))
	{
		fprintf(stderr,"Open input vad file %s failed!\n", vad_name);
		FILECLOSE(f_inwav);
		FILECLOSE(f_vad);
		exit(1);
	}
	else
	{
		fseek(f_vad, 0L, SEEK_END);
		vad_len = ftell(f_vad) / 2; /*vad file must be in Unix format*/
		if ((vad_len + 1) * FRAME_MOVE * hd_in.nbytespersample < flen)
		{
#ifdef INFO_OUT
			fprintf(stderr,"[Warning]: Vad file is shorter than wav\n");
#endif
		}
		else if ((vad_len - 1) * FRAME_MOVE * hd_in.nbytespersample > flen)
		{
#ifdef INFO_OUT
			fprintf(stderr,"[Warning]: Vad file is longer than wav\n");
#endif
		}
	}
	
	//SNR Est
	fseek(f_inwav, hd_len, SEEK_SET);
	fseek(f_vad, 0L, SEEK_SET);
	if (NULL == (data = (short *)malloc(FRAME_MOVE * hd_in.numchans * sizeof(short))))
	{
		fprintf(stderr,"[Error]: Allocate memory for audio buffer failed!\n");
		FILECLOSE(f_inwav);
		FILECLOSE(f_vad);
		exit(1);
	}

	//fprintf(stdout,"[Debug]: SHORT_MAX: %d\n",SHORT_MAX);
	while ((n_samples = fread(data, sizeof(short), FRAME_MOVE * hd_in.numchans, f_inwav)) > 0 && fgets(str, 10, f_vad) != NULL)
	{
		vadflag = atoi(str);
		memmove(fft, fft + FRAME_MOVE, (FRAME_LEN - FRAME_MOVE) * sizeof(fft[0]));
		if (vadflag)
		{
			frame_s++;
			for (i=0; i<FRAME_MOVE; i++)
			{
				signal += data[i * hd_in.numchans] * data[i * hd_in.numchans];
#ifdef HIGH_PASS
				fft[FRAME_LEN - FRAME_MOVE + i] = (float)data[i * hd_in.numchans] / SHORT_MAX;
#endif
			}
		}
		else
		{
			frame_n++;
			for (i=0; i<FRAME_MOVE; i++)
			{
				noise += data[i * hd_in.numchans] * data[i * hd_in.numchans];
#ifdef HIGH_PASS
				fft[FRAME_LEN - FRAME_MOVE + i] = (float)data[i * hd_in.numchans] / SHORT_MAX;
#endif
			}
		}
#ifdef HIGH_PASS
		HisfCommonAddWin(fft, kHammingWin512, FRAME_LEN, fft_win);
		HisfCommonFFTRealValue(fft_win, FRAME_LEN, LOG_FFT);
		HisfCommonGetPsd(fft_win, FRAME_LEN, psd);

		if (vadflag)
		{
			for (i=FL; i<FH; i++)	
			{
				signal_f += psd[i];
			}
		}
		else
		{
			for (i=FL; i<FH; i++)
			{
				noise_f += psd[i];
			}
		}
#endif
	}
	
	snr = 10.0f * log10f(signal * frame_n / (noise * frame_s + MIN_DIV0));
#ifdef HIGH_PASS
	snr_f = 10.0f * log10f(signal_f * frame_n / (noise_f * frame_s + MIN_DIV0));
#endif
	
	snr = MAX(MIN(snr,40.0f), -40.0f);
#ifdef HIGH_PASS
	snr_f = MAX(MIN(snr_f,40.0f), -40.0f);
#endif
#ifdef INFO_OUT
	fprintf(stdout,"[S]: %.1f\n", signal);
	fprintf(stdout,"[N]: %.1f\n", noise);
	fprintf(stdout,"[SNR]: %.1f dB\n", snr);
	fprintf(stdout,"[S_f]: %.1f\n", signal_f);
	fprintf(stdout,"[N_f]: %.1f\n", noise_f);
	fprintf(stdout,"[SNR_f]: %.1f dB\n", snr_f);	
#else
	fprintf(stdout,"%.1f\n", snr);
	#ifdef HIGH_PASS
	fprintf(stdout,"%.1f\n", snr_f);	
	#endif
#endif
	FILECLOSE(f_inwav);
	FILECLOSE(f_vad);
	if (data) free(data);
	return 0;
	
}
