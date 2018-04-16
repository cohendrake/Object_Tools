#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "SNREst.h"
#include "wav_head.h"

/*--------------------------------------------------------------------------*
* Function  sig_wave_stream_init                                            *
*                                                                           *                                                                           *
* -> 获取原始语音数据头                                                     *
*                                                                           *
*---------------------------------------------------------------------------*/
int sig_wave_stream_init(hdr_t *wave_head, char *filename)
{
	FILE *fd = NULL;
	int head_len = 0;

	// open input file 
	if (NULL == (fd = fopen(filename, "rb")))
	{
		fprintf(stderr, "sig_wave_stream_init(): cannot open file %s\n", filename);
		return(-1);
	}

	{
		if (1 != fread(wave_head->riff, 4, 1, fd)) return (-1); 
		if (1 != fread(&wave_head->totsize, 4, 1, fd)) return (-1);
		if (1 != fread(wave_head->wave, 4, 1, fd)) return (-1);
		if (1 != fread(wave_head->fmtstr, 4, 1, fd)) return (-1);
		if (1 != fread(&wave_head->fmtsize, 4, 1, fd)) return (-1);
		if (1 != fread(&wave_head->compressionCode, 2, 1, fd)) return (-1);
		if (1 != fread(&wave_head->numchans, 2, 1, fd)) return (-1);
		if (1 != fread(&wave_head->Fs, 4, 1, fd)) return (-1);
		if (1 != fread(&wave_head->nbytespersec, 4, 1, fd)) return (-1);
		if (1 != fread(&wave_head->nbytespersample, 2, 1, fd)) return (-1);
		if (1 != fread(&wave_head->nbitspersample, 2, 1, fd)) return (-1);
		fseek(fd, wave_head->fmtsize - 16, SEEK_CUR);
		if (1 != fread(wave_head->data, 4, 1, fd)) return (-1);

		if (0 == strncmp(wave_head->data, "fact", 4))
		{
			fseek(fd, 8, SEEK_CUR);
			if (1 != fread(wave_head->data, 4, 1, fd)) return (-1);
		}
		
		if (1 != fread(&wave_head->datsize, 4, 1, fd)) return (-1);
		head_len = ftell(fd);
		FILECLOSE(fd);
	}

	if (strncmp(wave_head->riff, "RIFF", 4) != 0 || strncmp(wave_head->wave, "WAVE", 4) || strncmp(wave_head->data, "data", 4))
	{
		fprintf(stderr, "sig_wave_stream_init(): stream %s not in WAVE format\n", filename);
		FILECLOSE(fd);
		return(-1);
	}

	return(head_len);
}

/*--------------------------------------------------------------------------*
* Function  sig_wave_head_disp                                              *
*                                                                           *                                                                           *
* -> 显示wav头信息                                                          *
*                                                                           *
*---------------------------------------------------------------------------*/
void sig_wave_head_disp(hdr_t *wave_head)
{
	if (wave_head)
	{
		fprintf(stdout,"***********inwav info**************\n");
		fprintf(stdout,"[Channels]: %u\n",wave_head->numchans);
		fprintf(stdout,"[SampleRate]: %u Hz\n",wave_head->Fs);
		fprintf(stdout,"[Datasize]: %u Bytes\n",wave_head->datsize);
		fprintf(stdout,"[Duration]: %.1f s\n",(float)wave_head->datsize / (float)wave_head->nbytespersec);
		fprintf(stdout,"***********************************\n");
	}
	else
	{
		fprintf(stderr, "Invalid wav head\n");
	}
}
