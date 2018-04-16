#ifndef WAV_HEAD_H_
#define WAV_HEAD_H_


typedef struct hdr {
	char riff[4];                          // String 'RIFF' without \0                
	unsigned int totsize;                  // Total file size - 8                     
	char wave[4];                          // String 'WAVE' without \0                
	char fmtstr[4];                        // String 'fmt_' without \0                
	unsigned int fmtsize;                  // Total size of format chunk (Usually 0x10)  from  compressionCode to nbitspersample         
	unsigned short compressionCode;        // Always 0x01                             
	unsigned short numchans;               // Number of channels                      
	unsigned int Fs;                       // Sample rate (in Hz)                     
	unsigned int nbytespersec;             // number of bytes/seconde                 
	unsigned short nbytespersample;        // number of bytes/sample                  
								           // (1=8 bit mono, 2=8 bit stereo or 16 bit mono, 4= 16 bit stereo)   
	unsigned short nbitspersample;         // number of bits/samples                  
	char data[4];                          // string 'data' without \0                
	unsigned int datsize;                  // number of data bytes (not samples)      
} hdr_t;

/* sig_wave_stream_init declaration */
extern int sig_wave_stream_init(hdr_t *wave_head, char *filename);

/* sig_wave_head_disp declaration */
extern void sig_wave_head_disp(hdr_t *wave_head);

#endif
