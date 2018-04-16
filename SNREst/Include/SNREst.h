#ifndef SNREst_H_
#define SNREst_H_

//#define INFO_OUT

//#define HIGH_PASS

#define FILECLOSE(a) {if (a) {fclose(a); a = NULL;}}

#define MAX_STR_SIZE  4096

#define FRAME_MOVE  256

#define FRAME_LEN 512

#define FRAME_HALF (FRAME_LEN / 2 + 1)

#define LOG_FFT 9

#define FL 6 

#define FH FRAME_HALF

#define MIN_DIV0 1e-8f

#ifndef SHORT_MAX

#define SHORT_MAX ((1 << (sizeof(short) * 8 - 1)) - 1) 

#endif

#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define MAX(a,b) ((a) > (b) ? (a) : (b))

#endif


