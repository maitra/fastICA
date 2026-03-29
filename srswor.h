#ifndef __H_SRSWOR__
#define __H_SRSWOR__
#include "constants.h"

int srswor(SIZE_T n, SIZE_T k, SIZE_T *y);
int srsworcond(SIZE_T n, SIZE_T k, SIZE_T *y, int (*call_back)(SIZE_T *,
	SIZE_T, void *data), void *data);
#endif
