#include <sys/time.h>
#include <sys/resource.h>
#include "utility.h"

// 現在時刻を返す
void get_current_time_by_sec(double* time)
{
	struct rusage t;
	struct timeval tv;
	getrusage(RUSAGE_SELF, &t);
	tv = t.ru_utime;
	*time = tv.tv_sec + (double)tv.tv_usec*1e-6;
}

// メモリ消費量の取得
void get_use_memory_size_from_mac(long int* mem)
{
	struct rusage t;
	getrusage(RUSAGE_SELF, &t);
	*mem = t.ru_maxrss;
}
