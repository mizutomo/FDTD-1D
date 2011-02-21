#include <sys/time.h>
#include <sys/resource.h>
#include "utility.h"

// 現在時刻を返す
double get_current_time_by_sec()
{
	struct rusage t;
	struct timeval tv;
	getrusage(RUSAGE_SELF, &t);
	tv = t.ru_utime;
	return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

// メモリ消費量の取得
long int get_use_memory_size_from_mac()
{
	struct rusage t;
	getrusage(RUSAGE_SELF, &t);
	return t.ru_maxrss;
}
