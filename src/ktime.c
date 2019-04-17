#if defined(WIN32) || defined(_WIN32)
#include <windows.h>

/*
 * timezone information is stored outside the kernel so tzp isn't used anymore.
 *
 * Note: this function is not for Win32 high precision timing purpose. See
 * elapsed_time().
 */
int gettimeofday(struct timeval * tp, struct timezone *tzp)
{
    FILETIME    file_time;
    SYSTEMTIME  system_time;
    ULARGE_INTEGER ularge;

    GetSystemTime(&system_time);
    SystemTimeToFileTime(&system_time, &file_time);
    ularge.LowPart = file_time.dwLowDateTime;
    ularge.HighPart = file_time.dwHighDateTime;

    tp->tv_sec = (long) ((ularge.QuadPart - epoch) / 10000000L);
    tp->tv_usec = (long) (system_time.wMilliseconds * 1000);

    return 0;
}

// taken from https://stackoverflow.com/questions/5272470/c-get-cpu-usage-on-linux-and-windows
double cputime()
{
	HANDLE hProcess = GetCurrentProcess();
	FILETIME ftCreation, ftExit, ftKernel, ftUser;
	SYSTEMTIME stKernel;
	SYSTEMTIME stUser;

	GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser);
	FileTimeToSystemTime(&ftKernel, &stKernel);
	FileTimeToSystemTime(&ftUser, &stUser);

	double kernelModeTime = ((stKernel.wHour * 60.) + stKernel.wMinute * 60.) + stKernel.wSecond * 1. + stKernel.wMilliseconds / 1000.;
	double userModeTime = ((stUser.wHour * 60.) + stUser.wMinute * 60.) + stUser.wSecond * 1. + stUser.wMilliseconds / 1000.;

	return kernelModeTime + userModeTime;
}
#else
#include <sys/resource.h>
#include <sys/time.h>

double cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

long peak_memory(void)
{
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);

    return r.ru_maxrss * 1024;
}

#endif /* WIN32 || _WIN32 */

double realtime(void)
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}


