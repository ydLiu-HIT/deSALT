#if defined(WIN32) || defined(_WIN32)
#include <windows.h>

struct timezone
{
  __int32  tz_minuteswest; /* minutes W of Greenwich */
  int  tz_dsttime;     /* type of dst correction */
};

/*
 * gettimeofday.c
 *    Win32 gettimeofday() replacement
 *    taken from PostgreSQL, according to
 *    https://stackoverflow.com/questions/1676036/what-should-i-use-to-replace-gettimeofday-on-windows
 *
 * src/port/gettimeofday.c
 *
 * Copyright (c) 2003 SRA, Inc.
 * Copyright (c) 2003 SKC, Inc.
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation for any purpose, without fee, and without a
 * written agreement is hereby granted, provided that the above
 * copyright notice and this paragraph and the following two
 * paragraphs appear in all copies.
 *
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE TO ANY PARTY FOR DIRECT,
 * INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
 * LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
 * DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE AUTHOR SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS
 * IS" BASIS, AND THE AUTHOR HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE,
 * SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
 */

/* FILETIME of Jan 1 1970 00:00:00. */
static const unsigned __int64 epoch = ((unsigned __int64) 116444736000000000ULL);

/*
 * timezone information is stored outside the kernel so tzp isn't used anymore.
 *
 * Note: this function is not for Win32 high precision timing purpose. See
 * elapsed_time().
 */
int gettimeofday(struct timeval * tp, struct timezone *tzp);

// taken from https://stackoverflow.com/questions/5272470/c-get-cpu-usage-on-linux-and-windows
double cputime();
#else
#include <sys/resource.h>
#include <sys/time.h>

double cputime(void);

long peak_memory(void);

#endif /* WIN32 || _WIN32 */

double realtime(void);
