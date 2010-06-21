/* -*- charset: utf-8 -*- */
/*$Id$*/

/* Copyright (c) 2009 Alexey Ozeritsky (Алексей Озерицкий)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "linal_util.h"

extern "C"
{

#ifdef WIN32
#include <windows.h>
	void set_fpe_except()
	{
		int cw = _controlfp (0, 0);
		cw &= ~ (EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL);
		_controlfp (cw, MCW_EM);
	}
#else


	/*
	For reference, the layout of the MXCSR register:
	FZ:RC:RC:PM:UM:OM:ZM:DM:IM:Rsvd:PE:UE:OE:ZE:DE:IE
	15 14 13 12 11 10  9  8  7   6   5  4  3  2  1  0

	And the layout of the 387 FPU control word register:
	Rsvd:Rsvd:Rsvd:X:RC:RC:PC:PC:Rsvd:Rsvd:PM:UM:OM:ZM:DM:IM
	15   14   13 12 11 10  9  8   7    6   5  4  3  2  1  0

	Where:
	Rsvd - Reserved
	FZ   - Flush to Zero
	RC   - Rounding Control
	PM   - Precision Mask
	UM   - Underflow Mask
	OM   - Overflow Mask
	ZM   - Zerodivide Mask
	DM   - Denormal Mask
	IM   - Invalid Mask
	PE   - Precision Exception
	UE   - Underflow Exception
	OE   - Overflow Exception
	ZE   - Zerodivide Exception
	DE   - Denormal Exception
	IE   - Invalid Exception
	X    - Infinity control (unused on 387 and higher)
	PC   - Precision Control

	Source: Intel Architecture Software Development Manual, Volume 1, Basic Architecture
	*/

	void set_fpe_except()
	{
		//http://www.website.masmforum.com/tutorials/fptute/fpuchap1.htm
		unsigned m = 0;
asm ("fstcw %0" : : "m" (*&m) );
		asm ("fwait");
		m ^= 0x1f; //turn on all exceptions!
		m |= (1 << 8);
		m |= (1 << 9);
		m ^= (1 << 8);
		m ^= (1 << 9);
		m |= (1 << 9); // double precision
asm ("fldcw %0" : : "m" (*&m) );
		asm ("fwait");

		m = 1 << 12;
asm ("ldmxcsr %0" : : "m" (*&m) );
	}
#endif

#include <time.h>

#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

#ifdef _MSC_VER
#include <windows.h>

	struct timezone
	{
		int  tz_minuteswest; /* minutes W of Greenwich */
		int  tz_dsttime;     /* type of dst correction */
	};

	int gettimeofday (struct timeval *tv, struct timezone *tz)
	{
		FILETIME ft;
		unsigned __int64 tmpres = 0;
		static int tzflag;

		if (NULL != tv)
		{
			GetSystemTimeAsFileTime (&ft);

			tmpres |= ft.dwHighDateTime;
			tmpres <<= 32;
			tmpres |= ft.dwLowDateTime;

			/*converting file time to unix epoch*/
			tmpres /= 10;  /*convert into microseconds*/
			tmpres -= DELTA_EPOCH_IN_MICROSECS;
			tv->tv_sec = (long) (tmpres / 1000000UL);
			tv->tv_usec = (long) (tmpres % 1000000UL);
		}

		if (NULL != tz)
		{
			if (!tzflag)
			{
				_tzset();
				tzflag++;
			}
			tz->tz_minuteswest = _timezone / 60;
			tz->tz_dsttime = _daylight;
		}

		return 0;
	}
#else
#include <sys/time.h>
#endif

	double get_full_time()
	{
		struct timeval tv;
		gettimeofday (&tv, 0);
		return (double) (tv.tv_sec * 100.0 + tv.tv_usec / 10000.0);
	}

} /* extern "C" */

