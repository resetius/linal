#include <float.h>

#ifdef WIN32
#include <windows.h>
extern "C" void set_fpe_except()
{
	int cw = _controlfp (0, 0);
	cw &= ~ (EM_OVERFLOW | EM_UNDERFLOW | EM_ZERODIVIDE | EM_DENORMAL | EM_INVALID);
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

extern "C" void set_fpe_except()
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
