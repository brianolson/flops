package main

import "fmt"

//import "os"
import "syscall"

/* ported to Go, Brian Olson, http://bolson.org/, 2014-10 */
/*--------------------- Start flops.c source code ----------------------*/

/*****************************/
/*          flops.c          */
/* Version 2.0,  18 Dec 1992 */
/*         Al Aburto         */
/*      aburto@nosc.mil      */
/*****************************/

/*
   Flops.c is a 'c' program which attempts to estimate your systems
   floating-point 'MFLOPS' rating for the FADD, FSUB, FMUL, and FDIV
   operations based on specific 'instruction mixes' (discussed below).
   The program provides an estimate of PEAK MFLOPS performance by making
   maximal use of register variables with minimal interaction with main
   memory. The execution loops are all small so that they will fit in
   any cache. Flops.c can be used along with Linpack and the Livermore
   kernels (which exersize memory much more extensively) to gain further
   insight into the limits of system performance. The flops.c execution
   modules also include various percent weightings of FDIV's (from 0% to
   25% FDIV's) so that the range of performance can be obtained when
   using FDIV's. FDIV's, being computationally more intensive than
   FADD's or FMUL's, can impact performance considerably on some systems.

   Flops.c consists of 8 independent modules (routines) which, except for
   module 2, conduct numerical integration of various functions. Module
   2, estimates the value of pi based upon the Maclaurin series expansion
   of atan(1). MFLOPS ratings are provided for each module, but the
   programs overall results are summerized by the MFLOPS(1), MFLOPS(2),
   MFLOPS(3), and MFLOPS(4) outputs.

   The MFLOPS(1) result is identical to the result provided by all
   previous versions of flops.c. It is based only upon the results from
   modules 2 and 3. Two problems surfaced in using MFLOPS(1). First, it
   was difficult to completely 'vectorize' the result due to the
   recurrence of the 's' variable in module 2. This problem is addressed
   in the MFLOPS(2) result which does not use module 2, but maintains
   nearly the same weighting of FDIV's (9.2%) as in MFLOPS(1) (9.6%).
   The second problem with MFLOPS(1) centers around the percentage of
   FDIV's (9.6%) which was viewed as too high for an important class of
   problems. This concern is addressed in the MFLOPS(3) result where NO
   FDIV's are conducted at all.

   The number of floating-point instructions per iteration (loop) is
   given below for each module executed:

   MODULE   FADD   FSUB   FMUL   FDIV   TOTAL  Comment
     1        7      0      6      1      14   7.1%  FDIV's
     2        3      2      1      1       7   difficult to vectorize.
     3        6      2      9      0      17   0.0%  FDIV's
     4        7      0      8      0      15   0.0%  FDIV's
     5       13      0     15      1      29   3.4%  FDIV's
     6       13      0     16      0      29   0.0%  FDIV's
     7        3      3      3      3      12   25.0% FDIV's
     8       13      0     17      0      30   0.0%  FDIV's

   A*2+3     21     12     14      5      52   A=5, MFLOPS(1), Same as
	   40.4%  23.1%  26.9%  9.6%          previous versions of the
						flops.c program. Includes
						only Modules 2 and 3, does
						9.6% FDIV's, and is not
						easily vectorizable.

   1+3+4     58     14     66     14     152   A=4, MFLOPS(2), New output
   +5+6+    38.2%  9.2%   43.4%  9.2%          does not include Module 2,
   A*7                                         but does 9.2% FDIV's.

   1+3+4     62      5     74      5     146   A=0, MFLOPS(3), New output
   +5+6+    42.9%  3.4%   50.7%  3.4%          does not include Module 2,
   7+8                                         but does 3.4% FDIV's.

   3+4+6     39      2     50      0      91   A=0, MFLOPS(4), New output
   +8       42.9%  2.2%   54.9%  0.0%          does not include Module 2,
						and does NO FDIV's.

   NOTE: Various timer routines are included as indicated below. The
	timer routines, with some comments, are attached at the end
	of the main program.

   NOTE: Please do not remove any of the printouts.

   EXAMPLE COMPILATION:
   UNIX based systems
       cc -DUNIX -O flops.c -o flops
       cc -DUNIX -DROPT flops.c -o flops
       cc -DUNIX -fast -O4 flops.c -o flops
       .
       .
       .
     etc.

   Al Aburto
   aburto@nosc.mil
*/

/***************************************************************/
/* Timer options. You MUST uncomment one of the options below  */
/* or compile, for example, with the '-DUNIX' option.          */
/***************************************************************/
/* #define Amiga       */
//#define UNIX
/* #define UNIX_Old    */
/* #define VMS         */
/* #define BORLAND_C   */
/* #define MSC         */
/* #define MAC         */
/* #define IPSC        */
/* #define FORTRAN_SEC */
/* #define GTODay      */
/* #define CTimer      */
/* #define UXPM        */
/* #define MAC_TMgr*/
/* #define PARIX       */
/* #define POSIX       */
/* #define WIN32       */
/* #define POSIX1      */
/***********************/

/*
#if (! defined(Amiga)) && (! defined(UNIX)) && (! defined(UNIX_Old)) && (! defined(VMS)) && (! defined(BORLAND_C)) && (! defined(MSC)) && (! defined(MAC)) && (! defined(IPSC)) && (! defined(FORTRAN_SEC)) && (! defined(GTODay)) && (! defined(CTimer)) && (! defined(UXPM)) && (! defined(MAC_TMgr)) && (! defined(PARIX)) && (! defined(POSIX)) && (! defined(WIN32)) && (! defined(POSIX1))
#warning using default UNIX timer
#define UNIX
#endif

int dtime( double* );

#include <stdio.h>
#include <math.h>
*/
/* 'Uncomment' the line below to run   */
/* with 'register double' variables    */
/* defined, or compile with the        */
/* '-DROPT' option. Don't need this if */
/* registers used automatically, but   */
/* you might want to try it anyway.    */
/* #define ROPT */

/* Variables needed for 'dtime()'.     */
var nulltime float64
var TimeArray []float64 = []float64{0, 0, 0}

var TLimit float64 /* Threshold to determine Number of    */
/* Loops to run. Fixed at 15.0 seconds.*/

var T [36]float64 /* Global Array used to hold timing    */
/* results and other information.      */

var sa, sb, sc, sd, one, two, three float64
var four, five, piref, piprg float64
var scale, pierr float64

const A0 = 1.0
const A1 = -0.1666666666671334
const A2 = 0.833333333809067E-2

var A3 float64 = 0.198412715551283E-3

const A4 = 0.27557589750762E-5

var A5 float64 = 0.2507059876207E-7

const A6 = 0.164105986683E-9

const B0 = 1.0
const B1 = -0.4999999999982
const B2 = 0.4166666664651E-1
const B3 = -0.1388888805755E-2
const B4 = 0.24801428034E-4
const B5 = -0.2754213324E-6
const B6 = 0.20189405E-8

const C0 = 1.0
const C1 = 0.99999999668
const C2 = 0.49999995173
const C3 = 0.16666704243
const C4 = 0.4166685027E-1
const C5 = 0.832672635E-2
const C6 = 0.140836136E-2
const C7 = 0.17358267E-3
const C8 = 0.3931683E-4

const D1 = 0.3999999946405E-1
const D2 = 0.96E-3
const D3 = 0.1233153E-5

const E2 = 0.48E-3
const E3 = 0.411051E-6

type long int64

func main() {
	// int argc, char** argv

	var s, u, v, w, x float64

	var loops, NLimit long
	var i, m, n long

	fmt.Printf("\n")
	fmt.Printf("   FLOPS Go Program (Float64 Precision), V2.0 18 Dec 1992\n\n")

	/****************************/
	loops = 15625 /* Initial number of loops. */
	/*     DO NOT CHANGE!       */
	/****************************/

	/****************************************************/
	/* Set Variable Values.                             */
	/* T[1] references all timing results relative to   */
	/* one million loops.                               */
	/*                                                  */
	/* The program will execute from 31250 to 512000000 */
	/* loops based on a runtime of Module 1 of at least */
	/* TLimit = 15.0 seconds. That is, a runtime of 15  */
	/* seconds for Module 1 is used to determine the    */
	/* number of loops to execute.                      */
	/*                                                  */
	/* No more than NLimit = 512000000 loops are allowed*/
	/****************************************************/

	T[1] = 1.0E+06 / float64(loops)

	TLimit = 15.0
	NLimit = 512000000

	piref = 3.14159265358979324
	one = 1.0
	two = 2.0
	three = 3.0
	four = 4.0
	five = 5.0
	scale = one

	fmt.Printf("   Module     Error        RunTime      MFLOPS\n")
	fmt.Printf("                            (usec)\n")
	/*************************/
	/* Initialize the timer. */
	/*************************/

	dtime(TimeArray)
	dtime(TimeArray)

	/*******************************************************/
	/* Module 1.  Calculate integral of df(x)/f(x) defined */
	/*            below.  Result is ln(f(1)). There are 14 */
	/*            float64 precision operations per loop     */
	/*            ( 7 +, 0 -, 6 *, 1 / ) that are included */
	/*            in the timing.                           */
	/*            50.0% +, 00.0% -, 42.9% *, and 07.1% /   */
	/*******************************************************/
	n = loops
	sa = 0.0
	x = one / float64(n)
	s = 0.0

	for sa < TLimit {
		n = 2 * n
		x = one / float64(n) /*********************/
		s = 0.0              /*  Loop 1.          */
		v = 0.0              /*********************/
		w = one

		dtime(TimeArray)
		for i = 1; i <= n-1; i++ {
			v = v + w
			u = v * x
			s = s + (D1+u*(D2+u*D3))/(w+u*(D1+u*(E2+u*E3)))
		}
		dtime(TimeArray)
		sa = TimeArray[1]

		if n == NLimit {
			break
		}
		/* fmt.Printf(" %10ld  %12.5lf\n",n,sa); */
	}

	scale = 1.0E+06 / float64(n)
	T[1] = scale
	fmt.Printf(" n=%d  sa=%f  scale=%f\n", n, sa, scale)

	/****************************************/
	/* Estimate nulltime ('for' loop time). */
	/****************************************/
	dtime(TimeArray)
	for i = 1; i <= n-1; i++ {
	}
	dtime(TimeArray)
	nulltime = T[1] * TimeArray[1]
	if nulltime < 0.0 {
		nulltime = 0.0
	}

	T[2] = T[1]*sa - nulltime

	sa = (D1 + D2 + D3) / (one + D1 + E2 + E3)
	sb = D1

	T[3] = T[2] / 14.0               /*********************/
	sa = x * (sa + sb + two*s) / two /* Module 1 Results. */
	sb = one / sa                    /*********************/
	n = long(float64(40000*long(sb)) / scale)
	//fmt.Printf( "x=%lf\ts=%lf\tsa=%lf\tsb=%lf\tscale=%lf\tn=%ld\n", x, s, sa, sb, scale, n );
	sc = sb - 25.2
	T[4] = one / T[3]
	/********************/
	/*  DO NOT REMOVE   */
	/*  THIS PRINTOUT!  */
	/********************/
	//fmt.Printf("     1   %13.4le  %10.4lf  %10.4lf\n",sc,T[2],T[4]);
	fmt.Printf("     1   %13.4e  %10.4f  %10.4f\n", sc, T[2], T[4])

	m = n

	/*******************************************************/
	/* Module 2.  Calculate value of PI from Taylor Series */
	/*            expansion of atan(1.0).  There are 7     */
	/*            float64 precision operations per loop     */
	/*            ( 3 +, 2 -, 1 *, 1 / ) that are included */
	/*            in the timing.                           */
	/*            42.9% +, 28.6% -, 14.3% *, and 14.3% /   */
	/*******************************************************/

	//fmt.Printf("starting module 2, m=%ld\n", m );
	s = -five /********************/
	sa = -one /* Loop 2.          */
	/********************/
	dtime(TimeArray)
	for i = 1; i <= m; i++ {
		s = -s
		sa = sa + s
	}
	dtime(TimeArray)
	T[5] = T[1] * TimeArray[1]
	if T[5] < 0.0 {
		T[5] = 0.0
	}

	sc = float64(m)

	u = sa  /*********************/
	v = 0.0 /* Loop 3.           */
	w = 0.0 /*********************/
	x = 0.0

	dtime(TimeArray)
	for i = 1; i <= m; i++ {
		s = -s
		sa = sa + s
		u = u + two
		x = x + (s - u)
		v = v - s*u
		w = w + s/u
	}
	dtime(TimeArray)
	T[6] = T[1] * TimeArray[1]

	T[7] = (T[6] - T[5]) / 7.0 /*********************/
	m = long(sa * x / sc)      /*  PI Results       */
	//fmt.Printf("sa=%lf\tx=%lf\tsc\%lf\tm=%ld\n", sa, x, sc, m );
	sa = four * w / five /*********************/
	sb = sa + five/v
	sc = 31.25
	piprg = sb - sc/(v*v*v)
	pierr = piprg - piref
	T[8] = one / T[7]
	/*********************/
	/*   DO NOT REMOVE   */
	/*   THIS PRINTOUT!  */
	/*********************/
	fmt.Printf("     2   %13.4e  %10.4f  %10.4f\n", pierr, T[6]-T[5], T[8])

	/*******************************************************/
	/* Module 3.  Calculate integral of sin(x) from 0.0 to */
	/*            PI/3.0 using Trapazoidal Method. Result  */
	/*            is 0.5. There are 17 float64 precision    */
	/*            operations per loop (6 +, 2 -, 9 *, 0 /) */
	/*            included in the timing.                  */
	/*            35.3% +, 11.8% -, 52.9% *, and 00.0% /   */
	/*******************************************************/

	//fmt.Printf("starting module 3, m=%ld\n", m );
	x = piref / (three * float64(m)) /*********************/
	s = 0.0                          /*  Loop 4.          */
	v = 0.0                          /*********************/

	dtime(TimeArray)
	for i = 1; i <= m-1; i++ {
		v = v + one
		u = v * x
		w = u * u
		s = s + u*((((((A6*w-A5)*w+A4)*w-A3)*w+A2)*w+A1)*w+one)
	}
	dtime(TimeArray)
	T[9] = T[1]*TimeArray[1] - nulltime

	u = piref / three
	w = u * u
	sa = u * ((((((A6*w-A5)*w+A4)*w-A3)*w+A2)*w+A1)*w + one)

	T[10] = T[9] / 17.0         /*********************/
	sa = x * (sa + two*s) / two /* sin(x) Results.   */
	sb = 0.5                    /*********************/
	sc = sa - sb
	T[11] = one / T[10]
	/*********************/
	/*   DO NOT REMOVE   */
	/*   THIS PRINTOUT!  */
	/*********************/
	fmt.Printf("     3   %13.4e  %10.4f  %10.4f\n", sc, T[9], T[11])

	/************************************************************/
	/* Module 4.  Calculate Integral of cos(x) from 0.0 to PI/3 */
	/*            using the Trapazoidal Method. Result is       */
	/*            sin(PI/3). There are 15 float64 precision      */
	/*            operations per loop (7 +, 0 -, 8 *, and 0 / ) */
	/*            included in the timing.                       */
	/*            50.0% +, 00.0% -, 50.0% *, 00.0% /            */
	/************************************************************/
	A3 = -A3
	A5 = -A5
	x = piref / (three * float64(m)) /*********************/
	s = 0.0                          /*  Loop 5.          */
	v = 0.0                          /*********************/

	dtime(TimeArray)
	for i = 1; i <= m-1; i++ {
		u = float64(i) * x
		w = u * u
		s = s + w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1) + one
	}
	dtime(TimeArray)
	T[12] = T[1]*TimeArray[1] - nulltime

	u = piref / three
	w = u * u
	sa = w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1) + one

	T[13] = T[12] / 15.0              /*******************/
	sa = x * (sa + one + two*s) / two /* Module 4 Result */
	u = piref / three                 /*******************/
	w = u * u
	sb = u * ((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w + A0)
	sc = sa - sb
	T[14] = one / T[13]
	/*********************/
	/*   DO NOT REMOVE   */
	/*   THIS PRINTOUT!  */
	/*********************/
	fmt.Printf("     4   %13.4e  %10.4f  %10.4f\n", sc, T[12], T[14])

	/************************************************************/
	/* Module 5.  Calculate Integral of tan(x) from 0.0 to PI/3 */
	/*            using the Trapazoidal Method. Result is       */
	/*            ln(cos(PI/3)). There are 29 float64 precision  */
	/*            operations per loop (13 +, 0 -, 15 *, and 1 /)*/
	/*            included in the timing.                       */
	/*            46.7% +, 00.0% -, 50.0% *, and 03.3% /        */
	/************************************************************/

	x = piref / (three * float64(m)) /*********************/
	s = 0.0                          /*  Loop 6.          */
	v = 0.0                          /*********************/

	dtime(TimeArray)
	for i = 1; i <= m-1; i++ {
		u = float64(i) * x
		w = u * u
		v = u * ((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w + one)
		s = s + v/(w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one)
	}
	dtime(TimeArray)
	T[15] = T[1]*TimeArray[1] - nulltime

	u = piref / three
	w = u * u
	sa = u * ((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w + one)
	sb = w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1) + one
	sa = sa / sb

	T[16] = T[15] / 29.0        /*******************/
	sa = x * (sa + two*s) / two /* Module 5 Result */
	sb = 0.6931471805599453     /*******************/
	sc = sa - sb
	T[17] = one / T[16]
	/*********************/
	/*   DO NOT REMOVE   */
	/*   THIS PRINTOUT!  */
	/*********************/
	fmt.Printf("     5   %13.4e  %10.4f  %10.4f\n", sc, T[15], T[17])

	/************************************************************/
	/* Module 6.  Calculate Integral of sin(x)*cos(x) from 0.0  */
	/*            to PI/4 using the Trapazoidal Method. Result  */
	/*            is sin(PI/4)^2. There are 29 float64 precision */
	/*            operations per loop (13 +, 0 -, 16 *, and 0 /)*/
	/*            included in the timing.                       */
	/*            46.7% +, 00.0% -, 53.3% *, and 00.0% /        */
	/************************************************************/

	x = piref / (four * float64(m)) /*********************/
	s = 0.0                         /*  Loop 7.          */
	v = 0.0                         /*********************/

	dtime(TimeArray)
	for i = 1; i <= m-1; i++ {
		u = float64(i) * x
		w = u * u
		v = u * ((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w + one)
		s = s + v*(w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one)
	}
	dtime(TimeArray)
	T[18] = T[1]*TimeArray[1] - nulltime

	u = piref / four
	w = u * u
	sa = u * ((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w + one)
	sb = w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1) + one
	sa = sa * sb

	T[19] = T[18] / 29.0        /*******************/
	sa = x * (sa + two*s) / two /* Module 6 Result */
	sb = 0.25                   /*******************/
	sc = sa - sb
	T[20] = one / T[19]
	/*********************/
	/*   DO NOT REMOVE   */
	/*   THIS PRINTOUT!  */
	/*********************/
	fmt.Printf("     6   %13.4e  %10.4f  %10.4f\n", sc, T[18], T[20])

	/*******************************************************/
	/* Module 7.  Calculate value of the definite integral */
	/*            from 0 to sa of 1/(x+1), x/(x*x+1), and  */
	/*            x*x/(x*x*x+1) using the Trapizoidal Rule.*/
	/*            There are 12 float64 precision operations */
	/*            per loop ( 3 +, 3 -, 3 *, and 3 / ) that */
	/*            are included in the timing.              */
	/*            25.0% +, 25.0% -, 25.0% *, and 25.0% /   */
	/*******************************************************/

	/*********************/
	s = 0.0 /* Loop 8.           */
	w = one /*********************/
	sa = 102.3321513995275
	v = sa / float64(m)

	dtime(TimeArray)
	for i = 1; i <= m-1; i++ {
		x = float64(i) * v
		u = x * x
		s = s - w/(x+w) - x/(u+w) - u/(x*u+w)
	}
	dtime(TimeArray)
	T[21] = T[1]*TimeArray[1] - nulltime
	/*********************/
	/* Module 7 Results  */
	/*********************/
	T[22] = T[21] / 12.0
	x = sa
	u = x * x
	sa = -w - w/(x+w) - x/(u+w) - u/(x*u+w)
	sa = 18.0 * v * (sa + two*s)

	m = -2000 * long(sa)
	m = long(float64(m) / scale)

	sc = sa + 500.2
	T[23] = one / T[22]
	/********************/
	/*  DO NOT REMOVE   */
	/*  THIS PRINTOUT!  */
	/********************/
	fmt.Printf("     7   %13.4e  %10.4f  %10.4f\n", sc, T[21], T[23])

	/************************************************************/
	/* Module 8.  Calculate Integral of sin(x)*cos(x)*cos(x)    */
	/*            from 0 to PI/3 using the Trapazoidal Method.  */
	/*            Result is (1-cos(PI/3)^3)/3. There are 30     */
	/*            float64 precision operations per loop included */
	/*            in the timing:                                */
	/*               13 +,     0 -,    17 *          0 /        */
	/*            46.7% +, 00.0% -, 53.3% *, and 00.0% /        */
	/************************************************************/

	x = piref / (three * float64(m)) /*********************/
	s = 0.0                          /*  Loop 9.          */
	v = 0.0                          /*********************/

	dtime(TimeArray)
	for i = 1; i <= m-1; i++ {
		u = float64(i) * x
		w = u * u
		v = w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1) + one
		s = s + v*v*u*((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w+one)
	}
	dtime(TimeArray)
	T[24] = T[1]*TimeArray[1] - nulltime

	u = piref / three
	w = u * u
	sa = u * ((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w + one)
	sb = w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1) + one
	sa = sa * sb * sb

	T[25] = T[24] / 30.0        /*******************/
	sa = x * (sa + two*s) / two /* Module 8 Result */
	sb = 0.29166666666666667    /*******************/
	sc = sa - sb
	T[26] = one / T[25]
	/*********************/
	/*   DO NOT REMOVE   */
	/*   THIS PRINTOUT!  */
	/*********************/
	fmt.Printf("     8   %13.4e  %10.4f  %10.4f\n", sc, T[24], T[26])

	/**************************************************/
	/* MFLOPS(1) output. This is the same weighting   */
	/* used for all previous versions of the flops.c  */
	/* program. Includes Modules 2 and 3 only.        */
	/**************************************************/
	T[27] = (five*(T[6]-T[5]) + T[9]) / 52.0
	T[28] = one / T[27]

	/**************************************************/
	/* MFLOPS(2) output. This output does not include */
	/* Module 2, but it still does 9.2% FDIV's.       */
	/**************************************************/
	T[29] = T[2] + T[9] + T[12] + T[15] + T[18]
	T[29] = (T[29] + four*T[21]) / 152.0
	T[30] = one / T[29]

	/**************************************************/
	/* MFLOPS(3) output. This output does not include */
	/* Module 2, but it still does 3.4% FDIV's.       */
	/**************************************************/
	T[31] = T[2] + T[9] + T[12] + T[15] + T[18]
	T[31] = (T[31] + T[21] + T[24]) / 146.0
	T[32] = one / T[31]

	/**************************************************/
	/* MFLOPS(4) output. This output does not include */
	/* Module 2, and it does NO FDIV's.               */
	/**************************************************/
	T[33] = (T[9] + T[12] + T[18] + T[24]) / 91.0
	T[34] = one / T[33]

	fmt.Printf("\n")
	fmt.Printf("   Iterations      = %10d\n", m)
	fmt.Printf("   NullTime (usec) = %10.4f\n", nulltime)
	fmt.Printf("   MFLOPS(1)       = %10.4f\n", T[28])
	fmt.Printf("   MFLOPS(2)       = %10.4f\n", T[30])
	fmt.Printf("   MFLOPS(3)       = %10.4f\n", T[32])
	fmt.Printf("   MFLOPS(4)       = %10.4f\n\n", T[34])

	//   return 0;
}

/*****************************************************/
/* Various timer routines.                           */
/* Al Aburto, aburto@nosc.mil, 18 Feb 1997           */
/*                                                   */
/* dtime(p) outputs the elapsed time seconds in p[1] */
/* from a call of dtime(p) to the next call of       */
/* dtime(p).  Use CAUTION as some of these routines  */
/* will mess up when timing across the hour mark!!!  */
/*                                                   */
/* For timing I use the 'user' time whenever         */
/* possible. Using 'user+sys' time is a separate     */
/* issue.                                            */
/*                                                   */
/* Example Usage:                                    */
/* [Timer options added here]                        */
/* float64 RunTime, TimeArray[3];                     */
/* main()                                            */
/* {                                                 */
/* dtime(TimeArray);                                 */
/* [routine to time]                                 */
/* dtime(TimeArray);                                 */
/* RunTime = TimeArray[1];                           */
/* }                                                 */
/* [Timer code added here]                           */
/*****************************************************/

/******************************/
/* Timer code.                */
/******************************/

/*****************************************************/
/*  UNIX dtime(). This is the preferred UNIX timer.  */
/*  Provided by: Markku Kolkka, mk59200@cc.tut.fi    */
/*  HP-UX Addition by: Bo Thide', bt@irfu.se         */
/*****************************************************/
/*
#ifdef UNIX
#include <sys/time.h>
#include <sys/resource.h>

#ifdef hpux
#include <sys/syscall.h>
#define getrusage(a,b) syscall(SYS_getrusage,a,b)
#endif

struct rusage rusage;

int dtime(p)
float64 p[];
{
 float64 q;

 q = p[2];

 getrusage(RUSAGE_SELF,&rusage);

 p[2] = (float64)(rusage.ru_utime.tv_sec);
 p[2] = p[2] + (float64)(rusage.ru_utime.tv_usec) * 1.0e-06;
 p[1] = p[2] - q;

 return 0;
}
#endif*/

func dtime(p []float64) {
	var ru syscall.Rusage
	err := syscall.Getrusage(syscall.RUSAGE_SELF, &ru)
	if err != nil {
		fmt.Print(err)
	}
	q := p[2]
	p[2] = float64(ru.Utime.Sec) + (float64(ru.Utime.Usec) * 1.0e-6)
	p[1] = p[2] - q
	//fmt.Printf("TODO: write dtime()\n")
}

/*------ End flops.c code, say good night Jan! (Sep 1992) ------*/
