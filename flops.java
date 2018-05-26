/*
	Trivial applet that displays a string - 4/96 PNL
*/

import java.awt.*;
import java.applet.Applet;
import java.util.Date;

public class flops extends Applet {
	static void maybePrintenv(java.io.PrintWriter out, String key) {
		String v = System.getProperty(key);
		if (v != null) {
			out.print(key);
			out.print(": ");
			out.print(v);
			out.println();
		}
	}

	TextArea ta = new TextArea();

	public void init() {
		add(ta);
		ta.setVisible(true);
		java.io.StringWriter os = new java.io.StringWriter();
		java.io.PrintWriter pw = new java.io.PrintWriter(os);
		test(pw, 2.0);
		ta.setText(os.toString());
	}
	
	
	public static void main( String[] argv ) {
		System.out.print("   FLOPS Java Program (Double Precision), V2.0 18 Dec 1992\n\n");
		System.out.print("   Module     Error        RunTime      MFLOPS\n");
		System.out.print("                            (usec)\n");
		test(new java.io.PrintWriter(new java.io.OutputStreamWriter(System.out)));
	}

	static void test(java.io.PrintWriter out) {
		test(out, 15.0);
	}
	static void test(java.io.PrintWriter out, double TLimit) {
		maybePrintenv(out, "java.vm.name");
		maybePrintenv(out, "java.vm.vendor");
		maybePrintenv(out, "java.vm.version");
		double nulltime, TimeArray[];   /* Variables needed for 'dtime()'.     */
		//double TLimit = 15;                   /* Threshold to determine Number of    */
						 /* Loops to run. Fixed at 15.0 seconds.*/

		double T[];                    /* Global Array used to hold timing    */
		double sa,sb,sc,sd,one = 1,two = 2,three = 3;
		double four = 4,five = 5,piref = 3.14159265358979324,piprg;
		double scale,pierr;

		double A0 = 1.0;
		double A1 = -0.1666666666671334;
		double A2 = 0.833333333809067E-2;
		double A3 = 0.198412715551283E-3;
		double A4 = 0.27557589750762E-5;
		double A5 = 0.2507059876207E-7;
		double A6 = 0.164105986683E-9;

		double B0 = 1.0;
		double B1 = -0.4999999999982;
		double B2 = 0.4166666664651E-1;
		double B3 = -0.1388888805755E-2;
		double B4 = 0.24801428034E-4;
		double B5 = -0.2754213324E-6;
		double B6 = 0.20189405E-8;

		double C0 = 1.0;
		double C1 = 0.99999999668;
		double C2 = 0.49999995173;
		double C3 = 0.16666704243;
		double C4 = 0.4166685027E-1;
		double C5 = 0.832672635E-2;
		double C6 = 0.140836136E-2;
		double C7 = 0.17358267E-3;
		double C8 = 0.3931683E-4;

		double D1 = 0.3999999946405E-1;
		double D2 = 0.96E-3;
		double D3 = 0.1233153E-5;

		double E2 = 0.48E-3;
		double E3 = 0.411051E-6;

		double s = 0, u,v,w, x = 0;
		long loops = 15625, NLimit = 512000000;
		long i, m, n;
		TimeArray = new double[3];
		T = new double[36];

		/* Initialize the timer. */

		dtime(TimeArray);
		dtime(TimeArray);
		scale = one;
   		
		T[1] = 1.0E+06/(double)loops;

		/* Module 1.  Calculate integral of df(x)/f(x) defined */
		/*            below.  Result is ln(f(1)). There are 14 */
		/*            double precision operations per loop     */
		/*            ( 7 +, 0 -, 6 *, 1 / ) that are included */
		/*            in the timing.                           */
		/*            50.0% +, 00.0% -, 42.9% *, and 07.1% /   */
		n = loops;
		sa = 0.0;

		while ( sa < TLimit ) {
			n = 2 * n;
			x = one / (double)n;
			s = 0.0;                                        /*  Loop 1.          */
			v = 0.0;
			w = one;

			dtime(TimeArray);
			for( i = 1 ; i <= n-1 ; i++ ) {
				v = v + w;
				u = v * x;
				s = s + (D1+u*(D2+u*D3))/(w+u*(D1+u*(E2+u*E3)));
			}
			dtime(TimeArray);
			sa = TimeArray[1];

			if ( n == NLimit ) break;
			/* printf(" %10ld  %12.5lf\n",n,sa); */
		}

		scale = 1.0E+06 / (double)n;
		T[1]  = scale;

		/* Estimate nulltime ('for' loop time). */
		dtime(TimeArray);
		for( i = 1 ; i <= n-1 ; i++ ) {
		}
		dtime(TimeArray);
		nulltime = T[1] * TimeArray[1];
		if ( nulltime < 0.0 ) nulltime = 0.0;

		T[2] = T[1] * sa - nulltime;

		sa = (D1+D2+D3)/(one+D1+E2+E3);
		sb = D1;

		T[3] = T[2] / 14.0;
		sa = x * ( sa + sb + two * s ) / two;           /* Module 1 Results. */
		sb = one / sa;
		n  = (long)( (double)( 40000 * (long)sb ) / scale );
		sc = sb - 25.2;
		T[4] = one / T[3];
		out.print ("     1   " + sc + "  " + T[2] + "  " + T[4] + '\n');
		m = n;

		/* Module 2.  Calculate value of PI from Taylor Series */
		/*            expansion of atan(1.0).  There are 7     */
		/*            double precision operations per loop     */
		/*            ( 3 +, 2 -, 1 *, 1 / ) that are included */
		/*            in the timing.                           */
		/*            42.9% +, 28.6% -, 14.3% *, and 14.3% /   */

		s  = -five;
		sa = -one;                                       /* Loop 2.          */
		dtime(TimeArray);
		for ( i = 1 ; i <= m ; i++ ) {
			s  = -s;
			sa = sa + s;
		}
		dtime(TimeArray);
		T[5] = T[1] * TimeArray[1];
		if ( T[5] < 0.0 ) T[5] = 0.0;

		sc   = (double)m;

		u = sa;
		v = 0.0;                                        /* Loop 3.           */
		w = 0.0;
		x = 0.0;

		dtime(TimeArray);
		for ( i = 1 ; i <= m ; i++) {
			s  = -s;
			sa = sa + s;
			u  = u + two;
			x  = x +(s - u);
			v  = v - s * u;
			w  = w + s / u;
		}
		dtime(TimeArray);
		T[6] = T[1] * TimeArray[1];

		T[7] = ( T[6] - T[5] ) / 7.0;
		m  = (long)( sa * x  / sc );                    /*  PI Results       */
		sa = four * w / five;
		sb = sa + five / v;
		sc = 31.25;
		piprg = sb - sc / (v * v * v);
		pierr = piprg - piref;
		T[8]  = one  / T[7];
		out.print ("     2   " + pierr + "  " + (T[6]-T[5]) + "  " + T[8] + '\n');

		/* Module 3.  Calculate integral of sin(x) from 0.0 to */
		/*            PI/3.0 using Trapazoidal Method. Result  */
		/*            is 0.5. There are 17 double precision    */
		/*            operations per loop (6 +, 2 -, 9 *, 0 /) */
		/*            included in the timing.                  */
		/*            35.3% +, 11.8% -, 52.9% *, and 00.0% /   */
		/*******************************************************/

		x = piref / ( three * (double)m );              /*********************/
		s = 0.0;                                        /*  Loop 4.          */
		v = 0.0;                                        /*********************/

		dtime(TimeArray);
		for( i = 1 ; i <= m-1 ; i++ ) {
			v = v + one;
			u = v * x;
			w = u * u;
			s = s + u * ((((((A6*w-A5)*w+A4)*w-A3)*w+A2)*w+A1)*w+one);
		}
		dtime(TimeArray);
		T[9]  = T[1] * TimeArray[1] - nulltime;

		u  = piref / three;
		w  = u * u;
		sa = u * ((((((A6*w-A5)*w+A4)*w-A3)*w+A2)*w+A1)*w+one);

		T[10] = T[9] / 17.0;                            /*********************/
		sa = x * ( sa + two * s ) / two;                /* sin(x) Results.   */
		sb = 0.5;                                       /*********************/
		sc = sa - sb;
		T[11] = one / T[10];
							  /*********************/
							  /*   DO NOT REMOVE   */
							  /*   THIS PRINTOUT!  */
							  /*********************/
		out.print("     3   " + sc + "  " + T[9] + "  " + T[11] + "\n");

		/************************************************************/
		/* Module 4.  Calculate Integral of cos(x) from 0.0 to PI/3 */
		/*            using the Trapazoidal Method. Result is       */
		/*            sin(PI/3). There are 15 double precision      */
		/*            operations per loop (7 +, 0 -, 8 *, and 0 / ) */
		/*            included in the timing.                       */
		/*            50.0% +, 00.0% -, 50.0% *, 00.0% /            */
		/************************************************************/
		A3 = -A3;
		A5 = -A5;
		x = piref / ( three * (double)m );              /*********************/
		s = 0.0;                                        /*  Loop 5.          */
		v = 0.0;                                        /*********************/

		dtime(TimeArray);
		for( i = 1 ; i <= m-1 ; i++ ) {
			u = (double)i * x;
			w = u * u;
			s = s + w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one;
		}
		dtime(TimeArray);
		T[12]  = T[1] * TimeArray[1] - nulltime;

		u  = piref / three;
		w  = u * u;
		sa = w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one;

		T[13] = T[12] / 15.0;                             /*******************/
		sa = x * ( sa + one + two * s ) / two;            /* Module 4 Result */
		u  = piref / three;                               /*******************/
		w  = u * u;
		sb = u * ((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w+A0);
		sc = sa - sb;
		T[14] = one / T[13];
							  /*********************/
							  /*   DO NOT REMOVE   */
							  /*   THIS PRINTOUT!  */
							  /*********************/
		out.print("     4   " + sc + "  " + T[12] + "  " + T[14] + "\n");

		/************************************************************/
		/* Module 5.  Calculate Integral of tan(x) from 0.0 to PI/3 */
		/*            using the Trapazoidal Method. Result is       */
		/*            ln(cos(PI/3)). There are 29 double precision  */
		/*            operations per loop (13 +, 0 -, 15 *, and 1 /)*/
		/*            included in the timing.                       */
		/*            46.7% +, 00.0% -, 50.0% *, and 03.3% /        */
		/************************************************************/

		x = piref / ( three * (double)m );              /*********************/
		s = 0.0;                                        /*  Loop 6.          */
		v = 0.0;                                        /*********************/

		dtime(TimeArray);
		for( i = 1 ; i <= m-1 ; i++ )
		{
		u = (double)i * x;
		w = u * u;
		v = u * ((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w+one);
		s = s + v / (w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one);
		}
		dtime(TimeArray);
		T[15]  = T[1] * TimeArray[1] - nulltime;

		u  = piref / three;
		w  = u * u;
		sa = u*((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w+one);
		sb = w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one;
		sa = sa / sb;

		T[16] = T[15] / 29.0;                             /*******************/
		sa = x * ( sa + two * s ) / two;                  /* Module 5 Result */
		sb = 0.6931471805599453;                          /*******************/
		sc = sa - sb;
		T[17] = one / T[16];
							  /*********************/
							  /*   DO NOT REMOVE   */
							  /*   THIS PRINTOUT!  */
							  /*********************/
		out.print("     5   " + sc + "  " + T[15] + "  " + T[17] + "\n");

		/************************************************************/
		/* Module 6.  Calculate Integral of sin(x)*cos(x) from 0.0  */
		/*            to PI/4 using the Trapazoidal Method. Result  */
		/*            is sin(PI/4)^2. There are 29 double precision */
		/*            operations per loop (13 +, 0 -, 16 *, and 0 /)*/
		/*            included in the timing.                       */
		/*            46.7% +, 00.0% -, 53.3% *, and 00.0% /        */
		/************************************************************/

		x = piref / ( four * (double)m );               /*********************/
		s = 0.0;                                        /*  Loop 7.          */
		v = 0.0;                                        /*********************/

		dtime(TimeArray);
		for( i = 1 ; i <= m-1 ; i++ )
		{
		u = (double)i * x;
		w = u * u;
		v = u * ((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w+one);
		s = s + v*(w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one);
		}
		dtime(TimeArray);
		T[18]  = T[1] * TimeArray[1] - nulltime;

		u  = piref / four;
		w  = u * u;
		sa = u*((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w+one);
		sb = w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one;
		sa = sa * sb;

		T[19] = T[18] / 29.0;                             /*******************/
		sa = x * ( sa + two * s ) / two;                  /* Module 6 Result */
		sb = 0.25;                                        /*******************/
		sc = sa - sb;
		T[20] = one / T[19];
							  /*********************/
							  /*   DO NOT REMOVE   */
							  /*   THIS PRINTOUT!  */
							  /*********************/
		out.print("     6   " + sc + "  " + T[18] + "  " + T[20] + "\n");


		/*******************************************************/
		/* Module 7.  Calculate value of the definite integral */
		/*            from 0 to sa of 1/(x+1), x/(x*x+1), and  */
		/*            x*x/(x*x*x+1) using the Trapizoidal Rule.*/
		/*            There are 12 double precision operations */
		/*            per loop ( 3 +, 3 -, 3 *, and 3 / ) that */
		/*            are included in the timing.              */
		/*            25.0% +, 25.0% -, 25.0% *, and 25.0% /   */
		/*******************************************************/

							  /*********************/
		s = 0.0;                                        /* Loop 8.           */
		w = one;                                        /*********************/
		sa = 102.3321513995275;
		v = sa / (double)m;

		dtime(TimeArray);
		for ( i = 1 ; i <= m-1 ; i++)
		{
		x = (double)i * v;
		u = x * x;
		s = s - w / ( x + w ) - x / ( u + w ) - u / ( x * u + w );
		}
		dtime(TimeArray);
		T[21] = T[1] * TimeArray[1] - nulltime;
							  /*********************/
							  /* Module 7 Results  */
							  /*********************/
		T[22] = T[21] / 12.0;                                  
		x  = sa;                                      
		u  = x * x;
		sa = -w - w / ( x + w ) - x / ( u + w ) - u / ( x * u + w );
		sa = 18.0 * v * (sa + two * s );

		m  = -2000 * (long)sa;
		m = (long)( (double)m / scale );

		sc = sa + 500.2;
		T[23] = one / T[22];
							  /********************/
							  /*  DO NOT REMOVE   */
							  /*  THIS PRINTOUT!  */
							  /********************/
		out.print("     7   " + sc + "  " + T[21] + "  " + T[23] + "\n");

		/************************************************************/
		/* Module 8.  Calculate Integral of sin(x)*cos(x)*cos(x)    */
		/*            from 0 to PI/3 using the Trapazoidal Method.  */
		/*            Result is (1-cos(PI/3)^3)/3. There are 30     */
		/*            double precision operations per loop included */
		/*            in the timing:                                */
		/*               13 +,     0 -,    17 *          0 /        */
		/*            46.7% +, 00.0% -, 53.3% *, and 00.0% /        */
		/************************************************************/

		x = piref / ( three * (double)m );              /*********************/
		s = 0.0;                                        /*  Loop 9.          */
		v = 0.0;                                        /*********************/

		dtime(TimeArray);
		for( i = 1 ; i <= m-1 ; i++ )
		{
		u = (double)i * x;
		w = u * u;
		v = w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one;
		s = s + v*v*u*((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w+one);
		}
		dtime(TimeArray);
		T[24]  = T[1] * TimeArray[1] - nulltime;

		u  = piref / three;
		w  = u * u;
		sa = u*((((((A6*w+A5)*w+A4)*w+A3)*w+A2)*w+A1)*w+one);
		sb = w*(w*(w*(w*(w*(B6*w+B5)+B4)+B3)+B2)+B1)+one;
		sa = sa * sb * sb;

		T[25] = T[24] / 30.0;                             /*******************/
		sa = x * ( sa + two * s ) / two;                  /* Module 8 Result */
		sb = 0.29166666666666667;                         /*******************/
		sc = sa - sb;
		T[26] = one / T[25];
							  /*********************/
							  /*   DO NOT REMOVE   */
							  /*   THIS PRINTOUT!  */
							  /*********************/
		out.print("     8   " + sc + "  " + T[24] + "  " + T[26] + "\n");

		/**************************************************/   
		/* MFLOPS(1) output. This is the same weighting   */
		/* used for all previous versions of the flops.c  */
		/* program. Includes Modules 2 and 3 only.        */
		/**************************************************/ 
		T[27] = ( five * (T[6] - T[5]) + T[9] ) / 52.0;
		T[28] = one  / T[27];

		/**************************************************/   
		/* MFLOPS(2) output. This output does not include */
		/* Module 2, but it still does 9.2% FDIV's.       */
		/**************************************************/ 
		T[29] = T[2] + T[9] + T[12] + T[15] + T[18];
		T[29] = (T[29] + four * T[21]) / 152.0;
		T[30] = one / T[29];

		/**************************************************/   
		/* MFLOPS(3) output. This output does not include */
		/* Module 2, but it still does 3.4% FDIV's.       */
		/**************************************************/ 
		T[31] = T[2] + T[9] + T[12] + T[15] + T[18];
		T[31] = (T[31] + T[21] + T[24]) / 146.0;
		T[32] = one / T[31];

		/**************************************************/   
		/* MFLOPS(4) output. This output does not include */
		/* Module 2, and it does NO FDIV's.               */
		/**************************************************/ 
		T[33] = (T[9] + T[12] + T[18] + T[24]) / 91.0;
		T[34] = one / T[33];


		out.print("\n");
		out.print("   Iterations      = " + m + "\n");
		out.print("   NullTime (usec) = " + nulltime + "\n");
		out.print("   MFLOPS(1)       = " + T[28] + "\n");
		out.print("   MFLOPS(2)       = " + T[30] + "\n");
		out.print("   MFLOPS(3)       = " + T[32] + "\n");
		out.print("   MFLOPS(4)       = " + T[34] + "\n\n");
		out.flush();
	}
	static void dtime( double p[] ) {
		double q = p[2];
		Date td = new Date();
		p[2] = td.getTime() / 1000.0;//the new time
		p[1] = p[2] - q;
	}
}
