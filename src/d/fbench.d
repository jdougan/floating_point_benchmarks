@safe 
/*

        John Walker's Floating Point Benchmark, derived from...

        Marinchip Interactive Lens Design System

                                     John Walker   December 1980

                            By John Walker
                       http://www.fourmilab.ch/

        This  program may be used, distributed, and modified freely as
        long as the origin information is preserved.

        This  is  a  complete  optical  design  raytracing  algorithm,
        stripped of its user interface and recast into portable C.  It
        not only determines execution speed on an  extremely  floating
        point   (including   trig   function)   intensive   real-world
        application, it  checks  accuracy  on  an  algorithm  that  is
        exquisitely  sensitive  to  errors.   The  performance of this
        program is typically far more  sensitive  to  changes  in  the
        efficiency  of  the  trigonometric  library  routines than the
        average floating point program.

        The benchmark may be compiled in two  modes.   If  the  symbol
        INTRIG  is  defined,  built-in  trigonometric  and square root
        routines will be used for all calculations.  Timings made with
        INTRIG  defined  reflect  the  machine's  basic floating point
        performance for the arithmetic operators.  If  INTRIG  is  not
        defined,  the  system  library  <math.h>  functions  are used.
        Results with INTRIG not defined reflect the  system's  library
        performance  and/or  floating  point hardware support for trig
        functions and square root.  Results with INTRIG defined are  a
        good  guide  to  general  floating  point  performance,  while
        results with INTRIG undefined indicate the performance  of  an
        application which is math function intensive.

        Special  note  regarding  errors in accuracy: this program has
        generated numbers identical to the last digit it  formats  and
        checks on the following machines, floating point
        architectures, and languages:

        Marinchip 9900    QBASIC    IBM 370 double-precision (REAL * 8) format

        IBM PC / XT / AT  Lattice C IEEE 64 bit, 80 bit temporaries
                          High C    same, in line 80x87 code
                          BASICA    "Double precision"
                          Quick BASIC IEEE double precision, software routines

        Sun 3             C         IEEE 64 bit, 80 bit temporaries,
                                    in-line 68881 code, in-line FPA code.

        MicroVAX II       C         Vax "G" format floating point

        Macintosh Plus    MPW C     SANE floating point, IEEE 64 bit format
                                    implemented in ROM.

        Inaccuracies  reported  by  this  program should be taken VERY
        SERIOUSLY INDEED, as the program has been demonstrated  to  be
        invariant  under  changes in floating point format, as long as
        the format is a recognised double precision  format.   If  you
        encounter errors, please remember that they are just as likely
        to  be  in  the  floating  point  editing   library   or   the
        trigonometric  libraries  as  in  the low level operator code.

        The benchmark assumes that results are basically reliable, and
        only tests the last result computed against the reference.  If
        you're running on  a  suspect  system  you  can  compile  this
        program  with  ACCURACY defined.  This will generate a version
        which executes as an infinite loop, performing the  ray  trace
        and checking the results on every pass.  All incorrect results
        will be reported.
*/

// We are timing manually, read input from the user
// to prompt them to hit the stopwatch.
// maybe this should be an option, but does anyone need to use a
// stopwatch on modern systems? default to false.
enum manualTiming = false;

// Should we print the calculation results?
// Essentially a replace ment for the ACCRACY option.
enum showResults = false;

enum internalMath = false;
static if(internalMath)
  {
    // use a cut down transcendantals package
    // this will test generap fp spped rather than\
    // transcendantal fns.
    import math_fbench;
  }
else 
  {
    // system library
    import core.stdc.math;
  }

// #define cot(x) (1.0 / tan(x));

// for C APIs
enum TRUE   = 1;
enum FALSE = 0;

enum max_surfaces = 10;

/*  Local variables  */

// buffer for reading the "press any key" input 
static char[132] tbfr;

static short current_surfaces;
static short paraxial;

static double clear_aperture;

static double aberr_lspher;
static double aberr_osc;
static double aberr_lchrom;

static double max_lspher;
static double max_osc;
static double max_lchrom;

static double radius_of_curvature;
static double object_distance;
static double ray_height;
static double axis_slope_angle;
static double from_index;
static double to_index;

static double[9] spectral_line;
static double[max_surfaces][5] s;
static double[2][2] od_sa;

/* Computed output of program goes here */
static char[80][8] outarr;         

int itercount;                     /* The iteration counter for the main loop
                                      in the program is made global so that
                                      the compiler should not be allowed to
                                      optimise out the loop over the ray
                                      tracing code. */

// takes aboout 5 seconds on my laptop.
enum ITERATIONS = 10000000;

int niter = ITERATIONS;            /* Iteration counter */

// Apparantly cotangent isn't standard in system math libraries?!
// Original is a macro, change to mixin
@nogc nothrow pure @safe double cot(double x)
{
    pragma(inline,true); 
    return (1.0 / tan(x)); 
}
template cot_macro(string x) {
    const char[] cot_macro = "(1.0 / tan(" ~ x ~ "))";
}

immutable static string[8] refarr = [
    /* Reference results.  These happen to
    be derived from a run on Microsoft
    Quick BASIC on the IBM PC/AT. */

    "   Marginal ray          47.09479120920   0.04178472683",
    "   Paraxial ray          47.08372160249   0.04177864821",
    "Longitudinal spherical aberration:        -0.01106960671",
    "    (Maximum permissible):                 0.05306749907",
    "Offense against sine condition (coma):     0.00008954761",
    "    (Maximum permissible):                 0.00250000000",
    "Axial chromatic aberration:                0.00448229032",
    "    (Maximum permissible):                 0.05306749907"
];

/* The  test  case  used  in  this program is the  design for a 4 inch
   achromatic telescope  objective  used  as  the  example  in  Wyld's
   classic  work  on  ray  tracing by hand, given in Amateur Telescope
   Making, Volume 3.  */
immutable static double[4][4] testcase = [
    [27.05, 1.5137, 63.6, 0.52],
    [-16.68, 1, 0, 0.138],
    [-16.68, 1.6164, 36.7, 0.38],
    [-78.1, 1, 0, 0]
];



/*
    Calculate passage through surface

    If  the variable PARAXIAL is true, the trace through the
    surface will be done using the paraxial  approximations.
    Otherwise,  the normal trigonometric trace will be done.

    This routine takes the following inputs:

    RADIUS_OF_CURVATURE         Radius of curvature of surface
                                being crossed.  If 0, surface is
                                plane.

    OBJECT_DISTANCE             Distance of object focus from
                                lens vertex.  If 0, incoming
                                rays are parallel and
                                the following must be specified:

    RAY_HEIGHT                  Height of ray from axis.  Only
                                relevant if OBJECT.DISTANCE == 0

    AXIS_SLOPE_ANGLE            Angle incoming ray makes with axis
                                at intercept

    FROM_INDEX                  Refractive index of medium being left

    TO_INDEX                    Refractive index of medium being
                                entered.

    The outputs are the following variables:

    OBJECT_DISTANCE             Distance from vertex to object focus
                                after refraction.

    AXIS_SLOPE_ANGLE            Angle incoming ray makes with axis
                                at intercept after refraction.

*/
@nogc nothrow @safe static void transit_surface() 
{
    double
            iang,               /* Incidence angle */
            rang,               /* Refraction angle */
            iang_sin,           /* Incidence angle sin */
            rang_sin,           /* Refraction angle sin */
            old_axis_slope_angle, sagitta;

    if (paraxial) {
        if (radius_of_curvature != 0.0) {
            if (object_distance == 0.0) {
                axis_slope_angle = 0.0;
                iang_sin = ray_height / radius_of_curvature;
            } else {
                iang_sin = ((object_distance -
                    radius_of_curvature) / radius_of_curvature) *
                    axis_slope_angle;
            }
            rang_sin = (from_index / to_index) * iang_sin;
            old_axis_slope_angle = axis_slope_angle;
            axis_slope_angle = axis_slope_angle + iang_sin - rang_sin;
            if (object_distance != 0.0) {
                ray_height = object_distance * old_axis_slope_angle;
            }
            object_distance = ray_height / axis_slope_angle;
            return;
        }
        object_distance = object_distance * (to_index / from_index);
        axis_slope_angle = axis_slope_angle * (from_index / to_index);
        return;
    }  

    if (radius_of_curvature != 0.0) {
        if (object_distance == 0.0) {
            axis_slope_angle = 0.0;
            iang_sin = ray_height / radius_of_curvature;
        } else {
            iang_sin = ((object_distance -
                radius_of_curvature) / radius_of_curvature) *
                sin(axis_slope_angle);
        }
        iang = asin(iang_sin);
        rang_sin = (from_index / to_index) * iang_sin;
        old_axis_slope_angle = axis_slope_angle;
        axis_slope_angle = axis_slope_angle +
            iang - asin(rang_sin);
        sagitta = sin((old_axis_slope_angle + iang) / 2.0);
        sagitta = 2.0 * radius_of_curvature*sagitta*sagitta;
        object_distance = ((radius_of_curvature * sin(
            old_axis_slope_angle + iang)) *
            //cot(axis_slope_angle)) + sagitta;
            mixin(cot_macro!("(axis_slope_angle)")) )  + sagitta;
        return;
    }

    rang = -asin((from_index / to_index) * sin(axis_slope_angle));
    object_distance = object_distance * ((to_index *
        cos(-rang)) / (from_index *
        cos(axis_slope_angle)));
    axis_slope_angle = -rang;
}

/*  Perform ray trace in specific spectral line  */
@nogc nothrow @safe static void trace_line(int line, double ray_h)
{
    int i;

    object_distance = 0.0;
    ray_height = ray_h;
    from_index = 1.0;

    for (i = 1; i <= current_surfaces; i++) {
        radius_of_curvature = s[i][1];
        to_index = s[i][2];
        if (to_index > 1.0) {
            to_index = to_index + ((spectral_line[4] -
                spectral_line[line]) /
                (spectral_line[3] - spectral_line[6])) *
                ((s[i][2] - 1.0) / s[i][3]);
        }
        transit_surface();
        from_index = to_index;
        if (i < current_surfaces) {
            object_distance = object_distance - s[i][4];
        }
    }
}

/*  Initialise when called the first time  */
int main(string[] args)
{
    // FIXME change to an appropriate conversion function in stdc?
    import std.conv: to ; 
    import std.string: toStringz ;
    import core.stdc.stdio;
    import core.stdc.stdlib;
    import core.stdc.string;


    int i, j, errors;
    ulong k;
    double od_fline, od_cline;

    spectral_line[1] = 7621.0;       /* A */
    spectral_line[2] = 6869.955;     /* B */
    spectral_line[3] = 6562.816;     /* C */
    spectral_line[4] = 5895.944;     /* D */
    spectral_line[5] = 5269.557;     /* E */
    spectral_line[6] = 4861.344;     /* F */
    spectral_line[7] = 4340.477;     /* G'*/
    spectral_line[8] = 3968.494;     /* H */

    /* Process the number of iterations argument, if one is supplied. */
    // printf("\n%s\n", toStringz(to!string(args)));
    if (args.length > 1) {
        // FIXME Blow up is there is an arg 1 and it is non-numeric
        niter = to!int(args[1]);
        if (args[1][1] == '-' || niter < 1) {
            printf("This is John Walker's floating point accuracy and\n");
            printf("performance benchmark program.  You call it with\n");
            printf("\nfbench <itercount>\n\n");
            printf("where <itercount> is the number of iterations\n");
            printf("to be executed.  Archival timings should be made\n");
            printf("with the iteration count set so that roughly five\n");
            printf("minutes of execution is timed.\n");
            return(0);
        }
    }

        /* Load test case into working array */

    clear_aperture = 4.0;
    current_surfaces = 4;
    for (i = 0; i < current_surfaces; i++) {
        for (j = 0; j < 4; j++) {
            s[i + 1][j + 1] = testcase[i][j];
        }
    }

    printf("Ready to begin John Walker's floating point accuracy\n");
    printf("and performance benchmark.  %d iterations will be made.\n\n",
       niter);

    printf("Measured run time in seconds should be divided by %.f\n", niter / 1000.0);
    printf("to normalise for reporting results.  For archival results,\n");
    printf("adjust iteration count so the benchmark runs about five minutes.\n\n");

    static if (manualTiming) {
        printf("Start the timer then press return to begin benchmark:");
        fgets(cast(char*)&tbfr, tbfr.length, core.stdc.stdio.stdin);
    }


    for (itercount = 0; itercount < niter; itercount++) {
        for (paraxial = 0; paraxial <= 1; paraxial++) {
            /* Do main trace in D light */
            trace_line(4, clear_aperture / 2.0);
            od_sa[paraxial][0] = object_distance;
            od_sa[paraxial][1] = axis_slope_angle;
        }
        paraxial = FALSE;

        /* Trace marginal ray in C */
        trace_line(3, clear_aperture / 2.0);
        od_cline = object_distance;

        /* Trace marginal ray in F */
        trace_line(6, clear_aperture / 2.0);
        od_fline = object_distance;

        aberr_lspher = od_sa[1][0] - od_sa[0][0];
        aberr_osc = 1.0 - (od_sa[1][0] * od_sa[1][1]) /
          (sin(od_sa[0][1]) * od_sa[0][0]);
        aberr_lchrom = od_fline - od_cline;
        max_lspher = sin(od_sa[0][1]);

        /* D light */
        max_lspher = 0.0000926 / (max_lspher * max_lspher);
        max_osc = 0.0025;
        max_lchrom = max_lspher;
    }

    static if (manualTiming) {
        printf("Stop the timer then press return:\x07");
        fgets(cast(char*)&tbfr, tbfr.length, core.stdc.stdio.stdin);
    }

    /* Now evaluate the accuracy of the results from the last ray trace */
    sprintf(cast(char*)&(outarr[0]), "%15s   %21.11f  %14.11f",
      "Marginal ray".ptr, od_sa[0][0], od_sa[0][1]);
    sprintf(cast(char*)&(outarr[1]), "%15s   %21.11f  %14.11f",
       cast(char*)"Paraxial ray", od_sa[1][0], od_sa[1][1]);
    sprintf(cast(char*)&(outarr[2]),
       "Longitudinal spherical aberration:      %16.11f",
       aberr_lspher);
    sprintf(cast(char*)&(outarr[3]),
       "    (Maximum permissible):              %16.11f",
       max_lspher);
    sprintf(cast(char*)&(outarr[4]),
       "Offense against sine condition (coma):  %16.11f",
       aberr_osc);
    sprintf(cast(char*)&(outarr[5]),
       "    (Maximum permissible):              %16.11f",
       max_osc);
    sprintf(cast(char*)&(outarr[6]),
       "Axial chromatic aberration:             %16.11f",
       aberr_lchrom);
    sprintf(cast(char*)&(outarr[7]),
       "    (Maximum permissible):              %16.11f",
       max_lchrom);

    static if (showResults) {
        //printf("\n");
        for (i = 0; i < 8; i++) {
            printf("%s\n", cast(char*)&(outarr[i]));
        }
    }

    /* Now compare the edited results with the master values from
           reference executions of this program. */
    errors = 0;
    for (i = 0; i < 8; i++) {
        if (strcmp(cast(char*)&(outarr[i]), (refarr[i]).ptr) != 0) {
            printf("\nError in results on line %d...\n", i + 1);
            printf("Expected:  \"%s\"\n", (refarr[i]).ptr);
            printf("Received:  \"%s\"\n", cast(char*)&(outarr[i]));
            printf("(Errors)    ");
            k = strlen(cast(char*)&refarr[i]);
            for (j = 0; j < k; j++) {
                printf("%c", refarr[i][j] == outarr[i][j] ? ' ' : '^');
                if (refarr[i][j] != outarr[i][j]) {
                    errors++;
                }
            }
            printf("\n");
        }
    }
    if (errors > 0) {
        char* suffix  = cast(char*)"s";
        if (errors >1) {
            suffix = cast(char*)"";
        }
        printf("\n%d error%s in results.  This is VERY SERIOUS.\n",
            errors, suffix);
    } else {
        printf("\nNo errors in results.\n");
    }
    return 0;
}
