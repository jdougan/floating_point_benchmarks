
#   Options for compiling and linking D fbench programs

#       For gdmd/GDC
DCC = gdmd
#       For CLANG
#CC = $(HOME)/linuxtools/clang+llvm-5.0.0-linux-x86_64-ubuntu16.04/bin/clang

#   gsms.GDC compiler common options
# this is specifically for gdc
DOPTS = -q,-O3 -inline 

#   Enable to use internal trigonometric functions in fbench
# FIXME what is the D equivalent?
#INTRIG = -DINTRIG

#   Iterations to run in timing tests
ITERATIONS_FBENCH_TEST = 20000000
ITERATIONS_FBENCH_ARCHIVE = 782608690

ITERATIONS_FFBENCH_TEST = 1000

PROGRAMS = fbench ffbench

all:	$(PROGRAMS)

clean:
	rm -f $(PROGRAMS) *.o *.out core* fbench ffbench *.map

allclean: clean
	# del subdirs, onece we have them

fbench: fbench.d math_fbench.d
	$(DCC) fbench.d  math_fbench.d $(DOPTS)

test_fbench: fbench
	echo -e "" | sh -c "time -p ./fbench $(ITERATIONS_FBENCH_TEST)"

time_fbench_archival: fbench
	echo -e "" | sh -c "time -p ./fbench $(ITERATIONS_FBENCH_ARCHIVE) ; time -p ./fbench $(ITERATIONS_FBENCH_ARCHIVE) ; time -p ./fbench $(ITERATIONS_FBENCH_ARCHIVE) "

ffbench: ffbench.d
	$(DCC) ffbench.d -debug -g -gs -map -op -w

test_ffbench: fbench
	echo -e "" | sh -c "time -p ./fbench $(ITERATIONS_FFBENCH_TEST)"
