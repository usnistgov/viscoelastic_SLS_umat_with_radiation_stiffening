---------- INSTALL GNU COMPILERS AND LIBRARIES ----------

Make sure to install the latest GNU Fortran and C/C++ compilers:

sudo apt install build-essential
sudo apt install gfortran
sudo apt install libopenblas-dev

Note, I am using gcc/gfortran version 9.4.0. If you also have Intel
compilers installed, then it would be a good idea to NOT activate their
environment variables. If you have some sort of command in your 
.bashrc file like "source /opt/intel/oneapi/setvars.sh", then comment
it out and restart your terminal session. This just adds an extra
level of security to ensure that the correct compilers are chosen
in the steps to come.


---------- COMPILE ARPACK ----------

ARPACK is a dependent library of CalculiX. At the time of writing, however, this
library was no longer being hosted by Rice University. An old copy of the library
can work, but I also found a more up-to-date version on GitHub called ARPACK-NG.

Download the repository:
https://github.com/opencollab/arpack-ng

The easiest method of compiling ARPACK-NG is through CMake, so make sure that
is installed first. There is an 32-bit version (LP64) and 64-bit (ILP64) version
of ARPACK-NG; we want to use the default LP64 version for CalculiX.

In the high-level directory of ARPACK-NG that contains the "CMakeLists.txt" file,
create a new directory called "build" and then go into this directory, via the
following commands:

mkdir build
cd build

Then, type of the following CMake command, assuming default installation paths:

cmake -D CMAKE_Fortran_COMPILER=/usr/bin/gfortran -D CMAKE_C_COMPILER=/usr/bin/gcc -D EXAMPLES=OFF -D MPI=OFF -D BUILD_SHARED_LIBS=OFF ../

If successful, the Configuration Summary should contain only paths to GNU
compilers and libraries, like this:

-- Configuration summary for arpack-ng-3.9.0:
   -- prefix: /usr/local
   -- MPI: OFF (ICB provided )
   -- ICB: OFF
   -- INTERFACE64: 0
   -- FC:      /usr/bin/gfortran
   -- FCFLAGS: -O3 -DNDEBUG -O3
   -- BLAS:
      -- link:    /usr/lib/x86_64-linux-gnu/libopenblas.so
   -- LAPACK:
      -- link:    /usr/lib/x86_64-linux-gnu/libopenblas.so
      -- link:    -lpthread
      -- link:    -lm
      -- link:    -ldl
-- Configuring done
-- Generating done
-- Build files have been written to: /home/nhm/calculix_build/ccx_gnu/arpack-ng-master/build

Now that CMake has created an appropriate Make file, it is time to actually compile
the library via the "make" command. While in the "build" directory containing the
newly generated "Makefile", simply type,

make

If successful, a new static library should be created, called "libarpack.a". To
make things clear, rename this file to "libarpack_GNU.a". You can also move it 
to a new directory if desired; the filepath will need to be recorded for later.



---------- COMPILE SPOOLES ----------

Download the Spooles 2.2 library:

https://netlib.org/linalg/spooles/spooles.2.2.html

Specifically, grab the entire repository called "spooles.2.2.tgz". On the
CalculiX website, there is also a correction patch to handle large input decks.
This file can be downloaded from here:

http://www.dhondt.de/ccx_2.20.SPOOLEScorrection.tar.bz2

Specifically, this replaces one file in Spooles located in the directory, 
"I2Ohash/src/PaxHeaders.2176/util.c"

By default, Spooles does not compile for multi-threaded applications, so a
couple of edits must be made in the "makefile":

1) Towards the end of the "lib :" target, uncomment the line,
"#cd MT               ; make lib". To do this, delete "#" and replace it
with a tab. Do not use spaces.

2) Towards the end of the "global :" target, uncomment the line,
"#cd MT/src             ; make -f makeGlobalLib". To do this, delete "#" and
replace it with a tab. Do not use spaces. 

Next, some edits must be made to the "Make.inc" file. This is to ensure that
the library is compiled with the correct compilers and flags:

1) Around Line 14/15, set: 
CC=gcc
Make sure to comment out the default "lang-4.0" compiler line.

2) Around Line 38, some additional C-flags need to be appended. I was
originally getting some linking errors later when trying to link Spooles
to Calculix, so I used the following line here:
CFLAGS = $(OPTLEVEL) -fopenmp -pthread -fPIE

3) Around Line 54, comment out the current THREAD_LIBS line and 
uncomment the last one so that it just reads as:
THREAD_LIBS = -lpthread

The Spooles library is now ready to be compiled. Use the following
command to compile a static library of Spooles,

make -f makefile lib

If all goes well, a static library called "spooles.a" should be
generated. 


---------- COMPILE CALCULIX ----------

Download the latest version of Calculix source code from:
http://www.dhondt.de/

Specifically, the download link for version 2.20 is:
http://www.dhondt.de/ccx_2.20.src.tar.bz2

Go to the "src" subdirectory and note that a Makefile has already been provided.
However, this default Makefile required some edits for me to make things work. 
So, a custom make file is provided, which was based on the provided 
Makefile_MT provided by feacluster.com. The custom Makefile is copied below
for completeness, which I have called "Makefile_GNU_MT":


# Update these paths to point to the directories that contain spooles.a and libarpack_GNU.a
# Do not leave a trialing slash in these path variables
SPOOLES_DIR = ../../spooles
ARPACK_DIR = ../../arpack-ng-master/build

CFLAGS = -Wall -O2 -fopenmp -I $(SPOOLES_DIR) -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DUSE_MT=1 -pthread -fPIE
FFLAGS = -Wall -O2 -fopenmp -pthread -fPIE

CC=gcc
FC=gfortran

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = ccx_2.20.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

LIBS = \
       $(SPOOLES_DIR)/MT/src/spoolesMT.a \
       $(SPOOLES_DIR)/spooles.a \
       $(ARPACK_DIR)/libarpack_GNU.a \
       -lpthread -lm -lopenblas

ccx_2.20_MT: $(OCCXMAIN) ccx_2.20_MT.a  $(LIBS)
	perl date.pl; $(CC) $(CFLAGS) -c ccx_2.20.c; $(FC) -fopenmp -Wall -O2 -no-pie -o $@ $(OCCXMAIN) ccx_2.20_MT.a $(LIBS)

ccx_2.20_MT.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?


Make sure to change the SPOOLES_DIR and ARPACK_DIR paths appropriately so that
they correctly point to the static Spooles and Arpack-NG libraries. I have
named this custom makefile, "Makefile_GNU_MT". Make sure it is saved in the
"src" directory, and then call the makefile,

make -f Makefile_GNU_MT

Upon successful completion, an executable will be created called "ccx_2.20_MT".
That should be it! Place this executable in whatever directory that you 
desire when it needs to be called to run a simulation. Also, this executable
is dynamically linked to the OpenBLAS libaries, so this executable may not
run if it is copied to a new Linux system which does not have these installed.
However, the machine that was used to successfully compile CalculiX will 
always have the OpenBLAS libraries in the appropriate PATH variables.


---------- TROUBLESHOOTING ----------

I ran into a few issues the first time I tried to compile
Calculix Version 2.20. First, the Perl script "date.pl" seemed to be trying
to edit files that did not exist. So, I changed it to the following:


#!/usr/bin/env perl

chomp($date=`date`);

# inserting the date into ccx_2.20.c

@ARGV="ccx_2.20.c";
$^I=".old";
while(<>){
    s/You are using an executable made on.*/You are using an executable made on $date\\n");/g;
    print;
}

# inserting the date into CalculiXstep.c

@ARGV="CalculiXstep.c";
$^I=".old";
while(<>){
    s/You are using an executable made on.*/You are using an executable made on $date\\n");/g;
    print;
}

# inserting the date into frd.c

@ARGV="frd.c";
$^I=".old";
while(<>){
    s/COMPILETIME.*/COMPILETIME       $date                    \\n\",p1);/g;
    print;
}

system "rm -f ccx_2.20.c.old";
system "rm -f CalculiXstep.c.old";
system "rm -f frd.c.old";


I also ran into an issue where a particular file that was called from the
"Makefile.inc" could not be found. I forgot which one it is, but it would seem
that CalculiX Version 2.20 shipped with one of the required files missing. I 
was able to download CalculiX Version 2.19 to get a copy of the missing file:

http://www.dhondt.de/ccx_2.19.src.tar.bz2


---------- ADD UMAT TO CCX INSTALLATION ----------

Once base CalculiX is installed, copy "umat.f" into:

<path-to-ccx>/ccx_2.20/src/

and re-compile (i.e., rerun "make -f Makefile_GNU_MT"). We tend to put the
resulting binary in /bin/ and call it something like "ccx_2.20_MT_SLS_umat".
