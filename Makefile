#
#	File name		: Incompressible Viscous
#	Date				: May 2019
#	Version			: beta
#	Author			: Adhika, Albert, Zakiy
#

DEST			= .

HDRS			= global.hpp \
						src/initialization.hpp \
						src/SIMPLE.hpp \

LIBS			=

INPS			=

COMPILER		= g++

OPT_STD 		= -std=c++11

OPT_O 			= -O3

OPTFLAG			= $(OPT_STD) $(OPT_O)

MAKEFILE		= Makefile


PROGRAM			= main

SRCS			= main.cpp \
						global.cpp \
						src/initialization.cpp \
						src/MomentumEq.cpp \
						src/pressure_correction.cpp \
						src/get_simple.cpp \

OBJS			= $(SRCS:.cpp=.o)

.cpp.o:
			$(COMPILER) $(OPTFLAG) -c $*.cpp -o $*.o


all:			$(PROGRAM)

$(PROGRAM):		$(OBJS) $(LIBS)
				@echo -n "Loading Program $(PROGRAM) ... "
				@$(COMPILER) $(OPTFLAG) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
				@echo "done"

clean:;			@rm -f $(SRCS:.cpp=.o) $(SRCS:.cpp=.il) $(SRCS:.cpp=.d) $(PROGRAM)
			@echo "-------------- CLEANED-----------------"
delete_data:
			@while [ -z "$$CONTINUE" ]; do \
			read -r -p "Do you really want to delete data files [y/N] ? " CONTINUE; \
			done ; \
			if [ $$CONTINUE != "y" ] && [ $$CONTINUE != "Y" ]; then \
			echo "Exiting." ; exit 1 ; \
			fi
			@rm -f output/*.dat output/*.csv
			@echo "-------------- DELETED-----------------"
