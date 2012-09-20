# Makefile to compile libsla

# define compiler etc
FF=gfortran
SOURCES=*.f
OBJECTS=*.o
LFLAGS=-fPIC -O -fno-second-underscore
LDFLAGS=-shared -fno-second-underscore
#LDFLAGS=-fno-second-underscore

all: libsla.so

# Link command to link shared library
libsla.so: compilation
	${FF} ${LDFLAGS} ${OBJECTS} -o libsla.so

compilation:
	${FF} -c ${LFLAGS} ${SOURCES}

clean:
	rm -rf ${OBJECTS} libsla.* 
