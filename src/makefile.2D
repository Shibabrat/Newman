ND=2

CC=mpic++
CFLAGS= -DND=$(ND) -c -W -Wall
LDFLAGS = -lgsl -lgslcblas
SOURCES=newman.cpp  parallel.cpp errors.cpp boundary.cpp tracers.cpp integrate.cpp velocity.cpp ftle.cpp parameters.cpp data.cpp eqnparser.cpp mapcoord.cpp

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=newman$(ND)D

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:		
	rm -f *.o

install:	$(EXECUTABLE)
		sudo mv $(EXECUTABLE) /usr/local/bin
