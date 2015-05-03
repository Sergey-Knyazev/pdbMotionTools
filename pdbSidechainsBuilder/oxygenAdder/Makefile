CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=oxygen_adder.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=oxygen_adder

all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(EXECUTABLE) $(OBJECTS)
