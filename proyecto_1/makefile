OBJECTS=proyecto_1.o
OUTPUT=proyecto_1

LDLIBS=-lX11 -lglut -lGLU -lGL -lm

$(OUTPUT): $(OBJECTS)
	cc $(CFLAGS) $(LDFLAGS) $(LDLIBS) -o $(OUTPUT) $(OBJECTS)

$(OBJECTS): proyecto_1.h

clean:
	rm -f *.o
	rm -f test
