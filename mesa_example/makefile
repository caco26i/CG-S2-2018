OBJECTS=mesa_example.o
OUTPUT=mesa_example

CFLAGS=-I/usr/local/Mesa-3.4/include
#LDLIBS=-lX11 -lglut -lMesaGLU -lMesaGL -lm -lXext -lXmu
LDLIBS=-lX11 -lglut -lGLU -lGL -lm -lXext -lXmu
LDFLAGS=-L/usr/local/Mesa-3.4/lib -L/usr/X11R6/lib

$(OUTPUT): $(OBJECTS)
	cc $(CFLAGS) $(LDFLAGS) $(LDLIBS) -o $(OUTPUT) $(OBJECTS)

$(OBJECTS): mesa_example.h

clean:
	rm -f *.o
	rm -f mesa_example
