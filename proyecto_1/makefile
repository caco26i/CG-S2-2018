OBJECTS=proyecto_1.o
OUTPUT=proyecto_1

LDLIBS=-lX11 -lglut -lGLU -lGL -lm -lIL

$(OUTPUT): $(OBJECTS)
	cc -o $(OUTPUT) $(OBJECTS) $(LDLIBS) 

$(OBJECTS): proyecto_1.h

clean:
	rm -f *.o
	rm -f proyecto_1
