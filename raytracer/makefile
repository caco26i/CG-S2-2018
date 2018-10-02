OBJECTS=raytracer.o
OUTPUT=raytracer

LDLIBS=-lX11 -lglut -lGLU -lGL -lm -lIL

$(OUTPUT): $(OBJECTS)
	cc -o $(OUTPUT) $(OBJECTS) $(LDLIBS) 

$(OBJECTS): raytracer.h

clean:
	rm -f *.o
	rm -f raytracer
