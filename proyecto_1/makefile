OBJECTS=test.c
OUTPUT=test

LDLIBS=-lX11 -lglut -lGLU -lGL -lm

$(OUTPUT): $(OBJECTS)
	cc -o $(OUTPUT) $(OBJECTS) $(LDLIBS)

$(OBJECTS): test.c

clean:
	rm -f *.o
	rm -f test
