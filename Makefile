CC=gcc
CFLAGS=-c -fPIC
LDFLAGS=-shared
LIBDIR=lib
OBJDIR=obj
SRCDIR=src

all: $(LIBDIR)/touchpadlib.so

$(LIBDIR)/touchpadlib.so: $(OBJDIR)/touchpadlib.o
	-@mkdir $(LIBDIR)
	$(CC) $(LDFLAGS) -Wl,-soname,touchpadlib.so -o $(LIBDIR)/touchpadlib.so $(OBJDIR)/touchpadlib.o

$(OBJDIR)/touchpadlib.o: $(SRCDIR)/touchpadlib.c $(SRCDIR)/touchpadlib.h
	-@mkdir $(OBJDIR)
	$(CC) $(CFLAGS) -o $(OBJDIR)/touchpadlib.o $(SRCDIR)/touchpadlib.c

clean:
	-@rm -f $(OBJDIR)/*.o $(LIBDIR)/touchpadlib.so 2>/dev/null || true
	-@rmdir $(LIBDIR) 2>/dev/null || true
	-@rmdir $(OBJDIR) 2>/dev/null || true
