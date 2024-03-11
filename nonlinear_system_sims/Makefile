COMPILER := g++

CCFLAGS := -std=c++11 -pthread -O3

ALL_LDFLAGS := -v
ALL_LDFLAGS += $(CCFLAGS)

# operating system
HOST_OS   := $(shell uname -s 2>/dev/null | tr "[:upper:]" "[:lower:]")
TARGET_OS ?= $(HOST_OS)
ifeq (,$(filter $(TARGET_OS),linux darwin qnx android))
    $(error ERROR - unsupported value $(TARGET_OS) for TARGET_OS!)
endif

# Common includes and paths for CUDA
INCLUDES  := -I/usr/local/Cellar/glew/2.1.0/include
INCLUDES  += -I/usr/local/Cellar/glfw/3.3/include
LIBRARIES := /usr/local/Cellar/glew/2.1.0/lib/libGLEW.2.1.0.dylib
LIBRARIES += /usr/local/Cellar/glfw/3.3/lib/libglfw.3.3.dylib

# Makefile include to help find GL Libraries

# OpenGL specific libraries
ifeq ($(TARGET_OS),darwin)
 # Mac OSX specific libraries and paths to include
 LIBRARIES += -L/System/Library/Frameworks/OpenGL.framework/Libraries
 LIBRARIES += -lGL -lGLU
 ALL_LDFLAGS += -Xlinker -framework -Xlinker GLUT
else
 LIBRARIES += $(GLLINK)
 LIBRARIES += -lGL -lGLU -lglut
endif

ifeq ($(SAMPLE_ENABLED),0)
EXEC ?= @echo "[@]"
endif

################################################################################

# Target rules
all: build

build: test

camera.o:camera.cpp
	$(EXEC) $(COMPILER) $(INCLUDES) $(CCFLAGS) -o $@ -c $<

shader.o:shader.cpp
	$(EXEC) $(COMPILER) $(INCLUDES) $(CCFLAGS) -o $@ -c $<

main.o:main.cpp
	$(EXEC) $(COMPILER) $(INCLUDES) $(CCFLAGS) -o $@ -c $<

test: main.o shader.o camera.o
	$(EXEC) $(COMPILER) $(INCLUDES) $(ALL_LDFLAGS) -o $@ $+ $(LIBRARIES)
	rm -r *.o

run: build
	$(EXEC) ./test

clean:
	rm -r test

clobber: clean
