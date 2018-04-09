CXX		:= g++
CXXFLAGS := -MMD -MP -std=gnu++14
BINDIR := bin
SRCDIR := src
OBJDIR := obj
INCLUDE := -I./include

LIBS := 
SRCS := $(wildcard $(SRCDIR)/*.cpp)
OBJS := $(addprefix $(OBJDIR)/, $(notdir $(SRCS:.cpp=.o)))
DEPS := $(OBJS:.o=.d)

ifeq ($(OS),Windows_NT)
INCLUDE	+= -IC:/Libraries/eigen-3.3.4 -IC:/Libraries/fftw-3.3.5-dll64
LINK := -LC:/Libraries/fftw-3.3.5-dll64 -lfftw3-3 -lfftw3f-3 -lfftw3l-3
EXECUTABLE	:= main.exe
RM := cmd //C del
else
INCLUDE	+= -I/usr/include/eigen3 -I/usr/include/fftw3
LINK := -lfftw3 -lfftw3f -lfftw3l -lGL -lGLEW -lglfw3 -lX11 -lXxf86vm -lXrandr -lpthread -lXi -ldl -lXinerama -lXcursor
EXECUTABLE	:= main
RM := rm -f
endif

all: $(BINDIR)/$(EXECUTABLE)

$(BINDIR)/$(EXECUTABLE): $(OBJS) $(LIBS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LINK)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $< $(LINK)

clean:
	$(RM) $(BINDIR)/$(EXECUTABLE) $(OBJS) $(DEPS)

-include $(DEPS)
