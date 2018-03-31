CXX		:= g++
CXXFLAGS := -g -MMD -MP -Wall -Wextra -Winit-self -Wno-missing-field-initializers -std=c++17
BINDIR := ./bin
SRCDIR := ./src
OBJDIR := ./obj
INCLUDE := -I./include

LIBS :=
SRCS := $(wildcard $(SRCDIR)/*.cpp)
OBJS := $(addprefix $(OBJDIR)/, $(notdir $(SRCS:.cpp=.o)))
DEPS := $(OBJS:.o=.d)
RM := rm -f

ifeq ($(OS),Windows_NT)
INCLUDE	+= -IC:/Libraries/eigen-3.3.4
EXECUTABLE	:= main.exe
else
INCLUDE	+= -I/usr/include/eigen3
EXECUTABLE	:= main
endif

all: $(BIN)/$(EXECUTABLE)

$(BIN)/$(EXECUTABLE): $(OBJS) $(LIBS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

clean:
	$(RM) $(BIN)/$(EXECUTABLE) $(OBJS) $(DEPS)

-include $(DEPS)
