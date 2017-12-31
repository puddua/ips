CC = g++
CFLAGS = -W -Wextra -Wall -std=c++11 -g
LIB = -larmadillo -lm 
TARGET = main
OBJS =  poly.o basis.o main.o
CCPLOT = gnuplot
PLOT = plot_res.gnu 


all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIB)

%.o: %.cpp %.h
	$(CC) $(CFLAGS) $< -c -o $@ $(LIB)

clean:
	rm -f $(OBJS)
	rm -f $(TARGET)

plot:
	$(CCPLOT) $(PLOT)
