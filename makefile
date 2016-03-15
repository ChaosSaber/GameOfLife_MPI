# the compiler: gcc for C program, define as g++ for C++
CC = mpicc

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  =  -Wall -std=c99 -lc -D _BSD_SOURCE 


# the build target executable:
TARGET2 = helloworld

TARGET = gameoflife
TARGET3 = gameoflife-VTKprinting

all: $(TARGET)

$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c


run: all
	mpirun -n 3 ./$(TARGET) -a -b -c
	
clean:
	$(RM) $(TARGET)
