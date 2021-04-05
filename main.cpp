#include <iostream>
#include "time.h"
#include <mpi.h>
#include <cmath>
#include <vector>
#include <iomanip>


class Timer {
public:
    Timer() {
        begin = clock();
    }

    ~Timer() {
        auto end = clock();
        auto elapsed_ms = (end - begin) * 1000 / CLOCKS_PER_SEC;
        std::cout << elapsed_ms << " ms" << std::endl;
    }

private:
    clock_t begin;
};

class MPI_RAI {
public:

    const int root = 0;
    const int tag = 0;
    int rank = 0;
    int size = 0;

    MPI_RAI(int argc, char *argv[]) {
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
    }

    ~MPI_RAI() {
        MPI_Finalize();
    }
};

class Calculation {
public:
    Calculation(MPI_RAI &MPI_, const int number_of_partitions) :
            MPI(MPI_),
            number_of_partitions_(number_of_partitions) {}

    std::vector<double> run();


private:
    MPI_RAI &MPI;
    const int end_ = MPI.size - 1;
    const int number_of_partitions_;

    std::pair<int, int> synchronize(const int rank, const int left, const int right);

    std::vector<double> calc(const int rank, const int cells_cnt);


};

std::pair<int, int> Calculation::synchronize(const int rank, const int own_left, const int own_right) {

    MPI_Status status;
    double other_left{0}, other_right{0};

    if (MPI.rank % 2 == 0) {
        if (MPI.rank != MPI.root) {
            MPI_Recv(&other_right, 1, MPI_DOUBLE, MPI.rank - 1, MPI.tag, MPI_COMM_WORLD, &status);
        }

        if (MPI.rank != end_) {
            MPI_Send(&own_right, 1, MPI_DOUBLE, MPI.rank + 1, MPI.tag, MPI_COMM_WORLD);
        }


        if (MPI.rank != MPI.root) {
            MPI_Send(&own_left, 1, MPI_DOUBLE, MPI.rank - 1, MPI.tag, MPI_COMM_WORLD);
        }

        if (MPI.rank != end_) {
            MPI_Recv(&other_left, 1, MPI_DOUBLE, MPI.rank + 1, MPI.tag, MPI_COMM_WORLD, &status);
        }

    } else if (MPI.rank % 2 == 1) {
        MPI_Recv(&other_right, 1, MPI_DOUBLE, MPI.rank - 1, MPI.tag, MPI_COMM_WORLD, &status);

        if (MPI.rank != end_) {
            MPI_Send(&own_right, 1, MPI_DOUBLE, MPI.rank + 1, MPI.tag, MPI_COMM_WORLD);
        }

        MPI_Send(&own_left, 1, MPI_DOUBLE, MPI.rank - 1, MPI.tag, MPI_COMM_WORLD);

        if (MPI.rank != end_) {
            MPI_Recv(&other_left, 1, MPI_DOUBLE, MPI.rank + 1, MPI.tag, MPI_COMM_WORLD, &status);
        }
    }

    return std::make_pair(other_left, other_right);
}

std::vector<double> Calculation::calc(const int rank, const int cells_cnt) {

}

std::vector<double> Calculation::run() {

    MPI_Status status;
    std::vector<double> result;
    result.reserve(number_of_partitions_);
    Timer *timer = new Timer;

    // распределить

    std::vector<int> cells_per_index(MPI.size, number_of_partitions_ / MPI.size);
    int division_reminder = number_of_partitions_ % MPI.size;
    for (int i = 0; division_reminder > 0; ++i, --division_reminder) {
        ++cells_per_index[i];
    }


    if (MPI.rank == MPI.root) {
        auto vec = calc(MPI.root, cells_per_index[MPI.root]);
        for (int i = 1; i < MPI.size; ++i) {
            std::vector<int> buff(cells_per_index[i]);
            MPI_Recv(&buff.front(), cells_per_index[i], MPI_DOUBLE, i, MPI.tag, MPI_COMM_WORLD, &status);
            std::copy(buff.begin(), buff.end(), std::back_inserter(result));
        }

    } else {
        auto vec = calc(MPI.rank, cells_per_index[MPI.rank]);
        MPI_Send(&vec.front(), vec.size(), MPI_DOUBLE, MPI.rank, MPI.tag, MPI_COMM_WORLD);
    }


    return result;
}


int main(int argc, char **argv) {

    MPI_RAI MPI(argc, argv);
    int number_of_partitions = 4;
    Calculation calc(MPI, 4);
    std::vector<double> result = calc.run();

    for(int i = 0; i < number_of_partitions; ++i){
        std::cout  << "U(" << std::setw( 10 )<< std::left << i <<", T) = "   << result[i] << std::endl;
    }


/*    if (MPI.rank % 2 == 0) {

        if (MPI.rank == MPI.root) {

        } else if (MPI.rank == MPI.size - 1) {

        } else {

        }

    } else if (MPI.rank % 2 == 1) {

        if (MPI.rank == MPI.size - 1) {

        } else {

        }
    }*/
}

#define PORTO_START 0.0 //integrate from, double
#define PORTO_FINISH 1.0 //integrate to, double
#define PORTO_COUNT 10000l //count of dots
#define PORTO_LEN (PORTO_FINISH - PORTO_START)
#define PORTO_T 0.001
#define PORTO_h (PORTO_LEN/PORTO_COUNT)
#define PORTO_t  0.3*PORTO_h*PORTO_h
#define PORTO_Graph_dot_count 100
#define PORTO_fileskip (PORTO_COUNT/PORTO_Graph_dot_count)
#define PORTO_FILEOUT true

class MPI_unit {
    const int Tag = 0;
    const int root = 0;
    int rank = 0, commSize = 0;

    void fun(double *arr, long long count);

    double left = 0.3, right = 1;
    double h = PORTO_h, t = PORTO_t, T = PORTO_T;
    long long dotcount = PORTO_COUNT, fileskip = PORTO_fileskip, T_n = T / t;

    inline double um(double u0, double u1, double u2) {
        return u1 + t / (h * h) * (u2 - 2 * u1 + u0);
    }

public:
    //int getrank() const {return rank;}
    //int getcommSize() const {return  commSize;}
    MPI_unit(int argc, char *argv[]) {
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    }

    ~MPI_unit() { MPI_Finalize(); }

    void run();

    void exchange(double *recL, double *recR, double sendL, double sendR);
};

void MPI_unit::exchange(double *recL, double *recR, double sendL, double sendR) {
    MPI_Status status;
    if (commSize == 1) {
        *recR = right;//might be function
        *recL = left;//might be function
        return;
    }
    if (rank == 0) {
        *recL = left;//might be function
        MPI_Send(&sendR, 1, MPI_DOUBLE, rank + 1, Tag, MPI_COMM_WORLD);
        MPI_Recv(recR, 1, MPI_DOUBLE, rank + 1, Tag, MPI_COMM_WORLD, &status);
    } else if (rank == commSize - 1) {
        *recR = right;//might be function
        if (rank % 2 == 0) {
            MPI_Send(&sendL, 1, MPI_DOUBLE, rank - 1, Tag, MPI_COMM_WORLD);
            MPI_Recv(recL, 1, MPI_DOUBLE, rank - 1, Tag, MPI_COMM_WORLD, &status);
        } else {
            MPI_Recv(recL, 1, MPI_DOUBLE, rank - 1, Tag, MPI_COMM_WORLD, &status);
            MPI_Send(&sendL, 1, MPI_DOUBLE, rank - 1, Tag, MPI_COMM_WORLD);
        }

    } else if (rank % 2 == 0) {
        MPI_Send(&sendL, 1, MPI_DOUBLE, rank - 1, Tag, MPI_COMM_WORLD);
        MPI_Recv(recL, 1, MPI_DOUBLE, rank - 1, Tag, MPI_COMM_WORLD, &status);
        MPI_Send(&sendR, 1, MPI_DOUBLE, rank + 1, Tag, MPI_COMM_WORLD);
        MPI_Recv(recR, 1, MPI_DOUBLE, rank + 1, Tag, MPI_COMM_WORLD, &status);
    } else {
        MPI_Recv(recR, 1, MPI_DOUBLE, rank + 1, Tag, MPI_COMM_WORLD, &status);
        MPI_Send(&sendR, 1, MPI_DOUBLE, rank + 1, Tag, MPI_COMM_WORLD);
        MPI_Recv(recL, 1, MPI_DOUBLE, rank - 1, Tag, MPI_COMM_WORLD, &status);
        MPI_Send(&sendL, 1, MPI_DOUBLE, rank - 1, Tag, MPI_COMM_WORLD);
    }
}


void MPI_unit::fun(double *arr, long long count) {
    for (long long j = 0; j < T_n; ++j) {
        exchange(arr - 1, arr + count, *(arr), *(arr + count - 1));
        double copy = arr[-1];
        //     val
        //      |
        //v0 - v1 - v2
        for (long long i = 0; i < count; ++i) {
            double c = arr[i];
            arr[i] = um(copy, arr[i], arr[i + 1]);
            copy = c;
        }
    }
}

std::pair<double, double> synchronize() {

}

void MPI_unit::run() {
    //fprintf(stderr, "run working, rank = %d\n", rank);

    if (rank == root) {
        Timer timer;
        //MPI root
        //fprintf(stderr, "I'm root, commSize is %d\n", commSize);//print config information
        MPI_Status status;
        long long count = PORTO_COUNT;
        double start = PORTO_START, finish = PORTO_FINISH;

        double *arr = new double[count + 2]{};//line array with whole picture to fill
        ++arr;
        long long partSize = count / commSize;//count_of_dots to each task
        long long shift = count % commSize;//remaining dots for distribution
        long long *msg = new long long[2 * commSize];//distribution massive//0 - number of start, 1 - count
        /*massive with start index on 2*i positions and
        count of lines on 2*i+1 positions. 0 is root,
        i = 1..commSize-1 are clients*/
        for (int i = root; i < shift; ++i) {
            msg[2 * i] = (partSize + 1) * i;
            msg[2 * i + 1] = partSize + 1;
        } //clients with count_of_dots partSize+1 to distribute remainig dots
        for (int i = shift; i < commSize; ++i) {
            msg[2 * i] = partSize * i + shift;
            msg[2 * i + 1] = partSize;
        } //clients with count_of_dots partSize
        for (int i = root + 1; i < commSize; ++i) {
            //fprintf(stderr, "I'm root, send to %d start %lld count %lld\n", i, msg[2*i], msg[2*i+1]);//printing configs
            MPI_Send(msg + 2 * i, 2, MPI_LONG_LONG, i, Tag, MPI_COMM_WORLD);
            //MPI_Send(arr + msg[2*i], msg[2*i+1], MPI_DOUBLE, i, Tag, MPI_COMM_WORLD);
            //sending tasks
        }
        //INTEGRATE
        fun(arr, msg[1]);
        //fprintf(stderr, "I'm %d, local res is\n", rank);

        if (PORTO_FILEOUT)
            for (int i = root + 1; i < commSize; ++i) {
                MPI_Recv(arr + msg[2 * i], msg[2 * i + 1], MPI_DOUBLE, i, Tag, MPI_COMM_WORLD, &status);

            }
        if (PORTO_FILEOUT) {
            FILE *file;
            file = fopen("data.txt", "w");
            for (long long i = 0; i < dotcount; i += fileskip)
                fprintf(file, "%lf %lf\n", i * h, arr[i]);
            fclose(file);
        }
        delete[] msg;
        --arr;
        delete[] arr;
        //printf("#####\nCommSize is %d \nDotcount is %ld\n", commSize, PORTO_COUNT);
    } else {
        //MPI client
        MPI_Status status;
        long long msg[2];//two ints to receive task
        //usleep(1000 + 100*rank);
        MPI_Recv(msg, 2, MPI_LONG_LONG, root, Tag, MPI_COMM_WORLD, &status);
        double *arr = new double[msg[1] + 2]{};
        ++arr;
        //MPI_Recv(arr, msg[1], MPI_LONG_LONG, root, Tag, MPI_COMM_WORLD, &status);
        //usleep(1000 + 100*rank);
        //fprintf(stderr, "I'm %d, start is %d, count is %d\n", rank, msg[0], msg[1]);

        long long count = PORTO_COUNT;
        double start = PORTO_START, finish = PORTO_FINISH;
        //INTEGRATE
        fun(arr, msg[1]);

        //fprintf(stderr, "I'm %d, local res is %f\n", rank, res);
        if (PORTO_FILEOUT)
            MPI_Send(arr, msg[1], MPI_DOUBLE, root, Tag, MPI_COMM_WORLD);//sending results
        --arr;
        delete[]arr;
    }

}

//compile:  mpic++ main.cpp
//run:      mpirun -np 8 ./a.out
//to know time run: time mpirun -np 8 ./a.out
