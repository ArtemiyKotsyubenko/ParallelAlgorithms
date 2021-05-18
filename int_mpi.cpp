#include <iostream>
#include <chrono>
#include <mpi.h>
#include <cmath>


class Timer {
public:
    Timer() {
        begin = std::chrono::steady_clock::now();
    }

    ~Timer() {
        auto end = std::chrono::steady_clock::now();
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
        std::cout << elapsed_ms.count() << " ms" << std::endl;
    }

private:
    std::chrono::time_point<std::chrono::_V2::steady_clock, std::chrono::duration<long int, std::ratio<1, 1000000000>>>
            begin;


};

class MPI_RAI {
public:

    const int main = 0;
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


constexpr long fineness = 21000000000;

double func(double x) {
    return sqrt(4 - x * x);
}

double calculate_integral(int num, int threads_cnt) {
    double sum = 0;
    double size = 2.0 / fineness;
    double mid = 1.0 / fineness;

    for (double x = mid + size * num; x < 2; x += threads_cnt * size) {
        sum += func(x) * size;
    }
    return sum;
}

int main(int argc, char *argv[]) {
    MPI_RAI MPI(argc, argv);
    Timer *timer = new Timer;

    if (MPI.rank == MPI.main) {// main thread
        MPI_Status status;
        double answer = calculate_integral(MPI.rank, MPI.size);

        if (MPI.size > 1) {
            for (int i = 1; i < MPI.size; ++i) {
                double result = 0;
                MPI_Recv(&result, 1, MPI_DOUBLE, i, MPI.tag, MPI_COMM_WORLD, &status);
                answer += result;
            }
        }

        delete timer;// call destructor in single thread
        std::cout << answer << std::endl;

    } else {
        double result = calculate_integral(MPI.rank, MPI.size);
        MPI_Send(&result, 1, MPI_DOUBLE, MPI.main, MPI.tag, MPI_COMM_WORLD);
    }

}
