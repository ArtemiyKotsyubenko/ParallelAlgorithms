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
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        std::cout << elapsed_ms.count() << " ns\\n";
    }

private:
    std::chrono::time_point<std::chrono::_V2::steady_clock, std::chrono::duration<long int, std::ratio<1, 1000000000>>>
            begin;


};

class MPI_RAI {
public:
    int rank = 0;
    int size = 0;
    const int main = 0;
    const int tag = 0;

    MPI_RAI(int argc, char **argv) {
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
    }

    ~MPI_RAI() {
        MPI_Finalize();
    }
};


double func(double x) {
    return sqrt(4 - x * x);
}


constexpr int fineness = 10000;
constexpr int thread_cnt = 1;

double calculate_integral(int num, int threads_cnt) {
    double sum = 0;
    double size = 2.0 / fineness;
    double mid = 1.0 / fineness;
    for (double x = mid + size * num; x < 2; x += threads_cnt * size) {
        //double x = mid + fineness * i; // (2 / fineness) / 2 + fineness * i
        sum += func(x) * size;
    }
    return sum;
}

int main(int argc, char *argv[]) {
    const int thread_cnt = atoi(argv[2]);
    MPI_RAI MPI(argc, argv);
    Timer timer;

    if (MPI.rank == MPI.main) {// main thread
        MPI_Status status;
        double answer = calculate_integral(MPI.rank, thread_cnt);

        if (thread_cnt > 1) {
            for (int i = 1; i < thread_cnt; ++i) {
                double result = 0;
                MPI_Recv(&result, 1, MPI_DOUBLE, i, MPI.tag, MPI_COMM_WORLD, &status);
                answer += result;
            }
        }
        std::cout << answer << '\n';

    } else {
        double result = calculate_integral(MPI.rank, thread_cnt);
        MPI_Send(&result, 1, MPI_DOUBLE, MPI.main, MPI.tag, MPI_COMM_WORLD);
    }

}