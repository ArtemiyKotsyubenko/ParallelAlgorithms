#include <iostream>
#include "time.h"
#include <mpi.h>
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

    MPI_RAI(int argc, char *argv[]) {// not effective
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
    Calculation(MPI_RAI &MPI_, const int number_of_partitions, const long long time) :
            MPI(MPI_),
            number_of_partitions_(number_of_partitions),
            time_(time) {}

    void run();


private:
    MPI_RAI &MPI;
    const int end_ = MPI.size - 1;
    const int number_of_partitions_;
    const long long time_;

    using sync_res = struct {
        double other_left;
        double other_right;
    };

    sync_res synchronize(const int rank, const double own_left, const double own_right);

    std::vector<double> calc(const int rank, const int cells_cnt);


};

Calculation::sync_res Calculation::synchronize(const int rank, const double own_left, const double own_right) {

    MPI_Status status;
    double other_left{0}, other_right{0};
    if (MPI.size == 1) {
        return {0, 0};
    }
    if (MPI.rank % 2 == 0) {

        if (MPI.rank == MPI.root) {
            MPI_Ssend(&own_right, 1, MPI_DOUBLE, MPI.rank + 1, MPI.tag, MPI_COMM_WORLD);
            MPI_Recv(&other_left, 1, MPI_DOUBLE, MPI.rank + 1, MPI.tag, MPI_COMM_WORLD, &status);

        } else if (MPI.rank == MPI.size - 1) {
            MPI_Ssend(&own_left, 1, MPI_DOUBLE, MPI.rank - 1, MPI.tag, MPI_COMM_WORLD);
            MPI_Recv(&other_right, 1, MPI_DOUBLE, MPI.rank - 1, MPI.tag, MPI_COMM_WORLD, &status);
        } else {
            MPI_Ssend(&own_left, 1, MPI_DOUBLE, MPI.rank + 1, MPI.tag, MPI_COMM_WORLD);
            MPI_Recv(&other_right, 1, MPI_DOUBLE, MPI.rank + 1, MPI.tag, MPI_COMM_WORLD, &status);

            MPI_Ssend(&own_right, 1, MPI_DOUBLE, MPI.rank - 1, MPI.tag, MPI_COMM_WORLD);
            MPI_Recv(&other_left, 1, MPI_DOUBLE, MPI.rank - 1, MPI.tag, MPI_COMM_WORLD, &status);
        }

    } else if (MPI.rank % 2 == 1) {
        if (MPI.rank == MPI.size - 1) {
            MPI_Recv(&other_right, 1, MPI_DOUBLE, MPI.rank - 1, MPI.tag, MPI_COMM_WORLD, &status);
            MPI_Ssend(&own_left, 1, MPI_DOUBLE, MPI.rank - 1, MPI.tag, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&other_left, 1, MPI_DOUBLE, MPI.rank - 1, MPI.tag, MPI_COMM_WORLD, &status);
            MPI_Ssend(&own_right, 1, MPI_DOUBLE, MPI.rank - 1, MPI.tag, MPI_COMM_WORLD);

            MPI_Recv(&other_right, 1, MPI_DOUBLE, MPI.rank + 1, MPI.tag, MPI_COMM_WORLD, &status);
            MPI_Ssend(&own_left, 1, MPI_DOUBLE, MPI.rank + 1, MPI.tag, MPI_COMM_WORLD);
        }
    }

    return {other_left, other_right};
}

inline std::vector<double> Calculation::calc(const int rank, const int cells_cnt) {
    std::vector<double> previous(cells_cnt, 0);
    std::vector<double>current(cells_cnt, 0);
    int cnt = 0;

    if (rank == MPI.root) {
        previous[0] = 1;
    }

    for (long long t = 0; t < time_; ++t) {

        sync_res s = synchronize(rank, previous.front(), previous.back());

        current[0] = (rank != MPI.root) ? previous[0] + 0.3 * (s.other_right - 2 * previous[0] + previous[1]) : 1;

        current[cells_cnt - 1] =
                previous[cells_cnt - 1] +
                0.3 * (previous[cells_cnt - 2] - 2 * previous[cells_cnt - 1] + s.other_left);



        for (int i = 1; i < cells_cnt - 1; ++i) {
            current[i] = previous[i] + 0.3 * (previous[i - 1] - 2 * previous[i] + previous[i + 1]);
        }

        current.swap(previous);
    }
    return previous;
}

void Calculation::run() {

    MPI_Status status;
    std::vector<double> result;
    result.reserve(number_of_partitions_);

    std::vector<int> cells_per_index(MPI.size, number_of_partitions_ / MPI.size);
    int division_reminder = number_of_partitions_ % MPI.size;
    for (int i = 0; division_reminder > 0; ++i, --division_reminder) {
        ++cells_per_index[i];
    }

    if (MPI.rank == MPI.root) {
        Timer* timer = new Timer();
        auto vec = calc(MPI.root, cells_per_index[MPI.root]);
        std::copy(vec.begin(), vec.end(), std::back_inserter(result));

        for (int i = 1; i < MPI.size; ++i) {
            std::vector<double> buff(cells_per_index[i]);

            MPI_Recv(&buff.front(), cells_per_index[i], MPI_DOUBLE, i, MPI.tag, MPI_COMM_WORLD, &status);
            std::copy(buff.begin(), buff.end(), std::back_inserter(result));
        }

        delete timer;

//        for (int i = 0; i < number_of_partitions_; ++i) {
//            std::cout << "U(" << std::setw(10) << std::left << i << ", T) = " << result[i] << std::endl;
//        }

    } else {
        auto vec = calc(MPI.rank, cells_per_index[MPI.rank]);
        MPI_Ssend(&vec.front(), vec.size(), MPI_DOUBLE, MPI.root, MPI.tag, MPI_COMM_WORLD);
    }


}


int main(int argc, char **argv) {

    MPI_RAI MPI(argc, argv);
    int number_of_partitions = atoi(argv[1]);
    long long time = std::strtod(argv[2], nullptr) * number_of_partitions * number_of_partitions / 0.3;
    Calculation calc(MPI, number_of_partitions, time);
    calc.run();
}// ???????????? : ?????????????????? ??????????????, ?????????????? ????????????
// ????????????????..

