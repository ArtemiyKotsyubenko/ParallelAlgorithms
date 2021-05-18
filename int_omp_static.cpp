#include <iostream>
#include <chrono>
#include <omp.h>
#include <cmath>
#include <vector>


constexpr long fineness = 100000;

double func(double x) {
    return sqrt(4 - x * x);
}


int main(int argc, char* argv[]) {
    std::vector<double> res(omp_get_max_threads(), 0);
    double size = 2.0 / fineness;
    double mid = 1.0 / fineness;

    double time_start, time_end;
    time_start = omp_get_wtime();
#pragma omp parallel
    {
#pragma omp for schedule (static)
        for (long i = 0; i < fineness; ++i) {
            res[omp_get_thread_num()] += func(mid + i * size) * size;
        }
    }

    double answer = 0;
    for (auto item : res) {
        answer += item;
    }
    std::cout << "time " << omp_get_wtime() - time_start << std::endl;

    std::cout << answer << std::endl;


}
