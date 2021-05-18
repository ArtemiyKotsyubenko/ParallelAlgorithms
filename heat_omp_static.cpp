#include <omp.h>
#include <iostream>
#include <vector>

constexpr int number_of_partitions = 10000;
constexpr long long time_ = 0.004 * number_of_partitions * number_of_partitions / 0.3;


int main(int argc, char** argv) {
    double time_start, time_end;
    time_start = omp_get_wtime();



    std::vector<double> curr_arr(number_of_partitions, 0);
    std::vector<double> prev_arr(number_of_partitions, 0);
    prev_arr[0] =  curr_arr[0] = 1;

#pragma omp parallel
    {
        for (long long t = 0; t < time_; ++t) {
//#pragma omp barrier
#pragma omp for schedule (static)
            for (long long i = 1; i < number_of_partitions - 1; i++) {
                curr_arr[i] = prev_arr[i] + 0.3 * (prev_arr[i - 1] - 2 * prev_arr[i] + prev_arr[i + 1]);
            }
//#pragma omp barrier
#pragma omp single
            {
                curr_arr.swap(prev_arr);
            };
        }
    }


    std::cout << "time " << omp_get_wtime() - time_start << std::endl;

    prev_arr.resize(prev_arr.size() - 1);
    for (auto item : prev_arr){
        std::cout << item << std::endl;
    }

    return 0;
}
