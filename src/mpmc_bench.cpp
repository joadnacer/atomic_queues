#include <cstddef>
#include <iostream>
#include <thread>
#include <chrono>
#include <iomanip>
#include <vector>
#include <functional> 

#include "atomic_queues.hpp"

constexpr uint64_t total_rounds = 100'000'000;
constexpr size_t capacity = 65536;

using jdz_queue_cmp_pow2 = typename jdz::MpmcQueue<uint64_t, capacity, jdz::EnforcePowerOfTwo>;
using jdz_queue_run_pow2 = typename jdz::MpmcQueue<uint64_t, 0, jdz::EnforcePowerOfTwo>;
using jdz_queue_cmp = typename jdz::MpmcQueue<uint64_t, capacity + 1, jdz::DoNotEnforcePowerOfTwo>;
using jdz_queue_run = typename jdz::MpmcQueue<uint64_t, 0, jdz::DoNotEnforcePowerOfTwo>;

struct QueueType {
    std::string name;
    std::function<uint64_t(int, size_t)> benchmark;
    std::size_t capacity;
};

template <typename T>
void mpmc_read_worker(T &queue, uint64_t rounds) {
    uint64_t round = 0;
    uint64_t val = 0;

    while (round++ < rounds) {
        queue.pop(val);
    }
}

template <typename T>
void mpmc_write_worker(T &queue, uint64_t rounds) {
    uint64_t round = 0;

    while (round++ < rounds) {
        queue.push(round);
    }
}

template <typename T>
auto bounded_mpmc_queue_bench(int num_threads, size_t size) {
    T queue(size);

    std::vector<std::thread> workers(num_threads * 2);

    const auto begin_time = std::chrono::system_clock::now();

    for (int i = 0; i < num_threads; i++) {
        new (&workers[i]) std::thread([&queue, num_threads]() {
            mpmc_write_worker(queue, total_rounds / num_threads);
        });
    }

    for (int i = 0; i < num_threads; i++) {
        new (&workers[num_threads + i]) std::thread([&queue, num_threads]() {
            mpmc_read_worker(queue, total_rounds / num_threads);
        });
    }

    for (int i = 0; i < num_threads * 2; i++) {
        workers[i].join();
    }

    const auto end_time = std::chrono::system_clock::now();

    return (end_time - begin_time).count();
}

void run_bench(int num_threads, const std::vector<QueueType>& queue_types) {
    fprintf(stdout, "=== Num Producers=%d - Num Consumers=%d ===\n", num_threads, num_threads);

    const int time_width = 15;
    const int label_width = 15;

    for (const auto& queue_type : queue_types) {
        auto time_us = queue_type.benchmark(num_threads, queue_type.capacity) / 1000;
        std::cout << std::setw(label_width) << std::left << queue_type.name
                  << std::setw(time_width) << std::right << time_us << "us | rounds = " << total_rounds << std::endl;
    }
}

int main(int argc, char *argv[]) {
    std::vector<QueueType> queue_types = {
        {"jdz-cmp-pow2", bounded_mpmc_queue_bench<jdz_queue_cmp_pow2>, capacity},
        {"jdz-run-pow2", bounded_mpmc_queue_bench<jdz_queue_run_pow2>, capacity},
        {"jdz-cmp", bounded_mpmc_queue_bench<jdz_queue_cmp>, capacity + 1},
        {"jdz-run", bounded_mpmc_queue_bench<jdz_queue_run>, capacity + 1}
    };

    for (int i = 1; i < argc; i++) {
        run_bench(std::atoi(argv[i]), queue_types);
    }

    return 0;
}