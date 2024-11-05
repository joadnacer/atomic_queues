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

using jdz_queue_cmp_pow2 = typename jdz::SpscQueue<uint64_t, capacity, jdz::EnforcePowerOfTwo, jdz::UseStackBuffer>;
using jdz_queue_run_pow2 = typename jdz::SpscQueue<uint64_t, 0, jdz::EnforcePowerOfTwo>;
using jdz_queue_cmp = typename jdz::SpscQueue<uint64_t, capacity + 1, jdz::DoNotEnforcePowerOfTwo, jdz::UseStackBuffer>;
using jdz_queue_run = typename jdz::SpscQueue<uint64_t, 0, jdz::DoNotEnforcePowerOfTwo>;

using jdz_mpmc_queue_cmp_pow2 = typename jdz::SpscQueue<uint64_t, capacity, jdz::EnforcePowerOfTwo, jdz::UseStackBuffer>;
using jdz_mpmc_queue_run_pow2 = typename jdz::SpscQueue<uint64_t, 0, jdz::EnforcePowerOfTwo>;
using jdz_mpmc_queue_cmp = typename jdz::SpscQueue<uint64_t, capacity + 1, jdz::DoNotEnforcePowerOfTwo, jdz::UseStackBuffer>;
using jdz_mpmc_queue_run = typename jdz::SpscQueue<uint64_t, 0, jdz::DoNotEnforcePowerOfTwo>;

struct QueueType {
    std::string name;
    std::function<uint64_t(size_t)> benchmark;
    std::size_t capacity;
};

template <typename T>
void spsc_read_worker(T &queue, uint64_t rounds) {
    uint64_t round = 0;
    uint64_t val = 0;

    while (round++ < rounds) {
        queue.pop(val);
    }
}

template <typename T>
void spsc_write_worker(T &queue, uint64_t rounds) {
    uint64_t round = 0;

    while (round++ < rounds) {
        queue.push(round);
    }
}

template <typename T>
auto bounded_spsc_queue_bench(size_t size) {
    T queue(size);

    std::vector<std::thread> workers(2);

    const auto begin_time = std::chrono::system_clock::now();

    new (&workers[0]) std::thread(spsc_write_worker<T>, std::ref(queue), total_rounds);
    new (&workers[1]) std::thread(spsc_read_worker<T>, std::ref(queue), total_rounds);

    workers[0].join();
    workers[1].join();

    const auto end_time = std::chrono::system_clock::now();

    return (end_time - begin_time).count();
}

int main() {
    fprintf(stdout, "=== Num Producers=1 - Num Consumers=1 ===\n");

    const int time_width = 15;
    const int label_width = 20;

    std::vector<QueueType> queue_types = {
        {"jdz-cmp-pow2", bounded_spsc_queue_bench<jdz_queue_cmp_pow2>, capacity},
        {"jdz-run-pow2", bounded_spsc_queue_bench<jdz_queue_run_pow2>, capacity},
        {"jdz-cmp", bounded_spsc_queue_bench<jdz_queue_cmp>, capacity + 1},
        {"jdz-run", bounded_spsc_queue_bench<jdz_queue_run>, capacity + 1},
        {"jdz-mpmc-cmp-pow2", bounded_spsc_queue_bench<jdz_mpmc_queue_cmp_pow2>, capacity},
        {"jdz-mpmc-run-pow2", bounded_spsc_queue_bench<jdz_mpmc_queue_run_pow2>, capacity},
        {"jdz-mpmc-cmp", bounded_spsc_queue_bench<jdz_mpmc_queue_cmp>, capacity + 1},
        {"jdz-mpmc-run", bounded_spsc_queue_bench<jdz_mpmc_queue_run>, capacity + 1}
    };

    for (const auto& queue_type : queue_types) {
        auto time_us = queue_type.benchmark(queue_type.capacity) / 1000;
        std::cout << std::setw(label_width) << std::left << queue_type.name
                  << std::setw(time_width) << std::right << time_us << "us | rounds = " << total_rounds << std::endl;
    }
}

