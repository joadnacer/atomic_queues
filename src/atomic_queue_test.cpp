#define CATCH_CONFIG_MAIN

#include <atomic>
#include <vector>
#include <random>
#include <thread>
#include <cassert>
#include <catch2/catch_template_test_macros.hpp>

#include "atomic_queues.hpp"

using namespace jdz;

#define RESET_CONSTRUCTORS constructor_count = 0; destructor_count = 0;

#define REQUIRE_CONSTRUCTORS(CONS, DES) REQUIRE(constructor_count == CONS); REQUIRE(destructor_count == DES);

#define TEST_CONSTRUCTORS(METHOD, NUM) RESET_CONSTRUCTORS METHOD; REQUIRE_CONSTRUCTORS(NUM, NUM) 

#define FUZZ_ROUNDS 1'000'000

static size_t constructor_count = 0;
static size_t destructor_count = 0;

class TestElement {
public:
    TestElement() noexcept : value(0) {
        constructor_count++;
    }

    TestElement(int value)  noexcept: value(value) {
        constructor_count++;
    }

    TestElement(const TestElement& other) noexcept {
        constructor_count++;
        value = other.value;
    }

    ~TestElement() noexcept {
        destructor_count++;
    }

private:
    int value;
};

static TestElement elem;

template <typename T>
void testDoubleFillAndEmpty(size_t capacity) {
    T queue(capacity);

    for (int i = 0; i < 2; i++) {
        REQUIRE(queue.empty());
        REQUIRE(queue.size() == 0);

        for (int j = 1; j <= capacity; j++)  {
            queue.emplace();
            REQUIRE(queue.size() == j);
        }

        for (int j = capacity - 1; j >= 0; j--) {
            queue.pop(elem);
            REQUIRE(queue.size() == j);
        }
    }
}

template <typename T>
void testAllMethods(size_t capacity) {
    assert(capacity >= 2 && "test requires capacity >= 2");

    T queue(capacity);

    REQUIRE(queue.empty());
    REQUIRE(queue.size() == 0);

    queue.emplace();
    queue.push(elem);

    REQUIRE(!queue.empty());
    REQUIRE(queue.size() == 2);

    queue.pop(elem);

    REQUIRE(queue.size() == 1);
    REQUIRE(queue.try_pop(elem) == true);
    REQUIRE(queue.size() == 0);
    REQUIRE(queue.try_pop(elem) == false);

    REQUIRE(queue.try_emplace() == true);
    REQUIRE(queue.try_push(elem) == true);
    REQUIRE(queue.size() == 2);

    for (size_t i = 2; i < capacity; i++) {
        queue.emplace();
    }

    REQUIRE(queue.size() == capacity);
    REQUIRE(queue.try_emplace() == false);
    REQUIRE(queue.try_push(elem) == false);
    REQUIRE(queue.size() == capacity);
}

template <typename T>
void fuzz_write_worker(T &queue, std::atomic<uint64_t> &global_sum, size_t rounds, size_t start_val) {
    uint64_t sum = 0;

    for (size_t i = start_val; i < rounds + start_val; i++) {
        queue.emplace(i);
        sum += i;
    }

    global_sum += sum;
}

template <typename T>
void fuzz_read_worker(T &queue, std::atomic<uint64_t> &global_sum, size_t rounds) {
    uint64_t sum = 0;

    for (size_t i = 0; i < rounds; i++) {
        uint64_t val;

        queue.pop(val);

        sum += val;
    }

    global_sum += sum;
}

template <typename T>
void fuzzTest(size_t num_threads, size_t capacity) {
    std::atomic<uint64_t> write_sum{0};
    std::atomic<uint64_t> read_sum{0};

    std::vector<std::thread> threads(num_threads * 2);

    T queue(capacity);

    for (size_t i = 0; i < num_threads; i++) {
        size_t r = i + num_threads;

        uint64_t rand_start = rand();

        threads[i] = std::thread(fuzz_write_worker<T>, std::ref(queue), std::ref(write_sum), FUZZ_ROUNDS/num_threads, rand_start);
        threads[r] = std::thread(fuzz_read_worker<T>, std::ref(queue), std::ref(read_sum), FUZZ_ROUNDS/num_threads);
    }

    for (size_t i = 0; i < num_threads * 2; i++) {
        threads[i].join();
    }

    REQUIRE(write_sum == read_sum);
}

TEMPLATE_TEST_CASE("MpmcQueueTest PowerOfTwo", "[unit][mpmcqueue]", 
        (MpmcQueue<TestElement, 0, EnforcePowerOfTwo>),
        (MpmcQueue<TestElement, 4, EnforcePowerOfTwo>),
        (MpmcQueue<TestElement, 4, EnforcePowerOfTwo, UseStackBuffer>),
        (MpmcQueue<TestElement, 0, DoNotEnforcePowerOfTwo>),
        (MpmcQueue<TestElement, 4, DoNotEnforcePowerOfTwo>),
        (MpmcQueue<TestElement, 4, DoNotEnforcePowerOfTwo, UseStackBuffer>)) {
    TEST_CONSTRUCTORS(testDoubleFillAndEmpty<TestType>(4), 8);

    TEST_CONSTRUCTORS(testAllMethods<TestType>(4), 6);
}

TEMPLATE_TEST_CASE("MpmcQueueTest NonPowerOfTwo", "[unit][mpmcqueue]", 
        (MpmcQueue<TestElement, 0, DoNotEnforcePowerOfTwo>),
        (MpmcQueue<TestElement, 5, DoNotEnforcePowerOfTwo>),
        (MpmcQueue<TestElement, 5, DoNotEnforcePowerOfTwo, UseStackBuffer>)) {
    TEST_CONSTRUCTORS(testDoubleFillAndEmpty<TestType>(5), 10);

    TEST_CONSTRUCTORS(testAllMethods<TestType>(5), 7);
}

TEMPLATE_TEST_CASE("SpscQueueTest PowerOfTwo", "[unit][spscqueue]", 
        (SpscQueue<TestElement, 0, EnforcePowerOfTwo>),
        (SpscQueue<TestElement, 4, EnforcePowerOfTwo>),
        (SpscQueue<TestElement, 4, EnforcePowerOfTwo, UseStackBuffer>),
        (SpscQueue<TestElement, 0, DoNotEnforcePowerOfTwo>),
        (SpscQueue<TestElement, 4, DoNotEnforcePowerOfTwo>),
        (SpscQueue<TestElement, 4, DoNotEnforcePowerOfTwo, UseStackBuffer>)) {
    TEST_CONSTRUCTORS(testDoubleFillAndEmpty<TestType>(4), 8);

    TEST_CONSTRUCTORS(testAllMethods<TestType>(4), 6);
}

TEMPLATE_TEST_CASE("SpscQueueTest NonPowerOfTwo", "[unit][spscqueue]", 
        (SpscQueue<TestElement, 0, DoNotEnforcePowerOfTwo>),
        (SpscQueue<TestElement, 5, DoNotEnforcePowerOfTwo>),
        (SpscQueue<TestElement, 5, DoNotEnforcePowerOfTwo, UseStackBuffer>)) {
    TEST_CONSTRUCTORS(testDoubleFillAndEmpty<TestType>(5), 10);

    TEST_CONSTRUCTORS(testAllMethods<TestType>(5), 7);
}

TEMPLATE_TEST_CASE("MpmcFuzzQueueTest PowerOfTwo", "[fuzz][mpmcqueue]", 
        (MpmcQueue<uint64_t, 0, EnforcePowerOfTwo>),
        (MpmcQueue<uint64_t, 8192, EnforcePowerOfTwo>),
        (MpmcQueue<uint64_t, 8192, EnforcePowerOfTwo, UseStackBuffer>),
        (MpmcQueue<uint64_t, 0, DoNotEnforcePowerOfTwo>),
        (MpmcQueue<uint64_t, 8192, DoNotEnforcePowerOfTwo>),
        (MpmcQueue<uint64_t, 8192, DoNotEnforcePowerOfTwo, UseStackBuffer>)) {
    fuzzTest<TestType>(1, 8192);
    fuzzTest<TestType>(4, 8192);
    fuzzTest<TestType>(8, 8192);
}

TEMPLATE_TEST_CASE("MpmcFuzzQueueTest PowerOfTwo MinSize", "[fuzz][mpmcqueue]", 
        (MpmcQueue<uint64_t, 0, EnforcePowerOfTwo>),
        (MpmcQueue<uint64_t, 2, EnforcePowerOfTwo>),
        (MpmcQueue<uint64_t, 2, EnforcePowerOfTwo, UseStackBuffer>),
        (MpmcQueue<uint64_t, 0, DoNotEnforcePowerOfTwo>),
        (MpmcQueue<uint64_t, 2, DoNotEnforcePowerOfTwo>),
        (MpmcQueue<uint64_t, 2, DoNotEnforcePowerOfTwo, UseStackBuffer>)) {
    fuzzTest<TestType>(1, 2);
    fuzzTest<TestType>(4, 2);
    fuzzTest<TestType>(8, 2);
}

TEMPLATE_TEST_CASE("MpmcFuzzQueueTest NonPowerOfTwo", "[fuzz][mpmcqueue]", 
        (MpmcQueue<uint64_t, 0, DoNotEnforcePowerOfTwo>),
        (MpmcQueue<uint64_t, 8193, DoNotEnforcePowerOfTwo>),
        (MpmcQueue<uint64_t, 8193, DoNotEnforcePowerOfTwo, UseStackBuffer>)) {
    fuzzTest<TestType>(1, 8193);
    fuzzTest<TestType>(4, 8193);
    fuzzTest<TestType>(8, 8193);
}

TEMPLATE_TEST_CASE("MpmcFuzzQueueTest NonPowerOfTwo MinSize", "[fuzz][mpmcqueue]", 
        (MpmcQueue<uint64_t, 0, DoNotEnforcePowerOfTwo>),
        (MpmcQueue<uint64_t, 3, DoNotEnforcePowerOfTwo>),
        (MpmcQueue<uint64_t, 3, DoNotEnforcePowerOfTwo, UseStackBuffer>)) {
    fuzzTest<TestType>(1, 3);
    fuzzTest<TestType>(4, 3);
    fuzzTest<TestType>(8, 3);
}

TEMPLATE_TEST_CASE("SpscQueueTest PowerOfTwo", "[fuzz][spscqueue]", 
        (SpscQueue<uint64_t, 0, EnforcePowerOfTwo>),
        (SpscQueue<uint64_t, 8192, EnforcePowerOfTwo>),
        (SpscQueue<uint64_t, 8192, EnforcePowerOfTwo, UseStackBuffer>),
        (SpscQueue<uint64_t, 0, DoNotEnforcePowerOfTwo>),
        (SpscQueue<uint64_t, 8192, DoNotEnforcePowerOfTwo>),
        (SpscQueue<uint64_t, 8192, DoNotEnforcePowerOfTwo, UseStackBuffer>)) {
    fuzzTest<TestType>(1, 8192);
}

TEMPLATE_TEST_CASE("SpscQueueTest PowerOfTwo MinSize", "[fuzz][spscqueue]", 
        (SpscQueue<uint64_t, 0, EnforcePowerOfTwo>),
        (SpscQueue<uint64_t, 2, EnforcePowerOfTwo>),
        (SpscQueue<uint64_t, 2, EnforcePowerOfTwo, UseStackBuffer>),
        (SpscQueue<uint64_t, 0, DoNotEnforcePowerOfTwo>),
        (SpscQueue<uint64_t, 2, DoNotEnforcePowerOfTwo>),
        (SpscQueue<uint64_t, 2, DoNotEnforcePowerOfTwo, UseStackBuffer>)) {
    fuzzTest<TestType>(1, 2);
}

TEMPLATE_TEST_CASE("SpscQueueTest NonPowerOfTwo", "[fuzz][spscqueue]", 
        (SpscQueue<uint64_t, 0, DoNotEnforcePowerOfTwo>),
        (SpscQueue<uint64_t, 8193, DoNotEnforcePowerOfTwo>),
        (SpscQueue<uint64_t, 8193, DoNotEnforcePowerOfTwo, UseStackBuffer>)) {
    fuzzTest<TestType>(1, 8193);
}

TEMPLATE_TEST_CASE("SpscQueueTest NonPowerOfTwo MinSize", "[fuzz][spscqueue]", 
        (SpscQueue<uint64_t, 0, DoNotEnforcePowerOfTwo>),
        (SpscQueue<uint64_t, 3, DoNotEnforcePowerOfTwo>),
        (SpscQueue<uint64_t, 3, DoNotEnforcePowerOfTwo, UseStackBuffer>)) {
    fuzzTest<TestType>(1, 3);
}