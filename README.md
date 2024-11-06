# atomic_queues

This repository contains a single header-file, atomic_queues.hpp, which contains bounded MpmcQueue and SpscQueue implementations for C++20. The SPSC queue appears to be the best performing SPSC queue, while the MPMC queue appears to be the best performing implementation for N threads <= N cores (and is by far the best performing at low contention).

If you require high performance at N threads > N cores, you may want to use a queue such as [moodycamel's ConcurrentQueue](https://github.com/cameron314/concurrentqueue). I am working on an original "composite" MMPMC queue implementation which will eventually be added to this repo and will hopefully be the top performing queue at high contention.

Note that performance varies greatly by system and the benchmarks in this README are not realistic, benchmark queues yourself if trying to optimize performance.

The implementation of the MPMC queue is based on [Dmitry Vyukov's bounded MPMC queue](https://www.1024cores.net/home/lock-free-algorithms/queues/bounded-mpmc-queue), and the implementation of the SPSC queue is based on [Erik Rigtor's SpscQueue](https://github.com/rigtorp/SPSCQueue). Modifications have been made to these implementations in order to improve performance and configurability.

# Usage and Configuration

Both queues' templates take the following parameters:
```c++
<typename T, size_t N = 0, ValidSizeConstraint SizeConstraint = jdz::EnforcePowerOfTwo, BufferType BufType = jdz::UseHeapBuffer>
```

- `typename T`: The type stored in the queue.
- `size_t N = 0`: the capacity of the queue - providing this as a template parameter allows the compiler to optimize the modulo operation, which may greatly improve performance at low contention (see benchmarks), and allows for compile-time checking of the parameter's validity*.
- `ValidSizeConstraint SizeConstraint = jdz::EnforcePowerOfTwo`: Constraining the capacity to power of two sizes allows for the use of a bitwise AND operation instead of a ternary operation for SPSC or a modulo for MPMC, which may improve performance. Options are `jdz::EnforcePowerOfTwo` and `jdz::DoNotEnforcePowerOfTwo`.
- `BufferType BufType = jdz::UseHeapBuffer`: Determines the buffer type - HeapBuffer or StackBuffer. Note that using StackBuffer may cause segfaults if the queue capacity is too large due to stack size limits. Pick an appropriate capacity, or adjust your stack size using ie `ulimit -s` on Linux. Options are `jdz::UseHeapBuffer` and `jdz::UseStackBuffer`.

The constructor's buffer_size argument must be provided if N is equal to 0, and follows the same restrictions as N*. Passing in a custom allocator allows for use of the queue with huge pages or in shared memory.

*Buffer capacity must be greater than 1, and must respect the power of two size constraint. If N is provided, then the buffer_size_ constructor argument must be either 0 or equal to N.

Both queues implement the following interface:

```c++
Queue(size_t buffer_size = N;
      std::allocator<jdz::details::Cell<T>> allocator = std::allocator<jdz::details::Cell<T>>());

template <typename... Args>
void emplace(Args &&... args) noexcept; // blocking 

void push(const T &data) noexcept
requires std::is_nothrow_copy_constructible_v<T>; // blocking

template <typename P>
void push(P &&data) noexcept
requires std::is_nothrow_constructible_v<T, P>; // blocking

template <typename... Args>
[[nodiscard]] bool try_emplace(Args &&... args) noexcept
requires std::is_nothrow_constructible_v<T, Args &&...>; // non-blocking

bool try_push(const T &data) noexcept
requires std::is_nothrow_copy_constructible_v<T>; // non-blocking

template <typename P>
[[nodiscard]] bool try_push(P &&data) noexcept
requires std::is_nothrow_constructible_v<T, P>; // non-blocking

void pop(T &v) noexcept; // blocking;
    
[[nodiscard]] bool try_pop(T &v) noexcept; // non-blocking

[[nodiscard]] size_t size() const noexcept;

[[nodiscard]] size_t capacity() const noexcept;

[[nodiscard]] bool empty() const noexcept;
```

The blocking methods perform better than the non-blocking methods when not-preempted. If number of threads is larger than the number of cores and preemption becomes an issue, you may achieve better performance with a retry loop over the `try*` methods, although I would recommend using a different library in this situation.

# Benchmarks
Benchmarks are included in the src/ folder of this repository, testing various versions of the queue implementations. These are currently not especially rigorous and will be improved.

Posted below are the results of benchmarking the implementations found in this repository against other notable implementations. Benchmarks are run using the blocking read/write methods where possible, ie `pop` and `push` for jdz queues. Note that these may perform differently than the non-blocking methods - ie, for vyukov style queues* (such as jdz implementations and rigtorp), these perform better as long as N threads <= N cores but perform worse than the non-blocking methods below that.

These benchmarks involve submitting small data types (uint64_t) between equal numbers of producers and consumers producing and consuming at max throughput - this is not a realistic benchmark. Contention in real systems is likely lower for same number of threads, and it is important to test your use case yourself.

I also plan to add benchmarks of the queues used as SPMC/MPSC, which may produce interesting results.

*Note that the original Vyukov implementation did not contain blocking methods.

## SPSC Benchmarks:
These benchmarks measure the ops/ms of one producer transmitting 1 billion uint64_t to one consumer. Benchmarks were run 5 times and averaged. Benchmarks are run with a capacity of 65536, or 65537 for non-power of 2 jdz trials. The queues from this repository are the best performing, followed by drogalis's queue.

Benchmarked queues are:

- `jdz-cmp-pow2`: jdz queue: `jdz::SpscQueue<uint64_t, 65536, jdz::EnforcePowerOfTwo, jdz::UseStackBuffer>` - compile-time power of 2 capacity.
- `jdz-run-pow2`: jdz queue: `jdz::SpscQueue<uint64_t, 0, jdz:::EnforcePowerOfTwo>` - runtime power of 2 capacity.
- `jdz-cmp`:      jdz queue: `jdz::SpscQueue<uint64_t, 65537, jdz::DoNotEnforcePowerOfTwo, jdz::UseStackBuffer>` - compile-time non-power of 2 capacity.
- `jdz-run`:      jdz queue: `jdz::SpscQueue<uint64_t, 0, jdz::DoNotEnforcePowerOfTwo>` - runtime non-power of 2 capacity.
- `dro`:          [Andrew Drogalis's SPSC-Queue](https://github.com/drogalis)
- `rigtorp`:      [Erik Rigtorp's SpscQueue](https://github.com/rigtorp/SPSCQueue).
- `atomic_queue`: [Maxim Egorushkin's atomic_queue](https://github.com/max0x7ba/atomic_queue) with SPSC=true.
- `deaod`:        [deaod's spsc_queue](https://github.com/Deaod/spsc_queue).
- `cml-rwcb`:     [moodycamel's BlockingReaderWriterCircularBuffer](https://github.com/cameron314/readerwriterqueue).
- `cml-rwq`:      [moodycamel's ReaderWriterQueue](https://github.com/cameron314/readerwriterqueue).

### x86_64 - Intel i7-11800H

![spscl](https://i.imgur.com/vQdPhrc.png)

## MPMC Benchmarks
These benchmarks show the throughput measured for one producer transmitting 100 million uint64_t to one consumer. Benchmarks were run 5 times and averaged. Benchmarks are run with a capacity of 8192, or 8193 for non-power of 2 jdz trials.

We can see clearly that moodycamel's queue is the best by far N threads > N cores, but performs less well below this.

Benchmarked queues are:

- `jdz-cmp-pow2`:  jdz queue: `jdz::MpmcQueue<uint64_t, 65536, jdz::EnforcePowerOfTwo, jdz::UseStackBuffer>` - compile-time power of 2 capacity.
- `jdz-run-pow2`:  jdz queue: `jdz::MpmcQueue<uint64_t, 0, jdz:::EnforcePowerOfTwo>` - runtime power of 2 capacity.
- `jdz-cmp`:       jdz queue: `jdz::MpmcQueue<uint64_t, 65537, jdz::DoNotEnforcePowerOfTwo, jdz::UseStackBuffer>` - compile-time non-power of 2 capacity.
- `jdz-run`:       jdz queue: `jdz::MpmcQueue<uint64_t, 0, jdz::DoNotEnforcePowerOfTwo>` - runtime non-power of 2 capacity.
- `rigtorp`:       [Erik Rigtorp's MpmcQueue](https://github.com/rigtorp/MPMCQueue).
- `atomic_queue`:  [Maxim Egorushkin's atomic_queue](https://github.com/max0x7ba/atomic_queue).
- `moodycamel`:    [moodycamel's ConcurrentQueue](https://github.com/cameron314/concurrentqueue).
- `es-mpmc`:       [Erez Strauss's lockfree_mpmc_queue](https://github.com/erez-strauss).
- `xenium-vyukov`: [Manuel Pöter's vyukov_bounded_queue](https://github.com/mpoeter/xenium/tree/master).

### x86_64 - Intel i7-11800H
![1p1cl](https://i.imgur.com/2aVkRSG.png)
![2p2cl](https://i.imgur.com/2jvkYWb.png)
![4p4cl](https://i.imgur.com/hjwKwZA.png)
![8p8cl](https://i.imgur.com/0ij0eo8.png)
![16p16cl](https://i.imgur.com/1ZoUIlb.png)
