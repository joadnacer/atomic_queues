/*
Copyright Joad Nacer

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the “Software”), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <atomic>
#include <new>
#include <cassert>
#include <array>

namespace jdz
{

struct EnforcePowerOfTwo {};

struct DoNotEnforcePowerOfTwo {};

namespace details
{

#if defined(__cpp_lib_hardware_interference_size)
static constexpr size_t cache_line = std::hardware_destructive_interference_size;
#else
static constexpr size_t cache_line = 64;
#endif

static constexpr bool is_power_of_two(size_t n) {
    return (n & (n-1)) == 0;
}

template <typename T>
concept IsEnforcePowerOfTwo = std::is_same_v<T, EnforcePowerOfTwo>;

template <typename T>
concept IsDoNotEnforcePowerOfTwo = std::is_same_v<T, DoNotEnforcePowerOfTwo>;

template <typename T>
concept IsValidSizeConstraint = IsEnforcePowerOfTwo<T> || IsDoNotEnforcePowerOfTwo<T>;

template <typename T>
concept AlignedToCacheLine = alignof(T) % cache_line == 0;

template <typename T>
concept SizeMultipleOfCacheLine = sizeof(T) % cache_line == 0;

template <typename T>
concept FalseSharingSafe = AlignedToCacheLine<T> && SizeMultipleOfCacheLine<T>;

template <size_t N, typename SizeConstraint>
concept OptionalPowerOfTwo = IsValidSizeConstraint<SizeConstraint>
                         && (IsDoNotEnforcePowerOfTwo<SizeConstraint> || is_power_of_two(N));

template <size_t N>
concept ZeroOrGreaterThanOne = N == 0 || N > 1;

template <typename T, IsValidSizeConstraint SizeConstraint, size_t N = 0>
class Buffer
{
private:
    T *buffer_;

    const size_t buffer_size_;
    const size_t buffer_mask_;

    std::allocator<T> allocator_;

public:
    explicit Buffer(const size_t buffer_size, const std::allocator<T> &allocator)
        : buffer_size_(buffer_size),
          buffer_mask_(buffer_size_ - 1),
          allocator_(allocator) {
            buffer_ = allocator_.allocate(buffer_size + 1);
        }

    ~Buffer() noexcept {
        allocator_.deallocate(buffer_, buffer_size_ + 1);
    }

    T& operator[](const size_t index) noexcept {
        if constexpr (IsEnforcePowerOfTwo<SizeConstraint>) {
            return buffer_[index & buffer_mask_]; 
        }
        else if constexpr (N != 0) {
            return buffer_[index % N];
        }
        else {
            return buffer_[index % buffer_size_];
        }
    }

    const T& operator[](const size_t index) const noexcept {
        if constexpr (IsEnforcePowerOfTwo<SizeConstraint>) {
            return buffer_[index & buffer_mask_];
        }
        else if constexpr (N != 0) {
            return buffer_[index % N];
        }
        else {
            return buffer_[index % buffer_size_];
        }
    }
};

template <bool HasSeq>
struct SeqField;

template <>
struct SeqField<true> {
    alignas(cache_line) std::atomic<size_t> seq;

    SeqField(size_t i) : seq(i) {}
};

template <>
struct SeqField<false> {
    SeqField(size_t) {}
};

template <bool HasSeq>
struct IsConstructedField;

template<>
struct IsConstructedField<true> {
    bool is_constructed;

    IsConstructedField(bool is_constructed) : is_constructed(is_constructed) {}
};

template<>
struct IsConstructedField<false> {
    #ifdef __aarch64__
    alignas(cache_line) bool is_constructed;
    #else
    bool is_constructed;
    #endif

    IsConstructedField(bool is_constructed) : is_constructed(is_constructed) {}
};

template <typename T>
using RawData = std::array<std::byte, sizeof(T)>;

template <typename T, bool HasSeq>
class Cell : public SeqField<HasSeq>, public IsConstructedField<HasSeq>
{
private:
    alignas(alignof(T)) RawData<T> data_;

public:
    Cell(size_t i) : SeqField<HasSeq>(i), IsConstructedField<HasSeq>(false) {}

    ~Cell() noexcept {
        if (this->is_constructed) destroy();
    }

    template <typename ...Args>
    void construct(Args &&...args) requires std::is_nothrow_constructible_v<T, Args&&...> {
        new (&data_) T(std::forward<Args>(args)...);

        this->is_constructed = true;
    }

    void destroy() noexcept {
        reinterpret_cast<T *>(&data_)->~T();

        this->is_constructed = false;
    }

    T &&move() noexcept {
        return reinterpret_cast<T &&>(data_);
    }
};

template <typename T, size_t N, typename SizeConstraint>
concept IsValidQueue = ZeroOrGreaterThanOne<N>
                    && IsValidSizeConstraint<SizeConstraint>
                    && OptionalPowerOfTwo<N, SizeConstraint>;

template <typename T, size_t N, typename SizeConstraint>
concept IsValidMpmcQueue = IsValidQueue<T, N, SizeConstraint> && FalseSharingSafe<Cell<T, true>>;

template <typename T, size_t N, typename SizeConstraint>
concept IsValidSpscQueue = IsValidQueue<T, N, SizeConstraint>;
 
template <typename DerivedImpl, bool UseSeq, typename T, size_t N = 0, IsValidSizeConstraint SizeConstraint = EnforcePowerOfTwo>
requires IsValidQueue<T, N, SizeConstraint>
class BaseQueue
{
private:
    using value_t = Cell<T, UseSeq>;

    using allocator_t = typename std::allocator<value_t>;

    using buffer_t = Buffer<value_t, SizeConstraint, N>;

protected:
    alignas(cache_line) buffer_t buffer_;

    alignas(cache_line) const size_t buffer_size_;

public:
    explicit BaseQueue(const size_t buffer_size = N, const allocator_t &allocator = allocator_t()) 
                                    : buffer_(buffer_size, allocator), buffer_size_(buffer_size) {

        if constexpr (N != 0) {
            assert(buffer_size == N
                && "Do not specify a constructor buffer size different from compile-time buffer size");
        }
        else {
            if (buffer_size <= 1) {
                throw std::invalid_argument("buffer_size should be greater than 1");
            }
            if (IsEnforcePowerOfTwo<SizeConstraint> && !details::is_power_of_two(buffer_size)) {
                throw std::invalid_argument("buffer_size should be a power of 2");
            }
        }

        static_assert(details::SizeMultipleOfCacheLine<BaseQueue>,
            "Queue size must be a multiple of a cache line to prevent false sharing with adjacent queues");

        for (size_t i = 0; i < buffer_size; i++) {
            new (&buffer_[i]) value_t(i);
        }
    }

    ~BaseQueue() noexcept {
        for (size_t i = 0; i < buffer_size_; i++) {
            buffer_[i].~Cell();
        }
    }

    BaseQueue(const BaseQueue &) = delete;
    BaseQueue &operator=(const BaseQueue &) = delete;

    template <typename... Args>
    void emplace(Args &&... args) noexcept
    requires std::is_nothrow_constructible_v<T, Args &&...> {
        static_cast<DerivedImpl*>(this)->emplace(std::forward<Args &&...>(args)...);
    }

    void push(const T &data) noexcept
    requires std::is_nothrow_copy_constructible_v<T> {
        static_cast<DerivedImpl*>(this)->emplace(data);
    }

    template <typename P>
    void push(P &&data) noexcept
    requires std::is_nothrow_constructible_v<T, P> {
        static_cast<DerivedImpl*>(this)->emplace(std::forward<P>(data));
    }

    template <typename... Args>
    [[nodiscard]] bool try_emplace(Args &&... args) noexcept
    requires std::is_nothrow_constructible_v<T, Args &&...> {
        return static_cast<DerivedImpl*>(this)->try_emplace(std::forward<Args &&...>(args)...);
    }

    [[nodiscard]] bool try_push(const T &data) noexcept
    requires std::is_nothrow_copy_constructible_v<T> {
        return static_cast<DerivedImpl*>(this)->try_emplace(data);
    }

    template <typename P>
    [[nodiscard]] bool try_push(P &&data) noexcept
    requires std::is_nothrow_constructible_v<T, P> {
        return static_cast<DerivedImpl*>(this)->try_emplace(std::forward<P>(data));
    }

    void pop(T &v) noexcept {
        static_cast<DerivedImpl*>(this)->pop(v);
    }
    
    [[nodiscard]] bool try_pop(T &v) noexcept {
        return static_cast<DerivedImpl*>(this)->try_pop(v);
    }

    /// Will return a negative value if there are one or more readers waiting on an empty queue
    [[nodiscard]] size_t size() const noexcept {
        return static_cast<const DerivedImpl*>(this)->get_enqueue_pos()
             - static_cast<const DerivedImpl*>(this)->get_dequeue_pos();
    }

    [[nodiscard]] size_t capacity() const noexcept {
        return buffer_size_;
    }

    [[nodiscard]] bool empty() const noexcept {
        return size() <= 0;
    }

};

static inline bool cas_add_one_relaxed(std::atomic<size_t> &atomic, size_t& val) noexcept {
    return atomic.compare_exchange_weak(
        val, val + 1, std::memory_order::relaxed, std::memory_order::relaxed);
}

} // namespace details

template <typename T, size_t N = 0, details::IsValidSizeConstraint SizeConstraint = EnforcePowerOfTwo>
requires details::IsValidMpmcQueue<T, N, SizeConstraint>
class MpmcQueue : public details::BaseQueue<MpmcQueue<T, N, SizeConstraint>, true, T, N, SizeConstraint>
{
private:
    using base_queue_t = details::BaseQueue<MpmcQueue<T, N, SizeConstraint>, true, T, N, SizeConstraint>;
    using allocator_t = typename std::allocator<details::Cell<T, true>>;

    alignas(details::cache_line) std::atomic<size_t> enqueue_pos_{0};
    alignas(details::cache_line) std::atomic<size_t> dequeue_pos_{0};

public:
    explicit MpmcQueue(const size_t buffer_size = N, const allocator_t &allocator = allocator_t()) 
                        : base_queue_t(buffer_size, allocator) {}

    template <typename... Args>
    void emplace(Args &&... args) noexcept
    requires std::is_nothrow_constructible_v<T, Args &&...> {
        size_t pos = enqueue_pos_.fetch_add(1, std::memory_order::relaxed);
        auto &cell = this->buffer_[pos];

        while (pos != cell.seq.load(std::memory_order::acquire));

        cell.construct(std::forward<Args>(args)...);
        cell.seq.store(pos+1, std::memory_order::release);
    }

    template <typename... Args>
    [[nodiscard]] bool try_emplace(Args &&... args) noexcept
    requires std::is_nothrow_constructible_v<T, Args &&...> {
        while (true) {
            size_t pos = enqueue_pos_.load(std::memory_order::relaxed);
            auto &cell = this->buffer_[pos];

            const size_t seq = cell.seq.load(std::memory_order::acquire);
            const int64_t diff = seq - pos;

            if (diff == 0 && details::cas_add_one_relaxed(enqueue_pos_, pos)) {
                cell.construct(std::forward<Args>(args)...);
                cell.seq.store(pos+1, std::memory_order::release);

                return true;
            }
            else if (diff < 0) {
                return false;
            }
        }
    }
    
    void pop(T &v) noexcept {
        const size_t pos = dequeue_pos_.fetch_add(1, std::memory_order::relaxed);
        auto &cell = this->buffer_[pos];

        while (pos + 1 != cell.seq.load(std::memory_order::acquire));

        v = cell.move();
        cell.destroy();

        cell.seq.store(pos + this->buffer_size_, std::memory_order::release);
    }

    [[nodiscard]] bool try_pop(T &v) noexcept {
        while (true) { 
            size_t pos = dequeue_pos_.load(std::memory_order::relaxed);
            auto &cell = this->buffer_[pos];

            const size_t seq = cell.seq.load(std::memory_order::acquire);
            const int64_t diff = seq - (pos + 1);

            if (diff == 0 && details::cas_add_one_relaxed(dequeue_pos_, pos)) {
                v = cell.move();
                cell.destroy();

                cell.seq.store(pos + this->buffer_size_, std::memory_order::release);

                return true;
            }
            else if (diff < 0) {
                return false;
            }
        }
    }

    [[nodiscard]] size_t get_enqueue_pos() const noexcept {
        return enqueue_pos_.load(std::memory_order::relaxed);
    }

    [[nodiscard]] size_t get_dequeue_pos() const noexcept {
        return dequeue_pos_.load(std::memory_order::relaxed);
    }
};

template <typename T, size_t N = 0, details::IsValidSizeConstraint SizeConstraint = EnforcePowerOfTwo>
requires details::IsValidSpscQueue<T, N, SizeConstraint>
class SpscQueue : public details::BaseQueue<SpscQueue<T, N, SizeConstraint>, false, T, N, SizeConstraint>
{
private:
    using base_queue_t = details::BaseQueue<SpscQueue<T, N, SizeConstraint>, false, T, N, SizeConstraint>;
    using allocator_t = typename std::allocator<details::Cell<T, false>>;

    alignas(details::cache_line) std::atomic<size_t> enqueue_pos_{0};
    size_t cached_enqueue_limit_ = 0;

    alignas(details::cache_line) std::atomic<size_t> dequeue_pos_{0};
    size_t cached_dequeue_limit_ = 0;

public:
    explicit SpscQueue(const size_t buffer_size = N, const allocator_t &allocator = allocator_t()) 
                        : base_queue_t(buffer_size, allocator) {}

    template <typename... Args>
    void emplace(Args &&... args) noexcept
    requires std::is_nothrow_constructible_v<T, Args &&...> {
        const size_t pos = enqueue_pos_.load(std::memory_order::relaxed);
        auto &cell = this->buffer_[pos];

        while (pos == cached_enqueue_limit_) {
            cached_enqueue_limit_ = dequeue_pos_.load(std::memory_order::acquire) + this->buffer_size_;
        }

        cell.construct(std::forward<Args>(args)...);
        enqueue_pos_.store(pos + 1, std::memory_order::release);
    }

    template <typename... Args>
    [[nodiscard]] bool try_emplace(Args &&... args) noexcept
    requires std::is_nothrow_constructible_v<T, Args &&...> {
        const size_t pos = enqueue_pos_.load(std::memory_order::relaxed);
        auto &cell = this->buffer_[pos];

        if (pos == cached_enqueue_limit_) {
            cached_enqueue_limit_ = dequeue_pos_.load(std::memory_order::acquire) + this->buffer_size_;

            if (pos == cached_enqueue_limit_) return false;
        }

        cell.construct(std::forward<Args>(args)...);
        enqueue_pos_.store(pos + 1, std::memory_order::release);

        return true;
    }
    
    void pop(T &v) noexcept {
        const size_t pos = dequeue_pos_.load(std::memory_order::relaxed);
        auto &cell = this->buffer_[pos];

        while (pos == cached_dequeue_limit_) {
            cached_dequeue_limit_ = enqueue_pos_.load(std::memory_order::acquire);
        }

        v = cell.move();
        cell.destroy();

        dequeue_pos_.store(pos + 1, std::memory_order::release);
    }

    [[nodiscard]] bool try_pop(T &v) noexcept {
        const size_t pos = dequeue_pos_.load(std::memory_order::relaxed);
        auto &cell = this->buffer_[pos];

        if (pos == cached_dequeue_limit_) {
            cached_dequeue_limit_ = enqueue_pos_.load(std::memory_order::acquire);

            if (pos == cached_dequeue_limit_) return false;
        }

        v = cell.move();
        cell.destroy();

        dequeue_pos_.store(pos + 1, std::memory_order::release);

        return true;
    }

    [[nodiscard]] size_t get_enqueue_pos() const noexcept {
        return enqueue_pos_.load(std::memory_order::relaxed);
    }

    [[nodiscard]] size_t get_dequeue_pos() const noexcept {
        return dequeue_pos_.load(std::memory_order::relaxed);
    }
};

} // namespace jdz