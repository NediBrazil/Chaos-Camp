#if _WIN64
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <cstdlib>
#include <algorithm>
#include <random>
#include <cassert>
#include <cstdio>
#include <cstring>

#if _WIN64
#define bassert(test) (!!(test) ? (void)0 : ((void)printf("-- Assertion failed at line %d: %s\n", __LINE__, #test), __debugbreak()))
#else
#define bassert assert
#endif

const int NOT_FOUND = -1;
const int NOT_SEARCHED = -2;

namespace
{

#include <stdint.h>

#if __linux__ != 0
#include <time.h>
    static uint64_t timer_nsec()
    {
#if defined(CLOCK_MONOTONIC_RAW)
        const clockid_t clockid = CLOCK_MONOTONIC_RAW;
#else
        const clockid_t clockid = CLOCK_MONOTONIC;
#endif
        timespec t;
        clock_gettime(clockid, &t);
        return uint64_t(t.tv_sec) * 1000000000ULL + uint64_t(t.tv_nsec);
    }

#elif defined(_WIN64) || defined(_WIN32)
#define NOMINMAX
#include <Windows.h>
    static struct TimerBase
    {
        LARGE_INTEGER freq;
        TimerBase() { QueryPerformanceFrequency(&freq); }
    } timerBase;

    static uint64_t timer_nsec()
    {
        LARGE_INTEGER t;
        QueryPerformanceCounter(&t);
        return 1000000000ULL * t.QuadPart / timerBase.freq.QuadPart;
    }

#elif __APPLE__ != 0
#include <mach/mach_time.h>
    static struct TimerBase
    {
        mach_timebase_info_data_t tb;
        TimerBase() { mach_timebase_info(&tb); }
    } timerBase;

    static uint64_t timer_nsec()
    {
        const uint64_t t = mach_absolute_time();
        return t * timerBase.tb.numer / timerBase.tb.denom;
    }

#else
    static uint64_t timer_nsec() { return 0; }
#endif

    template <typename T>
    T *alignedAlloc(size_t count, void *&unaligned)
    {
        const size_t bytes = count * sizeof(T);
        unaligned = malloc(bytes + 63);
        if (!unaligned)
            return nullptr;
        T *const aligned = reinterpret_cast<T *>((reinterpret_cast<uintptr_t>(unaligned) + 63) & ~uintptr_t(63));
        return aligned;
    }

    template <typename T>
    struct AlignedArrayPtr
    {
        void *allocated = nullptr;
        T *aligned = nullptr;
        int64_t count = -1;

        AlignedArrayPtr() = default;
        AlignedArrayPtr(int64_t count) { init(count); }
        void init(int64_t newCount)
        {
            bassert(newCount > 0);
            free(allocated);
            aligned = alignedAlloc<T>(newCount, allocated);
            count = newCount;
        }
        void memset(int value)
        {
            ::memset(aligned, value, sizeof(T) * count);
        }
        ~AlignedArrayPtr() { free(allocated); }
        T *get() { return aligned; }
        const T *get() const { return aligned; }
        operator T *() { return aligned; }
        operator const T *() const { return aligned; }
        int64_t getCount() const { return count; }
        const T *begin() const { return aligned; }
        const T *end() const { return aligned + count; }
        int operator[](int index) const { return aligned[index]; }
        int &operator[](int index) { return aligned[index]; }
        AlignedArrayPtr(const AlignedArrayPtr &) = delete;
        AlignedArrayPtr &operator=(const AlignedArrayPtr &) = delete;
    };

    typedef AlignedArrayPtr<int> AlignedIntArray;

    const char magic[] = ".BSEARCH";
    const int magicSize = sizeof(magic) - 1;

    bool storeToFile(const AlignedIntArray &hayStack, const AlignedIntArray &needles, const char *name)
    {
        FILE *file = fopen(name, "wb+");
        if (!file)
            return false;
        const int64_t sizes[2] = {hayStack.getCount(), needles.getCount()};
        fwrite(magic, 1, magicSize, file);
        fwrite(sizes, 1, sizeof(sizes), file);
        fwrite(hayStack.get(), sizeof(int), hayStack.getCount(), file);
        fwrite(needles.get(), sizeof(int), needles.getCount(), file);
        fclose(file);
        return true;
    }

    bool loadFromFile(AlignedIntArray &hayStack, AlignedIntArray &needles, const char *name)
    {
        FILE *file = fopen(name, "rb");
        if (!file)
            return false;
        char test[magicSize] = {0};
        int64_t sizes[2];
        int allOk = true;
        allOk &= magicSize == fread(test, 1, magicSize, file);
        if (strncmp(magic, test, magicSize))
        {
            printf("Bad magic constant in file [%s]\n", name);
            return false;
        }
        allOk &= sizeof(sizes) == fread(sizes, 1, sizeof(sizes), file);
        hayStack.init(sizes[0]);
        needles.init(sizes[1]);
        allOk &= hayStack.getCount() == int64_t(fread(hayStack.get(), sizeof(int), hayStack.getCount(), file));
        allOk &= needles.getCount() == int64_t(fread(needles.get(), sizeof(int), needles.getCount(), file));
        fclose(file);
        return allOk;
    }

    int verify(const AlignedIntArray &hayStack, const AlignedIntArray &needles, const AlignedIntArray &indices)
    {
        for (int c = 0; c < needles.getCount(); c++)
        {
            const int value = needles[c];
            int left = 0;
            int right = hayStack.getCount() - 1;
            int result = -1;
            while (left <= right)
            {
                int mid = left + (right - left) / 2;
                if (hayStack[mid] < value)
                {
                    left = mid + 1;
                }
                else if (hayStack[mid] > value)
                {
                    right = mid - 1;
                }
                else
                {
                    result = mid;
                    break;
                }
            }

            if (result == -1)
            {
                if (indices[c] != NOT_FOUND)
                    return c;
            }
            else
            {
                if (indices[c] != result)
                    return c;
            }
        }
        return -1;
    }
}

struct StackAllocator
{
    StackAllocator(uint8_t *ptr, int bytes)
        : totalBytes(bytes), data(ptr) {}

    template <typename T>
    T *alloc(int count)
    {
        const int size = count * sizeof(T);
        if (idx + size > totalBytes)
            return nullptr;
        uint8_t *start = data + idx;
        idx += size;
        return reinterpret_cast<T *>(start);
    }

    void freeAll() { idx = 0; }
    int maxBytes() const { return totalBytes; }
    int freeBytes() const { return totalBytes - idx; }

private:
    const int totalBytes;
    int idx = 0;
    uint8_t *data = nullptr;
};

static void binarySearch(const AlignedIntArray &hayStack, const AlignedIntArray &needles, AlignedIntArray &indices)
{
    for (int i = 0; i < needles.getCount(); i++)
    {
        int left = 0;
        int right = hayStack.getCount() - 1;
        int found = -1;
        while (left <= right)
        {
            int mid = left + (right - left) / 2;
            if (hayStack[mid] == needles[i])
            {
                found = mid;
                break;
            }
            else if (hayStack[mid] < needles[i])
            {
                left = mid + 1;
            }
            else
            {
                right = mid - 1;
            }
        }
        indices[i] = found;
    }
}

static void betterSearch(const AlignedIntArray &hayStack, const AlignedIntArray &needles, AlignedIntArray &indices, StackAllocator &allocator)
{
    const int n = hayStack.getCount();
    const int blockSize = 64;
    const int numBlocks = (n + blockSize - 1) / blockSize;

    struct Block
    {
        int minVal;
        int maxVal;
        int start;
        int end;
    };

    Block *blocks = allocator.alloc<Block>(numBlocks);
    if (!blocks)
    {
        binarySearch(hayStack, needles, indices);
        return;
    }
    for (int b = 0; b < numBlocks; b++)
    {
        int start = b * blockSize;
        int end = std::min(start + blockSize, n);
        blocks[b].start = start;
        blocks[b].end = end - 1;
        blocks[b].minVal = hayStack[start];
        blocks[b].maxVal = hayStack[end - 1];
    }

    for (int i = 0; i < needles.getCount(); i++)
    {
        int needle = needles[i];
        int found = -1;
        int L = 0, R = numBlocks - 1;

        while (L <= R)
        {
            int MID = (L + R) / 2;
            if (needle < blocks[MID].minVal)
            {
                R = MID - 1;
            }
            else if (needle > blocks[MID].maxVal)
            {
                L = MID + 1;
            }
            else
            {
                int left = blocks[MID].start;
                int right = blocks[MID].end;
                while (left <= right)
                {
                    int mid = left + (right - left) / 2;
                    if (hayStack[mid] == needle)
                    {
                        found = mid;
                        break;
                    }
                    else if (hayStack[mid] < needle)
                    {
                        left = mid + 1;
                    }
                    else
                    {
                        right = mid - 1;
                    }
                }
                if (found != -1)
                    break;
            }
        }

        indices[i] = found;
    }
}

int main()
{
    printf("+ Correctness tests ... \n");
    const int heapSize = 1 << 13;
    const int64_t searches = 400ll * (1 << 26);
    int testCaseCount = 0;
    for (int r = 0;; r++)
    {
        AlignedIntArray hayStack;
        AlignedIntArray needles;
        char fname[64] = {0};
        snprintf(fname, sizeof(fname), "%d.bsearch", r);
        if (!loadFromFile(hayStack, needles, fname))
            break;
        printf("Checking %s... ", fname);
        AlignedIntArray indices(needles.getCount());
        AlignedArrayPtr<uint8_t> heap(heapSize);
        StackAllocator allocator(heap, heapSize);
        indices.memset(NOT_SEARCHED);
        betterSearch(hayStack, needles, indices, allocator);
        if (verify(hayStack, needles, indices) != -1)
            return -1;
        indices.memset(NOT_SEARCHED);
        binarySearch(hayStack, needles, indices);
        if (verify(hayStack, needles, indices) != -1)
            return -1;
        printf("OK\n");
        ++testCaseCount;
    }
    printf("+ Speed tests ... \n");
    for (int r = 0; r < testCaseCount; r++)
    {
        AlignedIntArray hayStack;
        AlignedIntArray needles;
        char fname[64] = {0};
        snprintf(fname, sizeof(fname), "%d.bsearch", r);
        if (!loadFromFile(hayStack, needles, fname))
        {
            printf("Failed to load %s for speed test, continuing\n", fname);
            continue;
        }
        const int testRepeat = std::min<int64_t>(1000ll, searches / hayStack.getCount());
        printf("Running speed test for %s, %d repeats \n", fname, testRepeat);
        AlignedIntArray indices(needles.getCount());
        AlignedArrayPtr<uint8_t> heap(heapSize);
        StackAllocator allocator(heap, heapSize);
        uint64_t t0, t1;
        {
            indices.memset(NOT_SEARCHED);
            t0 = timer_nsec();
            for (int test = 0; test < testRepeat; ++test)
                binarySearch(hayStack, needles, indices);
            t1 = timer_nsec();
        }
        const double totalBinary = (double(t1 - t0) * 1e-9) / testRepeat;
        printf("\tbinarySearch time %f\n", totalBinary);
        {
            indices.memset(NOT_SEARCHED);
            t0 = timer_nsec();
            for (int test = 0; test < testRepeat; ++test)
                betterSearch(hayStack, needles, indices, allocator);
            t1 = timer_nsec();
        }
        const double totalBetter = (double(t1 - t0) * 1e-9) / testRepeat;
        printf("\tbetterSearch time %f\n", totalBetter);
        if (totalBetter < totalBinary)
            printf("Great success!\n");
    }
    return 0;
}
