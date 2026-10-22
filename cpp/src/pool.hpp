// Parallel infrastructure: an ordered, streaming thread pool over sequences.
//
//   run_ordered<Item>(read_batch, work, emit, jobs)
//     read_batch(std::vector<Item>&)  fills the next batch (left empty = end of input); calling thread
//     work(Item&)                     the per-sequence function; worker threads, any order
//     emit(std::vector<Item>&)        each finished batch, on the calling thread, always in input order
//
// Batch k+1 is read while batch k is being worked on, so at most two batches are in memory at once:
// memory is bounded by the batch size (the caller's choice, e.g. ~2 M nucleotides), not by the input.
// jobs <= 0 means one thread per hardware thread. An exception from work() stops the batch and is
// rethrown on the calling thread (no std::terminate from a worker); nothing of that batch is emitted.
#pragma once

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace coal {

inline int resolve_jobs(int jobs) {
    return jobs > 0 ? jobs : (int)std::max(1u, std::thread::hardware_concurrency());
}

// work(item) for every item of `batch` on up to `nthreads` threads; returns the number of threads used.
template <class Item, class Work>
int parallel_batch(std::vector<Item>& batch, Work& work, int nthreads) {
    const size_t n = batch.size();
    const int k = (int)std::min<size_t>((size_t)std::max(1, nthreads), n);
    if (k <= 1) {
        for (Item& it : batch) work(it);
        return n ? 1 : 0;
    }
    std::atomic<size_t> next{0};
    std::atomic<bool> stop{false};
    std::exception_ptr err;
    std::mutex m;
    auto body = [&]() {
        size_t i;
        while (!stop.load() && (i = next++) < n) {
            try {
                work(batch[i]);
            } catch (...) {
                std::lock_guard<std::mutex> g(m);
                if (!err) err = std::current_exception();
                stop = true;
            }
        }
    };
    std::vector<std::thread> threads;
    threads.reserve((size_t)k);
    try {
        for (int t = 0; t < k; ++t) threads.emplace_back(body);
    } catch (...) {  // could not start a thread: stop the others, then report
        stop = true;
        for (auto& th : threads) th.join();
        throw;
    }
    for (auto& th : threads) th.join();
    if (err) std::rethrow_exception(err);
    return k;
}

struct PoolStats {
    size_t batches = 0, items = 0;
    size_t max_live = 0;  // most items held at once (current + next batch)
    int max_threads = 0;  // most worker threads used for one batch
};

template <class Item, class Read, class Work, class Emit>
PoolStats run_ordered(Read&& read_batch, Work&& work, Emit&& emit, int jobs) {
    PoolStats st;
    const int nthreads = resolve_jobs(jobs);
    std::vector<Item> cur, nxt;
    read_batch(cur);
    while (!cur.empty()) {
        int used = 0;
        std::exception_ptr werr;
        std::thread worker([&]() {
            try {
                used = parallel_batch(cur, work, nthreads);
            } catch (...) {
                werr = std::current_exception();
            }
        });
        try {
            read_batch(nxt);
        } catch (...) {
            worker.join();
            throw;
        }
        worker.join();
        if (werr) std::rethrow_exception(werr);
        st.max_live = std::max(st.max_live, cur.size() + nxt.size());
        st.max_threads = std::max(st.max_threads, used);
        ++st.batches;
        st.items += cur.size();
        emit(cur);
        cur.swap(nxt);
        nxt.clear();
    }
    return st;
}

}  // namespace coal
