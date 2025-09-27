#include <iostream>
#include <vector>
#include <memory>
#include <thread>
#include <mutex>
#include <chrono>
#include <cstdlib>
#include <algorithm>

using namespace std;

struct Node
{
    int key;
    int value;
    unique_ptr<Node> left;
    unique_ptr<Node> right;
    Node(int k, int v) : key(k), value(v) {}
};

void bst_insert(unique_ptr<Node> &root, int key, int value)
{
    if (!root)
    {
        root = make_unique<Node>(key, value);
        return;
    }
    if (key == root->key)
        root->value = value;
    else if (key < root->key)
        bst_insert(root->left, key, value);
    else
        bst_insert(root->right, key, value);
}

unique_ptr<int> bst_search(const Node *root, int key)
{
    const Node *cur = root;
    while (cur)
    {
        if (key == cur->key)
            return make_unique<int>(cur->value);
        else if (key < cur->key)
            cur = cur->left.get();
        else
            cur = cur->right.get();
    }
    return nullptr;
}

class MultiThreadedBST
{
private:
    unique_ptr<Node> root;
    int num_threads;
    mutex tree_mutex;

public:
    MultiThreadedBST(int threads = thread::hardware_concurrency())
    {
        int max_threads = thread::hardware_concurrency();
        num_threads = threads > max_threads ? max_threads : threads;
    }

    void insert(const vector<pair<int, int>> &items)
    {
        int n = items.size();
        if (n == 0)
            return;
        int chunk_size = (n + num_threads - 1) / num_threads;

        vector<thread> threads;
        for (int t = 0; t < num_threads; ++t)
        {
            int start = t * chunk_size;
            int end = min((t + 1) * chunk_size, n);
            threads.emplace_back([this, &items, start, end]()
                                 {
                for (int i = start; i < end; ++i) {
                    lock_guard<mutex> lock(tree_mutex); 
                    bst_insert(root, items[i].first, items[i].second);
                } });
        }
        for (auto &th : threads)
            th.join();
    }

    void search(const vector<int> &keys, vector<unique_ptr<int>> &results)
    {
        int n = keys.size();
        results.resize(n);
        int chunk_size = (n + num_threads - 1) / num_threads;

        vector<thread> threads;
        for (int t = 0; t < num_threads; ++t)
        {
            int start = t * chunk_size;
            int end = min((t + 1) * chunk_size, n);
            threads.emplace_back([this, &keys, &results, start, end]()
                                 {
                for (int i = start; i < end; ++i)
                    results[i] = bst_search(root.get(), keys[i]); });
        }
        for (auto &th : threads)
            th.join();
    }
};

class SingleThreadedBST
{
private:
    unique_ptr<Node> root;

public:
    void insert(const vector<pair<int, int>> &items)
    {
        for (auto &p : items)
            bst_insert(root, p.first, p.second);
    }

    void search(const vector<int> &keys, vector<unique_ptr<int>> &results)
    {
        results.resize(keys.size());
        for (int i = 0; i < keys.size(); ++i)
            results[i] = bst_search(root.get(), keys[i]);
    }
};

void benchmark(int n, int num_threads)
{
    vector<pair<int, int>> data(n);
    vector<int> keys(n);

    for (int i = 0; i < n; ++i)
    {
        data[i] = {rand(), rand()};
        keys[i] = rand();
    }

    SingleThreadedBST sst;
    MultiThreadedBST mst(num_threads);

    auto start = chrono::high_resolution_clock::now();
    sst.insert(data);
    auto end = chrono::high_resolution_clock::now();
    cout << "Single-thread insert = "
         << chrono::duration<double>(end - start).count() << " s\n";

    start = chrono::high_resolution_clock::now();
    mst.insert(data);
    end = chrono::high_resolution_clock::now();
    cout << "Multi-thread insert = "
         << chrono::duration<double>(end - start).count() << " s\n";

    vector<unique_ptr<int>> res1, res2;

    start = chrono::high_resolution_clock::now();
    sst.search(keys, res1);
    end = chrono::high_resolution_clock::now();
    cout << "Single-thread search = "
         << chrono::duration<double>(end - start).count() << " s\n";

    start = chrono::high_resolution_clock::now();
    mst.search(keys, res2);
    end = chrono::high_resolution_clock::now();
    cout << "Multi-thread search = "
         << chrono::duration<double>(end - start).count() << " s\n";
}

int main()
{
    MultiThreadedBST tree(2);
    tree.insert({{3, 24}, {5, 55}, {15, 152}, {7, 70}});

    vector<unique_ptr<int>> results;
    vector<int> search_keys = {3, 7, 10, 15, 100};
    tree.search(search_keys, results);

    cout << "Search results:\n";
    for (size_t i = 0; i < results.size(); ++i)
    {
        if (results[i])
            cout << "Key " << search_keys[i] << " = " << *results[i] << "\n";
        else
            cout << "Key " << search_keys[i] << " not found\n";
    }

    vector<int> sizes = {10000, 50000, 100000, 1000000, 10000000};
    vector<int> threads = {1, 2, 4, 8, 16};

    for (int n : sizes)
    {
        for (int t : threads)
        {
            cout << "n=" << n << ", threads=" << t << "\n";
            benchmark(n, t);
            cout << "\n";
        }
    }
}
