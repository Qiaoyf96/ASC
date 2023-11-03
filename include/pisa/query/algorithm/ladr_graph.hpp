#pragma once

#include <fstream>
#include <string>

namespace pisa {

std::vector<std::string> split (const std::string &s, char delim) {
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}

struct ladr_graph {

    std::unordered_map<int, std::vector<int>> graph;
    bool visited[4096];
    int counts[4096];

    std::set<int> pos_docs;

    int m_num_clusters;

    explicit ladr_graph(const std::string& filename, int num) : m_num_clusters(num) {
        ifstream file(filename);

        std::string line;
        while (getline(file, line)) {
            auto v = split(line, '\t');
            auto clusters = split(v[1], ' ');
            for (auto &c : clusters) {
                graph[stoi(v[0])].push_back(stoi(c));
            }
        }
    }

    void init() {
        memset(visited, 0, sizeof(visited));
        memset(counts, 0, sizeof(counts));
        pos_docs.clear();
    }

    void init_counts() {
	memset(counts, 0, sizeof(counts));
	pos_docs.clear();
    }

    void add(int docid) {
        if (pos_docs.find(docid) == pos_docs.end()) {
            for (auto &d : graph[docid]) {
                counts[d] += 1;
            }
            pos_docs.insert(docid);
        }
    }

    int next_cluster() {
        int ret = -1;
        int max_count = 0;
        for (int c = 0; c < 4096; c++) {
            if (!visited[c] && counts[c] > max_count) {
                max_count = counts[c];
                ret = c;
            }
        }
        
        return ret;
    }

    bool is_visited(int clusterid) {
        return visited[clusterid];
    }

    void visit(int clusterid) {
        // printf("next c: %d\n", clusterid);
        visited[clusterid] = 1;
    }

};

}  // namespace pisa
