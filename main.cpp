#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
using namespace std;

class Vertex {
public:
    vector<Vertex*> adjacent_vertices;
    int seq_num;
    int position;
};

multimap<string, Vertex*> graph;
vector<pair<string, vector<pair<int, char>>>> sequences;
vector<pair<string, vector<pair<int, int>>>> quals;

void getFastaFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Could not open file " << filename << endl;
        return;
    }

    string line;
    pair<string, vector<pair<int, char>>> current;
    int counter = 1;
    while (getline(file, line)) {
        if (line[0] == '>') {
            if (!current.first.empty()) {
                sequences.emplace_back(current);
                current = pair<string, vector<pair<int, char>>>();
                counter = 1;
            }
            current.first = line.substr(1);
        } else {
            for(auto it: line) {
                current.second.emplace_back(pair{counter, it});
                counter++;
            }
        }
    }
    sequences.emplace_back(current);
}

void getQualsFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Could not open file " << filename << endl;
        return;
    }
    string line;
    pair<string, vector<pair<int, int>>> current;
    int counter = 1;

    while (getline(file, line)) {
        if (line[0] == '>') {
            if (!current.first.empty()) {
                quals.push_back(current);
                current = pair<string, vector<pair<int, int>>>();
                counter = 1;
            }
            current.first = line.substr(1);
        } else {
            istringstream iss(line);
            int temp;
            while(iss >> temp) {
                current.second.emplace_back(pair{counter, temp});
                counter++;
            }
        }
    }

    if (!current.first.empty()) {
        quals.push_back(current);
    }
}

void remove_nucleotides(int threshold) {
    int count = 0;
    for (auto it : quals) {
        int offset = 0;
        for(auto pairs : it.second) {
            if(pairs.second < threshold) {
                sequences.at(count).second.erase(sequences.at(count).second.begin()+pairs.first - 1 - offset);
                quals.at(count).second.erase(quals.at(count).second.begin()+pairs.first -1 - offset);
                offset++;
            }
        }
        count++;
    }
}

multimap<string, Vertex*> create_graph(int k) {
    //Tworzenie wierzchołków
    multimap<string, Vertex*> new_graph;
    int count_seq = 1;
    for(auto it: sequences) {
        if(k<=it.second.size()) {
            for(int n = 0; n < it.second.size() - k + 1; n++) {
                Vertex* new_vertex = new Vertex();
                string oligo;
                for(int o = 0; o < k; o++) {
                    oligo += it.second[n+o].second;
                }
                new_vertex->position = it.second[n].first;
                new_vertex->seq_num = count_seq;
                new_graph.emplace(oligo , new_vertex);
            }
        }

        count_seq++;
    }
    //Łączenie wierzchołków
    for (auto it = new_graph.begin(); it != new_graph.end(); ) {
        string oligo = it->first;
        auto range = new_graph.equal_range(oligo);
        vector<Vertex*> vertices;
        for (auto itr = range.first; itr != range.second; ++itr) {
            vertices.push_back(itr->second);
        }
        for (size_t i = 0; i < vertices.size(); ++i) {
            for (size_t j = i + 1; j < vertices.size(); ++j) {
                if (vertices[i]->seq_num != vertices[j]->seq_num && abs(vertices[i]->position - vertices[j]->position) < k*10) {
                    vertices[i]->adjacent_vertices.push_back(vertices[j]);
                    vertices[j]->adjacent_vertices.push_back(vertices[i]);
                }
            }
        }

        it = range.second;
    }

    return new_graph;
}

void output() {

    for (auto it : sequences) {
        cout<<it.first<<" "<<endl;
        for(auto pairs : it.second) {
            cout << pairs.first << "."<< pairs.second<<" ";
        }
        cout<<endl;
    }
    for (auto it : quals) {
        cout<<it.first<<" "<<endl;
        for(auto pairs : it.second) {
            cout << pairs.first << "."<< pairs.second<<" ";
        }
        cout<<endl<<"GRAPH:"<<endl;
    }
    for(auto it : graph) {
        cout<< it.first<<" seq:"<< it.second->seq_num<<" pos: "<<it.second->position<<endl;
    }
}

void print_graph(const multimap<string, Vertex*>& graph) {
    for (auto it = graph.begin(); it != graph.end(); ++it) {
        string oligo = it->first;
        Vertex* vertex = it->second;

        // Wypisz szczegóły o wierzchołku
        cout << oligo
             << " | Seq: " << vertex->seq_num
             << " | Position: " << vertex->position
             << " | Adjacent: ";

        // Wypisz wszystkie wierzchołki połączone z tym wierzchołkiem
        for (Vertex* adjacent : vertex->adjacent_vertices) {
            cout << "(Seq: " << adjacent->seq_num
                 << ", Pos: " << adjacent->position << ") ";
        }
        cout << endl;
    }
}

int main() {
    getFastaFromFile("C:/Users/micha/CLionProjects/MotivSearching/FastaInput");
    getQualsFromFile("C:/Users/micha/CLionProjects/MotivSearching/QualInput");

    remove_nucleotides(32);
    graph = create_graph(7);

    //output();
    print_graph(graph);
    return 0;
}
