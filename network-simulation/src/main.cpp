#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "consts.h"
#include "neural_network.h"

using namespace std;

vector<string> split(string& input, char delimiter)
{
    istringstream stream(input);
    string field;
    vector<string> result;
    while (getline(stream, field, delimiter)) {
        result.push_back(field);
    }
    return result;
}

int main() {
    // class has data paths
    Data_path data_path;

    // read the data of sampling time
    ifstream ifs_time(data_path.time);
    vector<double> time;
    string str_time;
    while(ifs_time >> str_time) {
        time.push_back(stod(str_time));
    }

    // read the data of input current
    ifstream ifs_input_currents(data_path.input_currents);
    vector< vector<double> > I_ext;
    string line;
    while (getline(ifs_input_currents, line)) {
        vector<string> str_vec = split(line, ',');
        vector<double> double_vec(str_vec.size());
        for (int i = 0; i < double_vec.size(); i++) {
            double_vec[i] = stod(str_vec[i]);
        }
        I_ext.push_back(double_vec);
    }

    long long N = I_ext[0].size();  // #(neurons)
    vector<Neuron> neurons(N);
    for (long long i = 0; i < N; i++) {
        Neuron neuron;
        neurons[i] = neuron;
    }

    // read the data of synaptic conductances
    ifstream ifs_synaptic_conductances(data_path.synaptic_conductances);
    vector< vector<double> > g_syn;
    while (getline(ifs_synaptic_conductances, line)) {
        vector<string> str_vec = split(line, ',');
        vector<double> double_vec(str_vec.size());
        for (long long i = 0; i < double_vec.size(); i++) {
            double_vec[i] = stod(str_vec[i]);
        }
        g_syn.push_back(double_vec);
    }

    // read the data of synaptic equilibrium potentails
    ifstream ifs_synaptic_equilibrium_potentials(data_path.synaptic_equilibrium_potentials);
    vector< vector<double> > V_syn;
    while (getline(ifs_synaptic_equilibrium_potentials, line)) {
        vector<string> str_vec = split(line, ',');
        vector<double> double_vec(str_vec.size());
        for (long long i = 0; i < double_vec.size(); i++) {
            double_vec[i] = stod(str_vec[i]);
        }
        V_syn.push_back(double_vec);
    }

    Neural_Network network(neurons, g_syn, V_syn);

    long long T = I_ext.size();  // #(sampling)
    string path_potentials = data_path.potentials;
    ofstream ofs_potential(path_potentials);
    for (long long t = 0; t < T; t++) {
        network.update_states(I_ext[t], t);
        vector<double> vec_V = network.get_V();
        for (long long i = 0; i < N; i++) {
            ofs_potential << fixed << setprecision(5) << vec_V[i];
            if (i < N - 1) {
                ofs_potential << ",";
            }
            else {
                ofs_potential << endl;
            }
        }
    }
    ofs_potential.close();

    return 0;
}
