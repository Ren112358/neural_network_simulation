#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>

#include "consts.h"

using namespace std;

// want to use dt in neuron.h

int main() {
    // consts
    Simulation_consts sim_consts;
    Synaptic_consts syn_consts;
    Data_path data_path;

    const double dt = sim_consts.dt;  // sampling interval
    const long long N = sim_consts.N;  // #(neuron)
    const long long T = sim_consts.T;  // #(sampling)

    // make and save samplig time
    double time = 0;
    ofstream ofs_time(data_path.time);
    for (long long i = 0; i < T; i++) {
        ofs_time << fixed << setprecision(5) << time << endl;
        time += dt;
    }
    ofs_time.close();

    // make input current
    static const double I_ext_max = sim_consts.I_ext_max;
    vector< vector<double> > I_ext;  // max input current
    for (long long i = 0; i < T; i++) {
        vector<double> tmp;
        for (long long j = 0; j < N; j++) {
            if (j == 0) {
                tmp.push_back(I_ext_max);
            }
            else {
                tmp.push_back(0);
            }
        }
        I_ext.push_back(tmp);
    }
    /*
    static const double PI = 3.141592653589793;
    double I_ext_max = 20.0;
    vector< vector<double> > I_ext;  // max input current
    for (long long i = 0; i < T; i++) {
        vector<double> tmp;
        for (long long j = 0; j < N; j++) {
            double numerator = 2 * PI * i;
            double denominator = T;
            double x = (numerator / denominator) - (j * PI / N);
            if (abs(sin(x)) > (1 / sqrt(2))) {
                tmp.push_back(I_ext_max);
            }
            else {
                tmp.push_back(0);
            }
        }
        I_ext.push_back(tmp);
    }
    */

    // save input current
    ofstream ofs_input_currents(data_path.input_currents);
    for (long long i = 0; i < T; i++) {
        for (long long j = 0; j < N; j++) {
            ofs_input_currents << fixed << setprecision(5) << I_ext[i][j];
            if (j < N - 1) {
                ofs_input_currents << ",";
            }
            else {
                ofs_input_currents << endl;
            }
        }
    }
    ofs_input_currents.close();

    // make synaptic conductance matrix and synaptic equilibrium potential matrix
    random_device rnd;
    mt19937 mt(rnd());
    uniform_int_distribution<> uniform_rand(0, 100);  // rand num [0, 100]
    vector< vector<double> > g_syn(N, vector<double>(N, 0.0));
    vector< vector<double> > V_syn(N, vector<double>(N, 0.0));

    for (long long i = 0; i < N; i++) {
        for (long long j = 0; j < N; j++) {
            if (uniform_rand(mt) < syn_consts.connection_ratio) {
                if (uniform_rand(mt) < syn_consts.excitatory_ratio) {
                    g_syn[i][j] = syn_consts.g_exc;
                    V_syn[i][j] = syn_consts.V_exc;
                }
                else {
                    g_syn[i][j] = syn_consts.g_inh;
                    V_syn[i][j] = syn_consts.V_inh;
                }
            }
        }
    }

    // no self connection
    for (long long i = 0; i < N; i++) {
        g_syn[i][i] = 0.0;
        V_syn[i][i] = 0.0;
    }

    // save synaptic conductance matrix
    ofstream ofs_synaptic_conductances(data_path.synaptic_conductances);
    for (long long i = 0; i < N; i++) {
        for (long long j = 0; j < N; j++) {
            ofs_synaptic_conductances << fixed << setprecision(5) << g_syn[i][j];
            if (j < N - 1) {
                ofs_synaptic_conductances << ",";
            }
            else {
                ofs_synaptic_conductances << endl;
            }
        }
    }
    ofs_synaptic_conductances.close();

    // save synaptic equilibrium potential matrix
    ofstream ofs_synaptic_equilibrium_potentials(data_path.synaptic_equilibrium_potentials);
    for (long long i = 0; i < N; i++) {
        for (long long j = 0; j < N; j++) {
            ofs_synaptic_equilibrium_potentials << fixed << setprecision(5) << V_syn[i][j];
            if (j < N - 1) {
                ofs_synaptic_equilibrium_potentials << ",";
            }
            else {
                ofs_synaptic_equilibrium_potentials << endl;
            }
        }
    }
    ofs_synaptic_equilibrium_potentials.close();

    return 0;
}
