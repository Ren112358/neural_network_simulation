#include <iostream>
#include <vector>
#include <cmath>

#include "neuron.h"

using namespace std;


class Neural_Network {
public:
    Neural_Network(
            vector<Neuron> neurons,
            vector< vector<double> > g_syn,
            vector< vector<double> > V_syn
            );
    void update_states(vector<double> I_ext, double time);
    vector<double> get_V(); // to record potentials
private:
    vector<Neuron> neurons;
    vector< vector<double> > g_syn;
    vector< vector<double> > V_syn;
    vector<double> t_f;  // arrivial time of a synaptic action potentials
    double tau;  // time constant
    double fire_threshold;  // threshold of neuron fire
    double heaviside_func(double x);
    vector<double> calc_synaptic_current(double time);
};


Neural_Network::Neural_Network(
        vector<Neuron> neurons,
        vector< vector<double> > g_syn,
        vector< vector<double> > V_syn
        ) {
    this->neurons = neurons;
    this->g_syn = g_syn;
    this->V_syn = V_syn;
    this->t_f = vector<double>(this->neurons.size(), 0);
    this->tau = 5;
    this->fire_threshold = 50;
}

double Neural_Network::heaviside_func(double x) {
    if (x > 0) {
        return 1.0;
    }
    else {
        return 0.0;
    }
}

vector<double> Neural_Network::calc_synaptic_current(double time) {
    vector<double> I_syn;
    vector<double> V = this->get_V();
    for (long long i = 0; i < this->neurons.size(); i++) {
        double I_tmp = 0;
        for (long long j = 0; j < this->neurons.size(); j++) {
            double arg1 = -(time - this->t_f[j]) / this->tau;
            double arg2 = time - this->t_f[j];
            double g_tmp = this->g_syn[i][j] * exp(arg1) * this->heaviside_func(arg2);
            I_tmp += g_tmp * (V[i] - this->V_syn[i][j]);
        }
        I_syn.push_back(I_tmp);
    }
    return I_syn;
}

void Neural_Network::update_states(vector<double> I_ext, double time) {
    vector<double> I_syn = calc_synaptic_current(time);
    for (long long i = 0; i < this->neurons.size(); i++) {
        this->neurons[i].update_state(I_ext[i], I_syn[i]);
    }
    vector<double> V = this->get_V();
    for (long long i = 0; i < V.size(); i++) {
        if (V[i] > this->fire_threshold) {
            this->t_f[i] = time;
        }
    }
}

vector<double> Neural_Network::get_V() {
    vector<double> vec_V;
    for (long long i = 0; i < this->neurons.size(); i++) {
        vec_V.push_back(this->neurons[i].get_V());
    }
    return vec_V;
}
