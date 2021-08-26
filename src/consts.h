#ifndef CONSTS_H
#define CONSTS_H

#include<string>

using namespace std;

class Simulation_consts {
    public:
        Simulation_consts();
        double dt;  // sampling time interval
        double N;  // #(neuron)
        long long T;  // #(sampling)
        double I_ext_max;  // max input current
};

class Synaptic_consts {
    public:
        Synaptic_consts();
        int connection_ratio;  // range [0, 100]
        int excitatory_ratio;  // excitatory for inhibitory ratio. range [0, 100]
        double g_exc;  // excitatory synaptic conductance
        double g_inh;  // inhibitory synaptic conductance
        double V_exc;  // excitatory synaptic equilibrium potential
        double V_inh;  // inhibitory synaptic equilibrium potential
};

class Data_path {
    public:
        Data_path();
        string time;
        string input_currents;
        string synaptic_equilibrium_potentials;
        string synaptic_conductances;
        string potentials;
};

Simulation_consts::Simulation_consts() {
    this->dt = 1e-2;
    this->N = 5;
    this->T = 100000;
    this->I_ext_max = 20.0;
}

Synaptic_consts::Synaptic_consts() {
    this->connection_ratio = 20;
    this->excitatory_ratio = 60;
    this->g_exc = 1e-1;
    this->g_inh = 1e-0;
    this->V_exc = 75.0;
    this->V_inh = 0.0;
}

Data_path::Data_path() {
    this->time = "../data/time.txt";
    this->input_currents = "../data/input_currents.txt";
    this->synaptic_equilibrium_potentials = "../data/synaptic_equilibrium_potentials.txt";
    this->synaptic_conductances = "../data/synaptic_conductances.txt";
    this->potentials = "../data/potentials.txt";
}
#endif
