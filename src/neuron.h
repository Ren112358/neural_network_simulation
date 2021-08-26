// ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1392413/pdf/jphysiol01442-0106.pdf
#include <iostream>
#include <cmath>

#include "consts.h"

using namespace std;

Simulation_consts sim_consts;
const double dt = sim_consts.dt;  // sampling interval

class Neuron {
public:
    Neuron();
    Neuron(
            double V,
            double m,
            double h,
            double n,
            double C,
            double g_Na,
            double g_K,
            double g_leak,
            double V_Na,
            double V_K,
            double V_leak
            );
    void update_V(double V, double m, double h, double n, double I_ext);  // for single neuron
    void update_V(double V, double m, double h, double n, double I_ext, double I_syn);  // for multiple neurons
    void update_m(double V, double m);
    void update_h(double V, double h);
    void update_n(double V, double n);
    void update_state(double I_ext);  // for single neuron
    void update_state(double I_ext, double I_syn);  // for multiple neurons
    double get_V();  // to record potential

private:
    double V;  // membrane potential
    double m;  // channel variable
    double h;  // channel variable
    double n;  // channel variable
    double C;  // membrane capacitance
    double V_Na;  // Sodium equilibrium potential
    double V_K;  // Potasium equilibrium potential
    double V_leak;  // leak current equilibrium potential
    double g_Na;  // Sodium conductance
    double g_K;  // Potasium conductance
    double g_leak;  // leak current conductance

    /* channel get open rate functions */
    double alpha_m(double V);
    double beta_m(double V);

    double alpha_h(double V);
    double beta_h(double V);

    double alpha_n(double V);
    double beta_n(double V);
};

Neuron::Neuron() {
    this->V = 0;
    this->m = 0;
    this->h = 0;
    this->n = 0;
    this->C = 1.0;
    this->V_Na = 115.0;
    this->V_K = -12.0;
    this->V_leak = 10.613;
    this->g_Na = 120.0;
    this->g_K = 36.0;
    this->g_leak = 0.3;
}

Neuron::Neuron(
        double V,
        double m,
        double h,
        double n,
        double C,
        double V_Na,
        double V_K,
        double V_leak,
        double g_Na,
        double g_K,
        double g_leak
        ) {
    this->V = V;
    this->m = m;
    this->h = h;
    this->n = n;
    this->C = C;
    this->V_Na = V_Na;
    this->V_K = V_K;
    this->V_leak = V_leak;
    this->g_Na = g_Na;
    this->g_K = g_K;
    this->g_leak = g_leak;
}

double Neuron::alpha_m(double V) {
    double numerator = 0.1 * (-V + 25);
    double denominator = exp((-V + 25) / 10) - 1;

    if (denominator < 1e-10) {
        return 1.0;
    }
    return numerator / denominator;
}

double Neuron::beta_m(double V) {
    return 4 * exp(-V / 18);
}

double Neuron::alpha_h(double V) {
    return 0.07 * exp(-V / 20);
}

double Neuron::beta_h(double V) {
    double numerator = 1.0;
    double denominator = exp((-V + 30) / 10) + 1;
    return numerator / denominator;
}

double Neuron::alpha_n(double V) {
    double numerator = 0.01 * (-V + 10);
    double denominator = exp((-V + 10) / 10) - 1;
    if (denominator < 1e-10) {
        return 0.1;
    }
    return numerator / denominator;
}

double Neuron::beta_n(double V) {
    return 0.125 * exp(-V / 80);
}

void Neuron::update_V(double V, double m, double h, double n, double I_ext) {
    double I_Na = this->g_Na * pow(m, 3) * h * (V - V_Na);
    double I_K = this->g_K * pow(n, 4) * (V - V_K);
    double I_leak = this->g_leak * (V - V_leak);
    double I_ion = I_Na + I_K + I_leak;
    double I = I_ext - I_ion;
    this->V = V + dt * (I / C);
}

void Neuron::update_V(double V, double m, double h, double n, double I_ext, double I_syn) {
    double I_Na = this->g_Na * pow(m, 3) * h * (V - V_Na);
    double I_K = this->g_K * pow(n, 4) * (V - V_K);
    double I_leak = this->g_leak * (V - V_leak);
    double I_ion = I_Na + I_K + I_leak;
    double I = I_ext - I_ion - I_syn;
    this->V = V + dt * (I / C);
}

void Neuron::update_m(double V, double m) {
    this->m = m + dt * ((1 - m) * alpha_m(V) - m * beta_m(V));
}

void Neuron::update_h(double V, double h) {
    this->h = h + dt * ((1 - h) * alpha_h(V) - h * beta_h(V));
}

void Neuron::update_n(double V, double n) {
    this->n = n + dt * ((1 - n) * alpha_n(V) - n * beta_n(V));
}

void Neuron::update_state(double I_ext) {
    double V = this->V;
    double m = this->m;
    double h = this->h;
    double n = this->n;
    update_V(V, m, h, n, I_ext);
    update_m(V, m);
    update_h(V, h);
    update_n(V, n);
}

void Neuron::update_state(double I_ext, double I_syn) {
    double V = this->V;
    double m = this->m;
    double h = this->h;
    double n = this->n;
    update_V(V, m, h, n, I_ext, I_syn);
    update_m(V, m);
    update_h(V, h);
    update_n(V, n);
}

double Neuron::get_V() {
    return this->V;
}
