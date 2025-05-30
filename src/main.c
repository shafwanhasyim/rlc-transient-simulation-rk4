#include <stdio.h>
#include <math.h>

#define L 1.0
#define R 0.0
#define C 0.25
#define E0 1.0
#define W 3.5 // omega
#define P 2.0 // p

// Fungsi sumber tegangan E(t)
double E(double t)
{
    return E0 * sin(W * t);
}

// Sistem ODE: y[0] = q, y[1] = i
void f(double t, double y[], double dydt[])
{
    dydt[0] = y[1];
    dydt[1] = (1.0 / L) * (E(t) - R * y[1] - y[0] / C);
}

// Satu langkah metode RK4
void rk4_step(double t, double y[], double h, double y_out[])
{
    int i;
    double k1[2], k2[2], k3[2], k4[2], yt[2];

    f(t, y, k1);
    for (i = 0; i < 2; i++)
        yt[i] = y[i] + h * k1[i] / 2.0;

    f(t + h / 2.0, yt, k2);
    for (i = 0; i < 2; i++)
        yt[i] = y[i] + h * k2[i] / 2.0;

    f(t + h / 2.0, yt, k3);
    for (i = 0; i < 2; i++)
        yt[i] = y[i] + h * k3[i];

    f(t + h, yt, k4);

    for (i = 0; i < 2; i++)
        y_out[i] = y[i] + (h / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
}

// Solusi analitik q(t) dari Eq. (28.11)
double q_analytic(double t)
{
    double term1 = (-E0 / (L * (P * P - W * W))) * (W / P) * sin(P * t);
    double term2 = (E0 / (L * (P * P - W * W))) * sin(W * t);
    return term1 + term2;
}

int main()
{
    double t0 = 0.0, tf = 20.0, h = 0.01; // gunakan tf lebih kecil supaya hasil jelas
    int N = (int)((tf - t0) / h);
    double t = t0;
    double y[2] = {0.0, 0.0}; // kondisi awal q=0, i=0
    double y_next[2];

    printf("t\tq_RK4\t\tq_Analytic\t\tDifference\n");
    for (int n = 0; n <= N; n++)
    {
        double q_num = y[0];
        double q_an = q_analytic(t);
        double diff = fabs(q_num - q_an);
        // Gunakan notasi ilmiah agar error kecil tetap terlihat
        printf("%.4f\t%.9f\t%.9f\t%.9e\n", t, q_num, q_an, diff);

        rk4_step(t, y, h, y_next);
        t += h;
        y[0] = y_next[0];
        y[1] = y_next[1];
    }

    return 0;
}
