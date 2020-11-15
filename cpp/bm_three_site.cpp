#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include "nlohmann/json.hpp"
using json = nlohmann::json;
using namespace std;


double dIxf_dt(double Ixf, double Iyf, double Ixb, double R2f, double wf, double kfb, double kbf)
{
    return (-R2f - kfb) * Ixf - wf * Iyf + kbf * Ixb;
}

double dIyf_dt(double Ixf, double Iyf, double Iyb, double R2f, double wf, double kfb, double kbf)
{
    return wf * Ixf + (-R2f - kfb) * Iyf + kbf * Iyb;
}

double dIzf_dt(double Izf, double Izb, double R1f, double Ieqf, double kfb, double kbf)
{
    return R1f * Ieqf + (-R1f - kfb) * Izf + kbf * Izb;
}


double dIxb_dt(double Ixf, double Ixb, double Iyb, double Ixe, double R2b, double wb, double kfb, double kbf, double kbe, double keb)
{
    return kfb * Ixf + (-R2b - kbf - kbe) * Ixb - wb * Iyb + keb * Ixe;
}

double dIyb_dt(double Iyf, double Ixb, double Iyb, double Iye, double R2b, double wb, double kfb, double kbf, double kbe, double keb)
{
    return kfb * Iyf + wb * Ixb + (-R2b - kbf - kbe) * Iyb + keb * Iye;
}

double dIzb_dt(double Izf, double Izb, double Ize, double R1b, double Ieqb, double kfb, double kbf, double kbe, double keb)
{
    return R1b * Ieqb + kfb * Izf + (-R1b - kbf - kbe) * Izb + keb * Ize;
}


double dIxe_dt(double Ixb, double Ixe, double Iye, double R2e, double we, double kbe, double keb)
{
    return kbe * Ixb + (-R2e - keb) * Ixe - we * Iye;
}

double dIye_dt(double Iyb, double Ixe, double Iye, double R2e, double we, double kbe, double keb)
{
    return kbe * Iyb + we * Ixe + (-R2e - keb) * Iye;
}

double dIze_dt(double Izb, double Ize, double R1e, double Ieqe, double kbe, double keb)
{
    return R1e * Ieqe + kbe * Izb + (-R1e - keb) * Ize;
}

void evolution(vector <double> &M, json system, double tend, bool sample)
{
    // system parameters
    double R1f = system["R1f"];
    double R2f = system["R2f"];
    double pf = system["pf"];
    double wf = 2 * M_PI * system["of_hz"];

    double R1b = system["R1b"];
    double R2b = system["R2b"];
    double pb = system["pb"];
    double wb = 2 * M_PI * system["ob_hz"];

    double R1e = system["R1e"];
    double R2e = system["R2e"];
    double pe = system["pe"];
    double we = 2 * M_PI * system["oe_hz"];

    double kbf = system["kbf"];
    double kfb = (pb / pf) * kbf;

    double keb = system["keb"];
    double kbe = (pe / pb) * keb;

    // time grid
    double dt = system["dt"];
    double n = tend / dt;
    

    // FID sampling
    double dw = system["DW"];
    double stride = dw / dt;


    // initialization
    double time = 0;

    double Ixf = 0; double Ixf_prev = 0; double Ixf1; double Ixf2; double Ixf3; double Ixf4;
    double Iyf = 0; double Iyf_prev = 0; double Iyf1; double Iyf2; double Iyf3; double Iyf4;
    double Izf = 0; double Izf_prev = 0; double Izf1; double Izf2; double Izf3; double Izf4;

    double Ixb = 0; double Ixb_prev = 0; double Ixb1; double Ixb2; double Ixb3; double Ixb4;
    double Iyb = 0; double Iyb_prev = 0; double Iyb1; double Iyb2; double Iyb3; double Iyb4;
    double Izb = 0; double Izb_prev = 0; double Izb1; double Izb2; double Izb3; double Izb4;

    double Ixe = 0; double Ixe_prev = 0; double Ixe1; double Ixe2; double Ixe3; double Ixe4;
    double Iye = 0; double Iye_prev = 0; double Iye1; double Iye2; double Iye3; double Iye4;
    double Ize = 0; double Ize_prev = 0; double Ize1; double Ize2; double Ize3; double Ize4;


    // initial conditions (after 90-deg pulse with (-y) phase -- all magnetization is along X)
    double Ieqf = pf;
    double Ieqb = pb;
    double Ieqe = pe;

    // Ixf = pf;
    // Ixb = pb;
    // Ixe = pe;

    Ixf = M[0];
    Iyf = M[1];
    Izf = M[2];
    Ixb = M[3];
    Iyb = M[4];
    Izb = M[5];
    Ixe = M[6];
    Iye = M[7];
    Ize = M[8];

    if (sample == true) {
        cout << time << " ";
        cout << Ixf << " " << Iyf << " " << Izf << " ";
        cout << Ixb << " " << Iyb << " " << Izb << " ";
        cout << Ixe << " " << Iye << " " << Ize << " ";
        cout << endl;
    }

    
    for (int i=1; i <= n; i++) {
        // time
        time = dt * i;


        // set previous step magnetization
        Ixf_prev = Ixf;
        Iyf_prev = Iyf;
        Izf_prev = Izf;
        Ixb_prev = Ixb;
        Iyb_prev = Iyb;
        Izb_prev = Izb;
        Ixe_prev = Ixe;
        Iye_prev = Iye;
        Ize_prev = Ize;


        // Y1
        Ixf1 = dt * dIxf_dt(Ixf_prev, Iyf_prev, Ixb_prev, R2f, wf, kfb, kbf);
        Iyf1 = dt * dIyf_dt(Ixf_prev, Iyf_prev, Iyb_prev, R2f, wf, kfb, kbf);
        Izf1 = dt * dIzf_dt(Izf_prev, Izb_prev, R1f, Ieqf, kfb, kbf);

        Ixb1 = dt * dIxb_dt(Ixf_prev, Ixb_prev, Iyb_prev, Ixe_prev, R2b, wb, kfb, kbf, kbe, keb);
        Iyb1 = dt * dIyb_dt(Iyf_prev, Ixb_prev, Iyb_prev, Iye_prev, R2b, wb, kfb, kbf, kbe, keb);
        Izb1 = dt * dIzb_dt(Izf_prev, Izb_prev, Ize_prev, R1b, Ieqb, kfb, kbf, kbe, keb);

        Ixe1 = dt * dIxe_dt(Ixb_prev, Ixe_prev, Iye_prev, R2e, we, kbe, keb);
        Iye1 = dt * dIye_dt(Iyb_prev, Ixe_prev, Iye_prev, R2e, we, kbe, keb);
        Ize1 = dt * dIze_dt(Izb_prev, Ize_prev, R1e, Ieqe, kbe, keb);


        // Y2
        Ixf2 = dt * dIxf_dt(Ixf_prev + Ixf1 / 2, Iyf_prev + Iyf1 / 2, Ixb_prev + Ixb1 / 2, R2f, wf, kfb, kbf);
        Iyf2 = dt * dIyf_dt(Ixf_prev + Ixf1 / 2, Iyf_prev + Iyf1 / 2, Iyb_prev + Iyb1 / 2, R2f, wf, kfb, kbf);
        Izf2 = dt * dIzf_dt(Izf_prev + Izf1 / 2, Izb_prev + Izb1 / 2, R1f, Ieqf, kfb, kbf);

        Ixb2 = dt * dIxb_dt(Ixf_prev + Ixf1 / 2, Ixb_prev + Ixb1 / 2, Iyb_prev + Iyb1 / 2, Ixe_prev + Ixe1 / 2, R2b, wb, kfb, kbf, kbe, keb);
        Iyb2 = dt * dIyb_dt(Iyf_prev + Iyf1 / 2, Ixb_prev + Ixb1 / 2, Iyb_prev + Iyb1 / 2, Iye_prev + Iye1 / 2, R2b, wb, kfb, kbf, kbe, keb);
        Izb2 = dt * dIzb_dt(Izf_prev + Izf1 / 2, Izb_prev + Izb1 / 2, Ize_prev + Ize1 / 2, R1b, Ieqb, kfb, kbf, kbe, keb);

        Ixe2 = dt * dIxe_dt(Ixb_prev + Ixb1 / 2, Ixe_prev + Ixe1 / 2, Iye_prev + Iye1 / 2, R2e, we, kbe, keb);
        Iye2 = dt * dIye_dt(Iyb_prev + Iyb1 / 2, Ixe_prev + Ixe1 / 2, Iye_prev + Iye1 / 2, R2e, we, kbe, keb);
        Ize2 = dt * dIze_dt(Izb_prev + Izb1 / 2, Ize_prev + Ize1 / 2, R1e, Ieqe, kbe, keb);


        // Y3
        Ixf3 = dt * dIxf_dt(Ixf_prev + Ixf2 / 2, Iyf_prev + Iyf2 / 2, Ixb_prev + Ixb2 / 2, R2f, wf, kfb, kbf);
        Iyf3 = dt * dIyf_dt(Ixf_prev + Ixf2 / 2, Iyf_prev + Iyf2 / 2, Iyb_prev + Iyb2 / 2, R2f, wf, kfb, kbf);
        Izf3 = dt * dIzf_dt(Izf_prev + Izf2 / 2, Izb_prev + Izb2 / 2, R1f, Ieqf, kfb, kbf);

        Ixb3 = dt * dIxb_dt(Ixf_prev + Ixf2 / 2, Ixb_prev + Ixb2 / 2, Iyb_prev + Iyb2 / 2, Ixe_prev + Ixe2 / 2, R2b, wb, kfb, kbf, kbe, keb);
        Iyb3 = dt * dIyb_dt(Iyf_prev + Iyf2 / 2, Ixb_prev + Ixb2 / 2, Iyb_prev + Iyb2 / 2, Iye_prev + Iye2 / 2, R2b, wb, kfb, kbf, kbe, keb);
        Izb3 = dt * dIzb_dt(Izf_prev + Izf2 / 2, Izb_prev + Izb2 / 2, Ize_prev + Ize2 / 2, R1b, Ieqb, kfb, kbf, kbe, keb);

        Ixe3 = dt * dIxe_dt(Ixb_prev + Ixb2 / 2, Ixe_prev + Ixe2 / 2, Iye_prev + Iye2 / 2, R2e, we, kbe, keb);
        Iye3 = dt * dIye_dt(Iyb_prev + Iyb2 / 2, Ixe_prev + Ixe2 / 2, Iye_prev + Iye2 / 2, R2e, we, kbe, keb);
        Ize3 = dt * dIze_dt(Izb_prev + Izb2 / 2, Ize_prev + Ize2 / 2, R1e, Ieqe, kbe, keb);


        // Y4
        Ixf4 = dt * dIxf_dt(Ixf_prev + Ixf3, Iyf_prev + Iyf3, Ixb_prev + Ixb3, R2f, wf, kfb, kbf);
        Iyf4 = dt * dIyf_dt(Ixf_prev + Ixf3, Iyf_prev + Iyf3, Iyb_prev + Iyb3, R2f, wf, kfb, kbf);
        Izf4 = dt * dIzf_dt(Izf_prev + Izf3, Izb_prev + Izb3, R1f, Ieqf, kfb, kbf);

        Ixb4 = dt * dIxb_dt(Ixf_prev + Ixf3, Ixb_prev + Ixb3, Iyb_prev + Iyb3, Ixe_prev + Ixe3, R2b, wb, kfb, kbf, kbe, keb);
        Iyb4 = dt * dIyb_dt(Iyf_prev + Iyf3, Ixb_prev + Ixb3, Iyb_prev + Iyb3, Iye_prev + Iye3, R2b, wb, kfb, kbf, kbe, keb);
        Izb4 = dt * dIzb_dt(Izf_prev + Izf3, Izb_prev + Izb3, Ize_prev + Ize3, R1b, Ieqb, kfb, kbf, kbe, keb);

        Ixe4 = dt * dIxe_dt(Ixb_prev + Ixb3, Ixe_prev + Ixe3, Iye_prev + Iye3, R2e, we, kbe, keb);
        Iye4 = dt * dIye_dt(Iyb_prev + Iyb3, Ixe_prev + Ixe3, Iye_prev + Iye3, R2e, we, kbe, keb);
        Ize4 = dt * dIze_dt(Izb_prev + Izb3, Ize_prev + Ize3, R1e, Ieqe, kbe, keb);


        // Y
        Ixf = Ixf_prev + (Ixf1 + 2 * Ixf2 + 2 * Ixf3 + Ixf4) / 6;
        Iyf = Iyf_prev + (Iyf1 + 2 * Iyf2 + 2 * Iyf3 + Iyf4) / 6;
        Izf = Izf_prev + (Izf1 + 2 * Izf2 + 2 * Izf3 + Izf4) / 6;

        Ixb = Ixb_prev + (Ixb1 + 2 * Ixb2 + 2 * Ixb3 + Ixb4) / 6;
        Iyb = Iyb_prev + (Iyb1 + 2 * Iyb2 + 2 * Iyb3 + Iyb4) / 6;
        Izb = Izb_prev + (Izb1 + 2 * Izb2 + 2 * Izb3 + Izb4) / 6;

        Ixe = Ixe_prev + (Ixe1 + 2 * Ixe2 + 2 * Ixe3 + Ixe4) / 6;
        Iye = Iye_prev + (Iye1 + 2 * Iye2 + 2 * Iye3 + Iye4) / 6;
        Ize = Ize_prev + (Ize1 + 2 * Ize2 + 2 * Ize3 + Ize4) / 6;


        if (((i % (int)stride) == 0) && (sample == true)) {
            cout << time << " ";
            cout << Ixf << " " << Iyf << " " << Izf << " ";
            cout << Ixb << " " << Iyb << " " << Izb << " ";
            cout << Ixe << " " << Iye << " " << Ize << " ";
            cout << endl;
        }
    }

    M[0] = Ixf;
    M[1] = Iyf;
    M[2] = Izf;
    M[3] = Ixb;
    M[4] = Iyb;
    M[5] = Izb;
    M[6] = Ixe;
    M[7] = Iye;
    M[8] = Ize;
}

void pulse_90(vector <double> &M, int phase)
{
    vector <double> M_prev;
    M_prev = M;
    if (phase == 0) {
        M[0] = M_prev[0]; M[1] = -M_prev[2]; M[2] = M_prev[1];
        M[3] = M_prev[3]; M[4] = -M_prev[5]; M[5] = M_prev[4];
        M[6] = M_prev[6]; M[7] = -M_prev[8]; M[8] = M_prev[7];
    } else if (phase == 1) {
        M[0] = M_prev[2]; M[1] = M_prev[1]; M[2] = -M_prev[0];
        M[3] = M_prev[5]; M[4] = M_prev[4]; M[5] = -M_prev[3];
        M[6] = M_prev[8]; M[7] = M_prev[7]; M[8] = -M_prev[6];
    } else if (phase == 2) {
        M[0] = M_prev[0]; M[1] = M_prev[2]; M[2] = -M_prev[1];
        M[3] = M_prev[3]; M[4] = M_prev[5]; M[5] = -M_prev[4];
        M[6] = M_prev[6]; M[7] = M_prev[8]; M[8] = -M_prev[7];
    } else if (phase == 3) {
        M[0] = -M_prev[2]; M[1] = M_prev[1]; M[2] = M_prev[0];
        M[3] = -M_prev[5]; M[4] = M_prev[4]; M[5] = M_prev[3];
        M[6] = -M_prev[8]; M[7] = M_prev[7]; M[8] = M_prev[6];
    }
}

void pulse_180(vector <double> &M, int phase)
{
    if ((phase == 0) || (phase == 2)) {
        M[0] = M[0]; M[1] = -M[1]; M[2] = -M[2];
        M[3] = M[3]; M[4] = -M[4]; M[5] = -M[5];
        M[6] = M[6]; M[7] = -M[7]; M[8] = -M[8];
    } else {
        M[0] = -M[0]; M[1] = M[1]; M[2] = -M[2];
        M[3] = -M[3]; M[4] = M[4]; M[5] = -M[5];
        M[6] = -M[6]; M[7] = M[7]; M[8] = -M[8];        
    }

}

void read_intensity(vector <double> M)
{
    for (int i = 0; i < M.size(); ++i) {
        if (i < (M.size() - 1)) {
            printf("%11.6f, ", M[i]);   
        } else {
            printf("%11.6f", M[i]);
        }
    }
    printf("\n");
}

void equlibrate(vector <double> &M)
{
    ifstream i("system.json");
    json system;
    i >> system;

    M[0] = 0;
    M[1] = 0;
    M[2] = (double)system["pf"];
    M[3] = 0;
    M[4] = 0;
    M[5] = (double)system["pb"];
    M[6] = 0;
    M[7] = 0;
    M[8] = (double)system["pe"];
}

void inv_rec(vector <double> M, json system, vector <double> delay)
{
    equlibrate(M);

    FILE * f = fopen("inv_rec.txt", "w");
    cout << "Simulating Inversion Recovery..." << endl;

    for (int i = 0; i < delay.size(); ++i) {
        pulse_180(M, 0);
        evolution(M, system, delay[i], false);
        
        fprintf(f, "%11.6f, %11.6f, %11.6f, %11.6f\n", delay[i], M[2], M[5], M[8]);

        equlibrate(M);
    }

    fclose(f);
}

void cpmg(vector <double> M, json system, vector <double> delay, double tau)
{
    double time = 0;
    int cur_delay = 0;
    char buf[100];
    double Ixyf, Ixyb, Ixye;

    equlibrate(M);
    
    pulse_90(M, 0);
    
    sprintf(buf, "cpmg_tau_%0.6f.txt", tau);
    FILE * f = fopen(buf, "w");
    cout << "Simulating CPMG for tau = " << tau << "..." << endl;

    while (time <= delay.back()) {
        if (time >= delay[cur_delay]) {
            Ixyf = sqrt(pow(M[0], 2) + pow(M[1], 2));
            Ixyb = sqrt(pow(M[3], 2) + pow(M[4], 2));
            Ixye = sqrt(pow(M[6], 2) + pow(M[7], 2));
            fprintf(f, "%11.6f, %11.6f, %11.6f, %11.6f\n", time, Ixyf, Ixyb, Ixye);
            cur_delay++;
        }

        evolution(M, system, tau / 2, false);
        time += tau / 2;

        pulse_180(M, 0);

        evolution(M, system, tau / 2, false);
        time += tau / 2;
    }
    fclose(f);
}

int main()
{
    ifstream i("system.json");
    json system;
    i >> system;

    vector <double> M;

    M.push_back(0);
    M.push_back(0);
    M.push_back((double)system["pf"]);
    M.push_back(0);
    M.push_back(0);
    M.push_back((double)system["pb"]);
    M.push_back(0);
    M.push_back(0);
    M.push_back((double)system["pe"]);

    vector <double> inv_rec_delay = system["inv_rec_delays"];
    vector <double> cpmg_delay = system["CPMG_delays"];
    vector <double> tau = system["CPMG_tau"];
    inv_rec(M, system, inv_rec_delay);
    for (int i = 0; i < tau.size(); ++i) {
        cpmg(M, system, cpmg_delay, tau[i]);
    }

    return 0;
}
