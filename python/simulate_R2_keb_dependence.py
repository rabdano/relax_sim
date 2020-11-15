import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from subprocess import call
import json
import numpy as np
import os
from shutil import rmtree, copy
from scipy.optimize import curve_fit


def inv_rec(t, A, R1):
    return A * (1.0 - 2.0 * np.exp(-t * R1))


def exp_dec(t, A, R2):
    return A * np.exp(-t * R2)


cpmg_taus = [0.003]

for cpmg_tau in cpmg_taus:
    print("\n\n\nCPMG tau =", cpmg_tau, "s")
    wd = f"keb_dependence_tau_{cpmg_tau:.5f}"
    exe = "bm_three_site.exe"
    with open("system.json", 'r') as f:
        pars = json.load(f)

    rmtree(wd, ignore_errors=True)
    os.makedirs(wd)


    #
    #  Simulate dependence of R2 on number of exchangeable sites
    #
    kebs = np.logspace(2, 5, 25)
    R1 = np.zeros_like(kebs, dtype="double")
    R2 = np.zeros_like(kebs, dtype="double")

    c = 0.5e-3
    tau_e = 0.71e-3  # Knispel, 1975, Chem. Phys. Lett, 32, 238
    # cmpg_tau = 0.05
    rho = 997.0479
    Mw = 18.01528
    nb = 384
    ne = 130

    for i, keb in enumerate(kebs):
        print("keb =", keb)
        pars["pf"] = (2.0 * rho - 2.0 * nb * c * Mw) / (2 * rho + ne * c * Mw)
        pars["pb"] = (2.0 * nb * c * Mw) / (2 * rho + ne * c * Mw)
        pars["pe"] = (ne * c * Mw) / (2 * rho + ne * c * Mw)
        pars["kbf"] = 1.0 / tau_e * pars["pf"] * np.log(2) / (pars["pf"] + pars["pb"])
        pars["keb"] = keb
        pars["CPMG_tau"] = [cpmg_tau]

        with open(wd + "/system.json", "w") as f:
            f.write(json.dumps(pars, sort_keys=True, indent=4))

        copy(exe, wd + "/" + exe)

        call([exe], cwd=wd)

        spin_lattice = np.genfromtxt(wd+ "/inv_rec.txt", delimiter=",")
        spin_spin = np.genfromtxt(wd + f"/cpmg_tau_{cpmg_tau:.6f}.txt", delimiter=",")


        # fit R1 data
        delays = spin_lattice[:, 0]
        mz = spin_lattice[:, 1] + spin_lattice[:, 2]

        popt_R1, pcov_R1 = curve_fit(f=inv_rec,
                                     xdata=delays,
                                     ydata=mz,
                                     p0=(1.0, 0.5),
                                     bounds=(0, [1.0, 1000]))

        Az_fit = popt_R1[0]
        R1_fit = popt_R1[1]
        R1[i] = popt_R1[1]
        print("R1 =", R1[i])

        delays_fit = np.linspace(delays[0], delays[-1], 1000)
        mz_fit = inv_rec(delays_fit, Az_fit, R1_fit)

        plt.figure(figsize=(16.0/2.54, 12.0/2.54))
        plt.plot(delays, mz, marker="o", linestyle="None", color="b", ms=5, label=r"R$\rm_1$ data")
        plt.plot(delays_fit, mz_fit, color="r", label=r"R$\rm_1$ fit")
        plt.xlabel("Time, s")
        plt.ylabel(r"M$\rm_z$")
        plt.legend()
        plt.savefig(wd + f"/R1_{keb:.2f}hz.png", bbox_inches="tight", dpi=300)
        plt.close()


        # fit R2 data
        delays = spin_spin[:, 0]
        mxy = spin_spin[:, 1] + spin_spin[:, 2]

        popt_R2, pcov_R2 = curve_fit(f=exp_dec,
                                     xdata=delays,
                                     ydata=mxy,
                                     p0=(1.0, 0.5),
                                     bounds=(0, [1.0, 1000]))

        Axy_fit = popt_R2[0]
        R2_fit = popt_R2[1]
        R2[i] = R2_fit
        print("R2 =", R2[i])

        delays_fit = np.linspace(delays[0], delays[-1], 1000)
        mxy_fit = exp_dec(delays_fit, Axy_fit, R2_fit)

        plt.figure(figsize=(16.0 / 2.54, 12.0 / 2.54))
        plt.plot(delays, mxy, marker="o", linestyle="None", color="b", ms=5, label=r"R$\rm_2$ data")
        plt.plot(delays_fit, mxy_fit, color="r", label=r"R$\rm_2$ fit")
        plt.xlabel("Time, s")
        plt.ylabel(r"M$\rm_{xy}$")
        plt.legend()
        plt.savefig(wd + f"/R2_{keb:.2f}hz.png", bbox_inches="tight", dpi=300)
        plt.close()


    plt.figure(figsize=(16.0 / 2.54, 12.0 / 2.54))
    plt.plot(kebs, R1, marker="o", color="b", ms=5, label=r"R$\rm_1$")
    plt.plot(kebs, R2, marker="o", color="r", ms=5, label=r"R$\rm_2$")
    plt.xlabel(r"k$\rm_{eb}$, 1/s")
    plt.ylabel(r"R$\rm_{1,2}$, s$^{\rm-1}$")
    plt.legend()
    plt.xscale('log')
    plt.savefig(wd + "/keb_dependence.png", bbox_inches="tight", dpi=300)
    plt.close()

    np.savetxt(wd + "/kebs_R1_R2.txt", np.vstack((kebs, R1, R2)).T, fmt="%.5f", delimiter=", ")
