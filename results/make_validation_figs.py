"""검증 그림: (1) 반암시적 vs 명시적 (R=1 kOhm, ~11 ns)  (2) v17 수렴SOR vs v18 vs v17 원본(60회)"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, json
plt.rcParams["font.family"] = ["Liberation Serif", "DejaVu Serif"]; plt.rcParams.update({"font.size": 11})
# ---- (1)
S = np.load("results/semi_11ns.npz"); E = np.load("results/expl_11ns.npz")
rs, re_ = S["rows"], E["rows"]
fig, ax = plt.subplots(2, 2, figsize=(12, 8))
for r, lb, c, ls in ((re_, f"explicit (fixed LU), {int(E['steps'])} substeps, {float(E['wall']):.0f} s", "k", "-"),
                     (rs, f"adaptive semi-implicit, {int(S['steps'])} substeps, {float(S['wall']):.0f} s", "tab:red", "--")):
    t = r[:, 0]*1e9
    ax[0, 0].plot(t, r[:, 2]/1e3, color=c, ls=ls, lw=1.8, label=lb)
    ax[0, 1].plot(t, r[:, 3], color=c, ls=ls, lw=1.8, label=lb)
    ax[1, 0].semilogy(t, r[:, 4], color=c, ls=ls, lw=1.8, label=lb)
    ax[1, 1].semilogy(t, r[:, 1]*1e12, color=c, ls=ls, lw=1.8, label=lb)
ax[0, 0].set(xlabel="t [ns]", ylabel="Needle voltage [kV]", title="Needle voltage (R = 1 kΩ, Vp = 30 kV)")
ax[0, 1].set(xlabel="t [ns]", ylabel="Discharge current I [A]", title="Discharge current (Sato)")
ax[1, 0].set(xlabel="t [ns]", ylabel="max $n_e$ [m$^{-3}$]", title="Peak electron density")
ax[1, 1].set(xlabel="t [ns]", ylabel="dt [ps]", title="Time step")
for a in ax.flat: a.grid(alpha=0.3); a.legend(fontsize=9)
fig.suptitle("Semi-implicit Poisson (adaptive) vs explicit scheme - needle-plane 10 mm, nx=161", fontsize=13)
fig.tight_layout(); fig.savefig("results/fig_semi_vs_explicit.png", dpi=120); plt.close(fig)
# ---- (2)
A = np.load("results/verify_converged.npy")
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
t = A[:, 0]
ax[0].semilogy(t, A[:, 3], "k-", lw=2.2, label="v17 + converged SOR (3000 it/substep)")
ax[0].semilogy(t, A[:, 4], "r--", lw=1.8, label="v18 (sparse direct)")
ax[0].semilogy(t, A[:, 5], "b:", lw=1.8, label="v17 original (60 it/substep)")
ax[0].set(xlabel="t [ns]", ylabel="max $n_e$ [m$^{-3}$]", title="Peak electron density")
ax[1].plot(t, A[:, 6], "k-", lw=2.2, label="v17 + converged SOR"); ax[1].plot(t, A[:, 7], "r--", lw=1.8, label="v18")
ax[1].plot(t, A[:, 8], "b:", lw=1.8, label="v17 original (60 it)")
ax[1].set(xlabel="t [ns]", ylabel="max |E| [kV/cm]", title="Peak field")
ax[2].semilogy(t, np.maximum(A[:, 9], 1e-6), "r-", lw=1.8)
ax[2].set(xlabel="t [ns]", ylabel="max |V(v17 conv.) - V(v18)| [V]", title="Potential difference (of 30 kV)")
for a in ax: a.grid(alpha=0.3)
ax[0].legend(fontsize=9); ax[1].legend(fontsize=9)
fig.suptitle("v18 reproduces v17 with a converged Poisson solver (identical dt sequence, 0-3 ns)", fontsize=13)
fig.tight_layout(); fig.savefig("results/fig_v17conv_vs_v18.png", dpi=120); plt.close(fig)
print("saved results/fig_semi_vs_explicit.png, results/fig_v17conv_vs_v18.png")
# 요약 수치
i = -1
print(f"semi vs explicit @11 ns: Va {rs[i,2]:.0f} vs {re_[i,2]:.0f} V ; I {rs[i,3]:.2f} vs {re_[i,3]:.2f} A ; ne {rs[i,4]:.2e} vs {re_[i,4]:.2e}")
# 브리징 시각
tb_s = rs[np.nonzero(rs[:,6] <= 1e-9)[0][0], 0]*1e9 if (rs[:,6]<=1e-9).any() else None
tb_e = re_[np.nonzero(re_[:,6] <= 1e-9)[0][0], 0]*1e9 if (re_[:,6]<=1e-9).any() else None
print("head reaches plate: semi", tb_s, "explicit", tb_e)
print("wall speedup semi/explicit to 11 ns:", float(E["wall"])/float(S["wall"]))
