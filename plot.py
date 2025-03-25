import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def read_log(filename):
    data = np.loadtxt(filename, skiprows=1, dtype=np.float64)
    t = data[:, 0]
    pos = data[:, 1:4]
    return t, pos

t_expl, expl = read_log('explicit.log')
t_impl, impl = read_log('implicit.log') 
t_adms, adms = read_log('adams.log')
t_adapt, adapt = read_log('adaptive.log')

plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(adapt[:, 0], adapt[:, 2], label='Adaptive', lw=2)
plt.plot(expl[:, 0], expl[:, 2], '--', label='Explicit')
plt.plot(impl[:, 0], impl[:, 2], '--', label='Implicit')
plt.plot(adms[:, 0], adms[:, 2], '--', label='Adams')
plt.title("Trajectory (XZ plane)")
plt.xlabel("X")
plt.ylabel("Z")
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
for name, t, pos in [('Explicit', t_expl, expl), 
                    ('Implicit', t_impl, impl),
                    ('Adams', t_adms, adms)]:
    interp_x = interp1d(t, pos[:, 0], fill_value="extrapolate")
    interp_y = interp1d(t, pos[:, 1], fill_value="extrapolate")
    
    x_interp = interp_x(t_adapt)
    y_interp = interp_y(t_adapt)
    
    error = np.sqrt((x_interp - adapt[:, 0])**2 + (y_interp - adapt[:, 1])**2)
    plt.plot(t_adapt, error, label=name)

plt.yscale('log')
plt.title("Position error (log scale)")
plt.xlabel("Time")
plt.ylabel("Error")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('trajectory_comparison.png', dpi=300)
