import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parametry wahadła
l = 0.055  # Długość wahadła (m)
g = 9.81  # Przyspieszenie grawitacyjne (m/s^2)
A = 0.010  # Amplituda oscylacji punktu zawieszenia (m)
m = 1 #masa wahadła
f = 20 #częstotliwość
omega = 2*np.pi*f  # prędkość kołowa (rad/s)
omega_theta0 = 0.0  # Początkowa prędkość kątowa (rad/s)
theta0_degrees = 180  # Początkowy kąt wychylenia (stopnie)
t_span = (0, 100)  # Czas od 0 do 100 sekun

# Równania ruchu dla wahadła Kapitzy
def equations(t, y):
    theta, omega_theta = y
    dtheta_dt = omega_theta
    domega_theta_dt = -(g/l) * np.sin(theta) - (A * omega**2 / l) * np.sin(theta) * np.cos(omega * t)
    return [dtheta_dt, domega_theta_dt]


theta0=np.deg2rad(theta0_degrees)
y0 = [theta0, omega_theta0]

# Zakres czasu
t_eval = np.linspace(t_span[0], t_span[1], 1000)

# Rozwiązanie równań różniczkowych
solution = solve_ivp(equations, t_span, y0, t_eval=t_eval)

# Konwersja kąta z radianów na stopnie
theta_degrees = np.degrees(solution.y[0])

# Wykres kąta od czasu
plt.figure(figsize=(10, 6))
plt.ylim(-200, 200)
plt.plot(solution.t, theta_degrees, label='Kąt θ(t)')
plt.xlabel('Czas (s)')
plt.ylabel('Kąt (stopnie)')
plt.title('Ruch wahadła Kapitzy')
plt.legend()
plt.grid(True)
plt.show()

# Potencjał efektywny
def effective_potential(theta):
    return -m * g * l * np.cos(theta) - (m * A**2 * omega**2 / 2) * np.cos(theta)**2

# Zakres kątów
theta = np.linspace(-1.2*np.pi, 1.2*np.pi, 1000)

# Obliczanie potencjału efektywnego dla każdego kąta
V_eff = effective_potential(theta)

# Konwersja kąta na stopnie
theta_deg = np.rad2deg(theta)

# Wykres
plt.figure(figsize=(10, 6))
plt.plot(theta_deg, V_eff, label=r'$V_{\text{eff}}(\theta)$')
plt.xlim(-200,200)
plt.xlabel('Kąt [stopnie]')
plt.ylabel('Potencjał efektywny [J]')
plt.title('Potencjał efektywny wahadła Kapitzy')
plt.legend()
plt.grid(True)
plt.show()

