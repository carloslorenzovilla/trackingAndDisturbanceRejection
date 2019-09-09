import matplotlib.pyplot as plt
import numpy as np
from control import *
from control import pzmap
from control.matlab import lsim
from mpldatacursor import datacursor


#Lab 7 Variables

K_b = 0.07 # V/rad/s = N-m/A
K_t = K_b
J = 6E-5 # kg-m^2
R = 2. # ohms
L = 2E-3 # Henries
B = 0.0004 # N-m/rad/sec
T_L = 0.# No Load Yet

tt = np.arange(0., 10.001, 1e-3)

#1

beta = K_t/R
delta = K_t/((B*R)+(K_b*K_t))
a = K_t/(J*R) #v(t), L=0
b = ((B*R) + (K_b*K_t))/(J*R) #denominator, L=0

#Proportional Controller
K_p = 1
num_c = K_p
den_c = 1
Gc_s = tf(num_c, den_c)

#Velocity Controller
RVC = series(Gc_s, tf([a], [1, b]))

#Position Controller
RPC = series(RVC, tf([1], [1, 0]))
M_s2 = feedback(RPC, 1)

U = np.ones((len(tt), 2))
idx = np.where(tt < 1)[0]
U[idx, 0] *= 0
U[:, 0] *= 1
U[:, 1] *= -T_L/beta

num_R = [1]
num_D = [den_c/num_c]
cPoly = [1]

num = [[num_R, num_D]]
den = [[cPoly, cPoly]]

theta_s2 = series(tf(num, den), M_s2)

y, t, x = lsim(theta_s2, U, tt, 0)

plt.figure()
line1 = plt.plot(t, y, markersize = 4, label=r'$\tau_L$ = 0.0')
datacursor(line1, draggable=True)

#2

T_L = 0.002

U2 = np.ones((len(tt), 2))
idx = np.where(tt < 1)[0]
U2[idx, 0] *= 0
U2[:, 0] *= 1
U2[:, 1] *= -T_L/beta

y2, t2, x2 = lsim(theta_s2, U2, tt, 0)

line1 = plt.plot(t2, y2, markersize = 4, label=r'$\tau_L$ = 0.002')
datacursor(line1, draggable=True)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.title(r'$\theta(t)$ Step Response via Transfer Function', fontsize=18)
plt.ylabel(r'$\theta(t)$ (rad)', fontsize=14)
plt.xlabel('Time (sec)', fontsize=14)
plt.grid()
plt.text(9, 1.1, 'Carlos Villa', fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
plt.legend(loc='best', fontsize=14)


#3 PI
T_L = 0.002
alpha = 1.31
K = 0.75

#Proportional Integrator
num_Gc_s = [K, K*alpha]
den_Gc_s = [1, 0]
Gc_s = tf(num_Gc_s, den_Gc_s)

#Rotational Velocity Control
RVC = series(Gc_s, tf([a], [1, b]))
print('RVC:\n', RVC)

#Rotational Position Control
RPC = series(RVC, tf([1], [1, 0]))
print('RPC:\n', RPC)

#Closed Loop
M_s3 = feedback(RPC, 1)
print('M_s3:\n', M_s3)

#MISO System
Ref = tf(1, 1)
print('Ref:\n', Ref)
Dist = tf(den_Gc_s, num_Gc_s)
print('Dist:\n', Dist)
Ref_plus_Dist = Ref - (Dist*(T_L/beta))
print('Ref_plus_Dist:\n', Ref_plus_Dist)

#Position
theta_s3 = series(Ref_plus_Dist, M_s3)
print('theta_s3:\n', theta_s3)

U3 = np.ones(len(tt))
idx3 = np.where(tt < 1)[0]
U3[idx3] *= 0
U3[:] *= 1

yS3, tS3, xS3 = lsim(theta_s3, U3, tt, 0)

plt.figure()
line1 = plt.plot(t, y, markersize = 4, label=r'Reference')
datacursor(line1, draggable=True)
line2 = plt.plot(tS3, yS3, markersize = 4, label=r'$\alpha$ = 1.31, K = 0.75')
datacursor(line2, draggable=True)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.title(r'PI $\theta(t)$ Step Response via Transfer Function', fontsize=18)
plt.ylabel(r'$\theta(t)$ (rad)', fontsize=14)
plt.xlabel('Time (sec)', fontsize=14)
plt.grid()
plt.text(9, 1.15, 'Carlos Villa', fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
plt.legend(loc='best', fontsize=14)

#4 PID

X = 3.1
zeta = 0.39
w_n = 5.2
the_nums = np.convolve([1, X], [1, 2*zeta*w_n, w_n**2])

T_L = 0.002

Ki = the_nums[3]/a
Kp = the_nums[2]/a
Kt = (the_nums[1]-b)/a

poles = np.roots([1, (b+(Kt*a)), a*Kp, a*Ki])

#Proportional Integrator
num_Gc_s = [Kp, Ki]
den_Gc_s = [1, 0]
Gc_s = tf(num_Gc_s, den_Gc_s)

#Rotational Velocity Control with tachometer (differentiator)
RVC = series(Gc_s, feedback(tf([a], [1, b]), Kt))
print('RVC:\n', RVC)

#Rotational Position Control
RPC = series(RVC, tf([1], [1, 0]))
print('RPC:\n', RPC)

#Closed Loop TF
M_s4 = feedback(RPC, 1)
print('M_s4:\n', M_s4)

#MISO System
Ref = tf(1, 1)
print('Ref:\n', Ref)
Dist = tf(den_Gc_s, num_Gc_s)
print('Dist:\n', Dist)
Ref_plus_Dist = Ref - (Dist*(T_L/beta))
print('Ref_plus_Dist:\n', Ref_plus_Dist)

#Position
theta_s4 = series(Ref_plus_Dist, M_s4)
print('theta_s4:\n', theta_s4)

U4 = np.ones(len(tt))
idx4 = np.where(tt < 1)[0]
U4[idx4] *= 0
U4[:] *= 1

yS, tS, xS = lsim(theta_s4, U4, tt, 0)

poles4 = np.roots([1, (b+(Kt*a)), a*Kp, a*Ki])

plt.figure()
line1 = plt.plot(t, y, markersize = 4, label=r'Reference')
datacursor(line1, draggable=True)
line2 = plt.plot(tS, yS, markersize = 4, label='Kp: %1.3f\nKi: %1.3f\nKt: %1.3f' % (Kp, Ki, Kt))
datacursor(line2, display = 'multiple', draggable=True)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.title(r'PID $\theta(t)$ Step Response via Transfer Function', fontsize=18)
plt.ylabel(r'$\theta(t)$ (rad)', fontsize=14)
plt.xlabel('Time (sec)', fontsize=14)
plt.grid()
plt.text(9, 1.15, 'Carlos Villa', fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
plt.legend(loc='best', fontsize=14)
plt.text(2, 0.3, poles4, fontsize=12, bbox=dict(facecolor='white', alpha=0.5))


print('Kp: %1.3f' % Kp)
print('Ki: %1.3f' % Ki)
print('Kt: %1.3f' % Kt)

#5 1% Settle Time 1s < t_s < 2s

X = 4.2
zeta = 0.7
w_n = 5.8
the_nums2 = np.convolve([1, X], [1, 2*zeta*w_n, w_n**2])

T_L = 0.002

Ki = the_nums2[3]/a
Kp = the_nums2[2]/a
Kt = (the_nums2[1]-b)/a

#Proportional Integrator
num_Gc_s = [Kp, Ki]
den_Gc_s = [1, 0]
Gc_s = tf(num_Gc_s, den_Gc_s)

#Rotational Velocity Control with tachometer (differentiator)
RVC = series(Gc_s, feedback(tf([a], [1, b]), Kt))
print('RVC:\n', RVC)

#Rotational Position Control
RPC = series(RVC, tf([1], [1, 0]))
print('RPC:\n', RPC)

#Closed Loop TF
M_s5 = feedback(RPC, 1)
print('M_s5:\n', M_s5)

#MISO System
Ref = tf(1, 1)
print('Ref:\n', Ref)
Dist = tf(den_Gc_s, num_Gc_s)
print('Dist:\n', Dist)
Ref_plus_Dist = Ref - (Dist*(T_L/beta))
print('Ref_plus_Dist:\n', Ref_plus_Dist)

#Position
theta_s5 = series(Ref_plus_Dist, M_s5)
print('theta_s5:\n', theta_s5)

U5 = np.ones(len(tt))
idx5 = np.where(tt < 1)[0]
U5[idx5] *= 0
U5[:] *= 1

yS, tS, xS = lsim(theta_s5, U5, tt, 0)

poles5 = np.roots([1, (b+(Kt*a)), a*Kp, a*Ki])

plt.figure()
line1 = plt.plot(t, y, markersize = 4, label=r'Reference')
datacursor(line1, draggable=True)
line2 = plt.plot(tS, yS, markersize = 4, label='Kp: %1.3f\nKi: %1.3f\nKt: %1.3f' % (Kp, Ki, Kt))
datacursor(line2, draggable=True)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.title(r'PID $\theta(t)$ Step Response via Transfer Function', fontsize=18)
plt.ylabel(r'$\theta(t)$ (rad)', fontsize=14)
plt.xlabel('Time (sec)', fontsize=14)
plt.text(9, 1.15, 'Carlos Villa', fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
plt.legend(loc='best', fontsize=14)
plt.text(2, 0.3, poles5, fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

plt.grid()

print('Kp: %1.3f' % Kp)
print('Ki: %1.3f' % Ki)
print('Kt: %1.3f' % Kt)

#6

#Proportional Controller
K_p = 1

Gc_s = tf(K_p, K_p)

#Velocity Controller
RVC = series(Gc_s, tf([a], [1, b]))

#Position Controller
RPC = series(RVC, tf([1], [1, 0]))
M_s6 = feedback(RPC, 1)

input6 = tf([[[1.], [K_p]]], [[[1.], [K_p]]])

theta_s6 = series(input6, M_s6)

tt6 = np.arange(0, 15, 1e-3)
Wn = 35 - 3 #Letter C = 3
T_L = 0.02

U6 = np.ones((len(tt6), 2))
idx = np.where(tt < 2)[0]
U6[idx, 0] *= 0
U6[:, 0] *= 1
U6[:, 1] *= -T_L/beta*np.sin(Wn*tt6)

y6, t6, x6 = lsim(theta_s6, U6, tt6, 0)

plt.figure()
plt.plot(tt6, U6[:, 0], markersize = 4, label=r'Reference')
plt.plot(t6, y6, markersize = 4, label=r'$\tau_L$ = $0.02sin(32t)$')
plt.tick_params(axis='both', which='major', labelsize=12)
plt.title(r'$\theta(t)$ Step Response via Transfer Function', fontsize=18)
plt.ylabel(r'$\theta(t)$ (rad)', fontsize=14)
plt.xlabel('Time (sec)', fontsize=14)
plt.grid()
plt.text(13.75, 1.35, 'Carlos Villa', fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
plt.legend(loc='best', fontsize=14)

#7

tt7 = np.arange(0, 15, 1e-3)
wn = 35 - 3 #Letter C = 3
T_L = 0.02

U7 = np.ones((len(tt7), 2))
idx = np.where(tt < 2)[0]
U7[idx, 0] *= 0
U7[:, 0] *= 1
U7[:, 1] *= -T_L/beta*np.sin(wn*tt7)

alpha = -1.31
K = -64.5
num7 = [a*K, a*K*alpha]
denm7 = [1, b, wn**2, ((wn**2)*b)+(a*K), a*K*alpha]

poles7 = np.roots(denm7)
M_s7 = tf(num7, denm7)

input7 = tf([[[1.], [1., 0, wn**2]]], [[[1.], [K, K*alpha]]])
theta_s7 = series(input7, M_s7)

y7, t7, x7 = lsim(theta_s7, U7, tt7, 0)

plt.figure()
plt.plot(tt7, U7[:, 0], markersize = 4, label=r'Reference')
plt.plot(t7, y7, markersize = 4, label='K: %1.1f' '\n' r'$\alpha$: %1.2f' % (K, alpha))
plt.tick_params(axis='both', which='major', labelsize=12)
plt.title(r'$\theta(t)$ Step Response via Transfer Function', fontsize=18)
plt.ylabel(r'$\theta(t)$ (rad)', fontsize=14)
plt.xlabel('Time (sec)', fontsize=14)
plt.text(13.75, 1.35, 'Carlos Villa', fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
plt.text(3, -2, poles7, fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
plt.legend(loc='best', fontsize=14)
plt.grid()

print('K: %1.1f' % K)
print('alpha: %1.2f' % alpha)

plt.show()

# -(((b**2)*alpha) + ((Wn**2)*b))/a
# -(((Wn**2)*b)+(a*K))/(b**2)
