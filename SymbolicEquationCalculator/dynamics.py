from sympy import Matrix, symbols, cos, sin, sec, Function, diff, simplify
from sympy.physics.mechanics import dynamicsymbols

def sq(p):
    return p.dot(p.T)

def shorten(str):
    return str.replace("Derivative", "d").replace("sin", "s").replace("cos", "c").replace("sec", "sc").replace("**","^")

def MProduct(M, i, j, k):
    return M[convert(i) - 1] * M[convert(j) - 1] * M[convert(k) - 1] 

def convert(i):
    return int(((( i - i % 10 ) / 10) - 1) * 4) + i % 10

m1, m2, m3, m4, I = symbols("m1 m2 m3 m4 I")
θ1, θ2, θ3, Φ, ψ = dynamicsymbols("θ1 θ2 θ3 Φ ψ")
l1, l2, l3, l4 = symbols("l1 l2 l3 l4")
λ1, λ2, λ3, λ4 = symbols("λ1 λ2 λ3 λ4")
g, t = symbols("g t")


# p1 = λ1 * (Matrix(3, 1, [cos(θ1), 0, sin(θ1)]))

# p2 = (
#     l1 * (Matrix(3, 1, [cos(θ1), 0, sin(θ1)])) + 
#     λ2 * (Matrix(3, 1, [cos(θ2), 0, sin(θ2)])))

# p3 = (
#     l1 * (Matrix(3, 1, [cos(θ1), 0, sin(θ1)])) + 
#     l2 * (Matrix(3, 1, [cos(θ2), 0, sin(θ2)])) + 
#     λ3 * (Matrix(3, 1, [cos(θ3), 0, sin(θ3)])))

# p4 = (
#     l1 * (Matrix(3, 1, [cos(θ1), 0, sin(θ1)])) + 
#     l2 * (Matrix(3, 1, [cos(θ2), 0, sin(θ2)])) + 
#     l3 * (Matrix(3, 1, [cos(θ2), 0, sin(θ2)])) + 
#     λ4 * l4 * sec(Φ) * Matrix(3, 1, [cos(θ3), 0, sin(θ3)]) + 
#     λ4 * (Matrix(3, 1, [0 , cos(Φ), 0])))

# # print(p1)
# # print(p2)
# # print(p3)
# # print(p4)

# K1 = 0.5 * m1 * sq(diff(p1, t))
# K2 = 0.5 * m2 * sq(diff(p2, t))
# K3 = 0.5 * m3 * sq(diff(p3, t))
# K4 = 0.5 * m4 * sq(diff(p4, t)) + 0.5 * I * diff(ψ) * diff(ψ)

# U1 = m1 * g * p1[2]
# U2 = m2 * g * p2[2]
# U3 = m3 * g * p3[2]
# U4 = m4 * g * p4[2]

# L = K1 + K2 + K3 + K4 - U1 - U2 - U3 - U4

# # ω1 = diff(θ1, t)
# # ω2 = diff(θ2, t)
# # ω3 = diff(θ3, t)
# # ω4 = diff(Φ, t)

# # T1 = diff(L, ω1)
# # diff(L, θ1)

# T1 = diff(diff(L, diff(θ1, t)), t) - diff(L, θ1)
# T2 = diff(diff(L, diff(θ2, t)), t) - diff(L, θ2)
# T3 = diff(diff(L, diff(θ3, t)), t) - diff(L, θ3)
# T4 = diff(diff(L, diff(Φ, t)), t) - diff(L, Φ)

# t1 = shorten(str(T1))
# t2 = shorten(str(T2))
# t3 = shorten(str(T3))
# t4 = shorten(str(T4))

# print("T1: " + t1)
# print("\n")
# print(T2: " + t2)
# print("\n")
# print(T3: " + t3)
# print("\n")
# print(T4: " + t4)
# print("\n")


# CHAT GPT GENERATED MAX MATRIX RECOMPUTE BY SENDING INTO GPT PREVIOUS T1-T4 EQUATIONS
cos12 = cos(θ1 - θ2)
cos13 = cos(θ1 - θ3)
cos_phi_term = cos(Φ - θ1 + θ3)
m234 = m2 + m3 + m4

# Matrix entries
M11 = l1**2 * m234 + m1 * λ1**2 + l1 * l2 * m234 * cos12 + l1 * l3 * m3 * cos13 + l1 * l4 * m4 * cos_phi_term
M12 = l1**2 * m234 + l1 * l2 * m234 * cos12
M13 = l1 * l3 * m3 * cos13 + l1 * l4 * m4 * cos_phi_term
M14 = l1 * l4 * m4 * cos_phi_term

M22 = l2**2 * m234 + m2 * λ2**2 + l1 * l2 * m234 * cos12
M33 = l3**2 * m3 + l4**2 * m4 + m3 * λ3**2 + m4 * λ4**2 + l1 * l3 * m3 * cos13 + l1 * l4 * m4 * cos_phi_term
M34 = l1 * l4 * m4 * cos_phi_term
M44 = l4**2 * m4 + m4 * λ4**2 + l1 * l4 * m4 * cos_phi_term

M = Matrix([
    [M11, M12, M13, M14],
    [M12, M22, 0,    0],
    [M13, 0,   M33,  M34],
    [M14, 0,   M34,  M44]
])
# CHAT GPT GENERATED MASS MATRIX



q = Matrix([
    θ1, θ2, θ3, Φ
])

dq = Matrix([
    [diff(θ1, t), diff(θ2, t), diff(θ3, t), diff(Φ, t)]
])

C = Matrix([
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0]
])
n = 4

for i in range(4):
    for j in range(4):
        C[i * 4 + j] = sum(0.5 * (diff(M[i * 4 + j], q[k]) + diff(M[i * 4 + k], q[j]) - diff(M[k * 4 + j], q[i])) * dq[k] for k in range(n))



# CHAT GPT GENERATED GRAVITY MATRIX RECOMPUTE BY SENDING INTO GPT PREVIOUS t1-t4 EQUATIONS
# Reused terms
m234 = m2 + m3 + m4

# Gravity vector components
G1 = g * l1 * m234 * cos(θ1) + g * m1 * λ1 * cos(θ1)
G2 = g * l2 * m234 * cos(θ2) + g * m2 * λ2 * cos(θ2)
G3 = (g * l3 * m3 * cos(θ3) + g * l4 * m4 * cos(θ3 + Φ)
      + g * m3 * λ3 * cos(θ3) + g * m4 * λ4 * cos(θ3 + Φ))
G4 = g * l4 * m4 * cos(θ3 + Φ) + g * m4 * λ4 * cos(θ3 + Φ)

# Combine into gravity vector
G = Matrix([G1, G2, G3, G4])

print("M: " + shorten(str(M)) + "\n")

print("C: " + shorten(str(C)) + "\n")

print("G: " + shorten(str(G)))
# i = 13
# print(int((i % 10 + (( i - i % 10 ) / 10) * 4)))

# for i in range(1,5):
#     for j in range(1,5):
#         print(convert(i * 10 + j))

# mdet = (
#     M[0] * (                       
#         MProduct(M, 22, 33, 44) + MProduct(M, 23, 34, 42) + MProduct(M, 24, 32, 43) - MProduct(M, 24, 33, 42) - MProduct(M, 23, 32, 44) - MProduct(M, 22, 34, 43))
#     - M[4] * (                     
#         MProduct(M, 12, 33, 44) + MProduct(M, 13, 34, 42) + MProduct(M, 14, 32, 43) - MProduct(M, 14, 33, 42) - MProduct(M, 13, 32, 44) - MProduct(M, 12, 34, 43))
#     - M[8] * (                     
#         MProduct(M, 12, 23, 44) + MProduct(M, 13, 24, 42) + MProduct(M, 14, 22, 43) - MProduct(M, 14, 23, 42) - MProduct(M, 13, 22, 44) - MProduct(M, 12, 24, 43))
#     - M[12] * (                    
#         MProduct(M, 12, 23, 34) + MProduct(M, 13, 24, 32) + MProduct(M, 14, 22, 33) - MProduct(M, 14, 23, 32) - MProduct(M, 13, 22, 34) - MProduct(M, 12, 24, 33))
# )

# m11 = (MProduct(M, 22, 33, 44) + MProduct(M, 23, 34, 42) + MProduct(M, 24, 32, 43)
#     - MProduct(M, 24, 33, 42) - MProduct(M, 23, 32, 44) - MProduct(M, 22, 34, 43))

# m12 = (MProduct(M, 21, 33, 44) + MProduct(M, 23, 34, 41) + MProduct(M, 24, 31, 43)
#     - MProduct(M, 24, 33, 41) - MProduct(M, 23, 31, 44) - MProduct(M, 21, 34, 43))

# m13 = (MProduct(M, 21, 32, 44) + MProduct(M, 22, 34, 41) + MProduct(M, 24, 31, 42)
#     - MProduct(M, 24, 32, 41) - MProduct(M, 22, 31, 44) - MProduct(M, 21, 34, 42))

# m14 = (MProduct(M, 21, 32, 43) + MProduct(M, 22, 33, 41) + MProduct(M, 23, 31, 42)
#     - MProduct(M, 23, 32, 41) - MProduct(M, 22, 31, 43) - MProduct(M, 21, 33, 42))

# m21 = (MProduct(M, 12, 33, 44) + MProduct(M, 13, 34, 42) + MProduct(M, 14, 32, 43)
#     - MProduct(M, 14, 33, 42) - MProduct(M, 13, 32, 44) - MProduct(M, 12, 34, 43))

# m22 = (MProduct(M, 11, 33, 44) + MProduct(M, 13, 34, 41) + MProduct(M, 14, 31, 43)
#     - MProduct(M, 14, 33, 41) - MProduct(M, 13, 31, 44) - MProduct(M, 11, 34, 43))

# m23 = (MProduct(M, 11, 32, 44) + MProduct(M, 12, 34, 41) + MProduct(M, 14, 31, 42)
#     - MProduct(M, 14, 32, 41) - MProduct(M, 12, 31, 44) - MProduct(M, 11, 34, 42))

# m24 = (MProduct(M, 11, 32, 43) + MProduct(M, 12, 33, 41) + MProduct(M, 13, 31, 42)
#     - MProduct(M, 13, 32, 41) - MProduct(M, 12, 31, 43) - MProduct(M, 11, 33, 42))

# m31 = (MProduct(M, 12, 23, 44) + MProduct(M, 13, 24, 42) + MProduct(M, 14, 22, 43)
#     - MProduct(M, 14, 23, 42) - MProduct(M, 13, 22, 44) - MProduct(M, 12, 24, 43))

# m32 = (MProduct(M, 11, 23, 44) + MProduct(M, 13, 24, 41) + MProduct(M, 14, 21, 43)
#     - MProduct(M, 14, 23, 41) - MProduct(M, 13, 21, 44) - MProduct(M, 11, 24, 43))

# m33 = (MProduct(M, 11, 22, 44) + MProduct(M, 12, 24, 41) + MProduct(M, 14, 21, 42)
#     - MProduct(M, 14, 22, 41) - MProduct(M, 12, 21, 44) - MProduct(M, 11, 24, 42))

# m34 = (MProduct(M, 11, 22, 43) + MProduct(M, 12, 23, 41) + MProduct(M, 13, 21, 42)
#     - MProduct(M, 13, 22, 41) - MProduct(M, 12, 21, 43) - MProduct(M, 11, 23, 42))

# m41 = (MProduct(M, 12, 23, 34) + MProduct(M, 13, 24, 32) + MProduct(M, 14, 22, 33)
#     - MProduct(M, 14, 23, 32) - MProduct(M, 13, 22, 34) - MProduct(M, 12, 24, 33))

# m42 = (MProduct(M, 11, 23, 34) + MProduct(M, 13, 24, 31) + MProduct(M, 14, 21, 33)
#     - MProduct(M, 14, 23, 31) - MProduct(M, 13, 21, 34) - MProduct(M, 11, 24, 33))

# m43 = (MProduct(M, 11, 22, 34) + MProduct(M, 12, 24, 31) + MProduct(M, 14, 21, 32)
#     - MProduct(M, 14, 22, 31) - MProduct(M, 12, 21, 34) - MProduct(M, 11, 24, 32))

# m44 = (MProduct(M, 11, 22, 33) + MProduct(M, 12, 23, 31) + MProduct(M, 13, 21, 32)
#     - MProduct(M, 13, 22, 31) - MProduct(M, 12, 21, 33) - MProduct(M, 11, 23, 32))

# mInv = Matrix([
#     [m11, m12, m13, m14],
#     [m21, m22, m23, m24],
#     [m31, m32, m33, m34],
#     [m41, m42, m43, m44]
# ])

# mInv = (1/mdet) * mInv

# actualDq = dq = Matrix([
#     [diff(θ1, t)],
#     [diff(θ2, t)],
#     [diff(θ3, t)],
#     [diff(θ3, t)]
# ])


# E = simplify(mInv * (-G - (C* actualDq)))

# dE = E.jacobian(q)

# print(C * actualDq)


# m12 = MProduct(M, 21, 33, 44) + MProduct(M, 23, 34, 41) + MProduct(M, 24, 31, 43) - MProduct(M, 24, 33, 41) - MProduct(M, 23, 31, 44) - MProduct(M, 21, 34, 43)
# m13 = MProduct(M, 21, 32, 44) + MProduct(M, 22, 34, 41) + MProduct(M, 24, 32, 41) - MProduct(M, 24, 32, 41) - MProduct(M, 23, 31, 44) - MProduct(M, 21, 34, 43)

# print(mInv)

# print(mdet)

# print(E)

# print("dE1: " + E.jacobian(q))

# print("dE2: " + E.jacobian(dq))
