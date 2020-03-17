# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

from numpy import array, zeros, loadtxt, sin, cos, deg2rad
from re import findall
from .constants import gamma, gamma_m1, R_gas

def P2U(P):
    U = zeros(4)
    U[0] = P[0]
    U[1] = P[0] * P[1]
    U[2] = P[0] * P[2]
    U[3] = P[3] / gamma_m1 + 0.5 * P[0] * (P[1] ** 2 + P[2] ** 2)
    return U

# Hàm P2F: tính dòng qua mặt (công thức ở bài 18)
def P2F(P, side):
    n = side.normal    #vector pháp tuyến đơn vị của mặt
    vn = n.dot(P[1:3]) #vận tốc vuông góc bề mặt V.n
    F = zeros(4)
    F[0] = P[0] * vn
    F[1] = F[0] * P[1] + P[3] * n[0]
    F[2] = F[0] * P[2] + P[3] * n[1]
    F[3] = F[0] * (P[3] / P[0] * gamma / gamma_m1 + 0.5 * (P[1] ** 2 + P[2] ** 2))
    return F * side.area

# Hàm U2P: xác định biến biên thủy P từ biến bảo toàn U
def U2P(U, P):
    P[0] = U[0]
    P[1] = U[1] / U[0]
    P[2] = U[2] / U[0]
    P[3] = (U[3] - 0.5 * P[0] * (P[1] ** 2 + P[2] ** 2)) * gamma_m1


# hàm xác định khối lượng riêng phụ thuộc nhiệt độ và áp suất
# theo phương trình trạng thái
def rho(T, p):
    return p/(R_gas*T)

def Temperature(P):
    return P[3]/(R_gas*P[0])

# vận tốc âm thanh
def VSound(P):
    return (gamma * P[3] / P[0]) ** 0.5

# hàm tính số mach
def Mach(P):
    a = (gamma * P[3] / P[0]) ** 0.5
    u = (P[1] * P[1] + P[2] * P[2]) ** 0.5
    M = u / a
    return M

def Pt2P(M=0.0, Tt=293.15, pt=101.325, alf=0.0):
    m = 1.0/(1.0+M*M*gamma_m1/2.0)
    T = Tt*m
    p = pt*(m**(gamma/gamma_m1))
    r = rho(T, p)
    a = (gamma*p/r)**0.5
    V = M*a
    alf = deg2rad(alf)
    return array([r, V*cos(alf), V*sin(alf), p])

def P2Pt(P):
    r, u, v, p = P[0], P[1], P[2], P[3]
    M = Mach(P)
    m = (1.0+M*M*gamma_m1/2.0)
    T = Temperature(P)
    Tt = T*m
    pt = p*(m**(gamma/gamma_m1))
    rt = r*(m**(1.0/gamma_m1))
    return rt, pt, Tt

def Pm2P(M=0.0, T=293.15, p=101.325, alf=0.0):
    r = rho(T, p)
    a = (gamma*p/r)**0.5
    V = M*a
    alf = deg2rad(alf)
    return array([r, V*cos(alf), V*sin(alf), p])

# Hàm đọc lưới gồm nhiều block
def import_mesh(file_name):
    print('\nImport mesh from: %s\n' % file_name)
    zone_n = 0
    zone_names = []
    zone_id = []
    with open(file_name, 'r') as f:
        for line in f:
            if line[:4] == 'ZONE': # ví dụ: ZONE T="1", I=60, J=40
                zone_n += 1
                # lấy giá trị số thứ tự Block T, số điểm lưới IxJ (NixNj)
                ints = list(map(int, findall(r'\d+', line)))
                zone_names.append(str(ints[0]))
                zone_id.append([ints[1], ints[2]])

    # đọc tọa độ các điểm lưới bằng hàm loadtxt, bỏ 3 hàng đầu
    # dùng reshape để chuyển mảng về 3 chiều
    try: nodes = loadtxt(file_name, usecols=(0, 1), delimiter=' ', skiprows=3, comments=['Z'])
    except: nodes = loadtxt(file_name, usecols=(0, 1), delimiter=',', skiprows=3, comments=['Z'])

    zone_nodes = []
    node_start = 0
    for i in range(zone_n):
        Ni, Nj = zone_id[i][0], zone_id[i][1]
        node_end = node_start + Ni*Nj
        nodes_in = nodes[node_start: node_end].reshape((Nj, Ni, 2))
        node_start = node_end
        zone_nodes.append(nodes_in)

    return zone_n, zone_names, zone_nodes