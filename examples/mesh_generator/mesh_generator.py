# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

import  numpy as np

# Hàm xuất lưới
def export_mesh(nodes_list, file_name):
    f = open(file_name, 'w')
    f.write('TITLE = "vncfd python"\n')
    f.write('VARIABLES = "X", "Y"\n')

    for n, nodes in enumerate(nodes_list):
        Nj, Ni = nodes.shape[0], nodes.shape[1]
        f.write('ZONE T="%d", I= %d, J= %d\n' % (n+1, Ni, Nj))
        for j in range(Nj):
            for i in range(Ni):
                f.write('%f %f\n' % (nodes[j, i, 0], nodes[j, i, 1]))
    f.close()

def generate_mesh_1D(Nj=2, Ni=501):
    # Kích thước vùng tính toán
    ly, lx = 0.01, 1.0
    # Tạo mảng 3 chiều để chứa tọa độ các điểm lưới
    nodes = np.zeros((Nj, Ni, 2))
    # tọa độ x tại các điểm lưới
    x = np.linspace(0, lx, Ni)
    dy = ly
    # xác định tọa độ (x, y) của từng điểm
    for j in range(Nj):
        for i in range(Ni):
            nodes[j, i, 0] = x[i]
            nodes[j, i, 1] = j*dy
    export_mesh([nodes], '../test1D/test1D.mesh')
    return

generate_mesh_1D()

# Tạo lưới bài toán Mach 2
def generate_mesh_M2(Nj=41, Ni=101):
    # Kích thước vùng tính toán
    ly, lx = 4.0, 10.0
    # Tạo mảng 3 chiều để chứa tọa độ các điểm lưới 
    nodes = np.zeros((Nj, Ni, 2))
    # tọa độ x tại các điểm lưới
    dx = lx / Ni
    x = np.linspace(0, lx, Ni)
    # tọa độ y của biên dưới
    y0 = np.zeros(Ni)
    # index i tương ứng vị trí x = 2, 4 trên biên dưới
    i2 = int(2./dx)
    i4 = int(4./dx)

    y0[i2:i4] = (x[i2:i4]-2.)*np.tan(np.pi/12)
    y0[i4:] = 2.0*np.tan(np.pi/12)

    # khoảng cách dy giữa hai điểm cùng cột
    dy = np.array([(ly-y)/(Nj-1) for y in y0])
    
    # xác định tọa độ (x, y) của từng điểm 
    for j in range(Nj):
        for i in range(Ni):
            nodes[j, i, 0] = x[i]
            nodes[j, i, 1] = y0[i]+j*dy[i]

    export_mesh([nodes], 'data/mach_2_1blk.dat')
    export_mesh([nodes[:, :Ni // 2 + 1, :], nodes[:, Ni // 2:, :]], 'data/mach_2_2blk.dat')
    return

def generate_mesh_wedge(Nj=51, Ni=101):
    # Kích thước vùng tính toán
    ly, lx = 10.0, 10.0
    # Tạo mảng 3 chiều để chứa tọa độ các điểm lưới
    nodes = np.zeros((Nj, Ni, 2))
    # tọa độ x tại các điểm lưới
    dx = lx / Ni
    x = np.linspace(0, lx, Ni)
    # tọa độ y của biên dưới
    y0 = np.zeros(Ni)
    # index i tương ứng vị trí x = 2, 4 trên biên dưới
    i3 = int(3. / dx)
    y0[i3:] = (x[i3:]-3.0) * np.tan(np.pi / 12)
    lyi = np.array([(ly - y) for y in y0])

    ratio = 1.1
    dy0i = np.array([li/(ratio**(Nj-1)-1.0)*(ratio-1.0) for li in lyi])
    # xác định tọa độ (x, y) của từng điểm
    for j in range(Nj):
        for i in range(Ni):
            nodes[j, i, 0] = x[i]
            if j == 0: nodes[j, i, 1] = y0[i]
            else: nodes[j, i, 1] = nodes[j-1, i, 1] + dy0i[i]*(ratio**(j-1))

    # export_mesh([nodes], '../wedge/mach_2_1blk.dat')
    export_mesh([nodes[:, :i3+20 + 1, :], nodes[:, i3+20:, :]], '../wedge/wedge.mesh')
    return

# generate_mesh_wedge()

def generate_mesh_backstep(dx=0.1, dy=0.1):
    x1 = np.arange(-11.0, 0.0+dx, dx)
    y1 = np.arange(1.0, 9.0+dy, dy)
    X1, Y1 = np.meshgrid(x1, y1)
    nodes1 = np.zeros((X1.shape[0], X1.shape[1], 2))
    nodes1[:, :, 0] = X1
    nodes1[:, :, 1] = Y1
    x2 = np.arange(0.0, 10.0 + dx, dx)
    y2 = np.arange(0.0, 9.0 + dy, dy)
    X2, Y2 = np.meshgrid(x2, y2)
    nodes2 = np.zeros((X2.shape[0], X2.shape[1], 2))
    nodes2[:, :, 0] = X2
    nodes2[:, :, 1] = Y2
    export_mesh([nodes1, nodes2], 'data/backstep_2blk.dat')
    return

# generate_mesh_backstep(dx=0.1, dy=0.1)

# Tạo lưới bài toán dòng chảy quanh hình trụ
def generate_mesh_cylinder(Nj=41, Ni=61):
    # chia góc 2Pi ra thành Ni điểm lưới
    alpha = np.linspace(0.0, -2 * np.pi, Ni)

    # bán kính tại các điểm lưới
    r = np.zeros(Nj)
    r[0] = 1.0  # bán kính hình trụ
    dr = 1e-2  # kích thước ('dộ dày') ô lưới đầu tiên sát biên (bề mặt trụ)
    ratio = 1.2  # tỷ lệ tăng kích thước ô lưới
    for j in range(1, Nj):
        r[j] = r[j - 1] + dr
        dr *= ratio

        # tạo độ điểm lưới
    nodes = np.zeros((Nj, Ni, 2))
    for j in range(Nj):
        for i in range(Ni):
            nodes[j, i, 0] = r[j] * np.cos(alpha[i])
            nodes[j, i, 1] = r[j] * np.sin(alpha[i])
    export_mesh([nodes], 'data/cylinder_mesh.dat')
    return nodes
