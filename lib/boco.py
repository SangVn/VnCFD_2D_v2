# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

from .constants import gamma, gamma_m1
from .functions import Mach, VSound, P2F
from .fluxes import flux_roe_fortran

P_freestream = [293.15, 0.0, 0.0, 101325] # thông số dòng tự do
p_exit = 101325.0 # áp suất tại mặt outlet

def set_boco_const(Pfree=None, pexit=None):
    global  P_freestream, p_exit
    P_freestream = Pfree
    p_exit = pexit

# Hàm sign_ic: xác định `thêm hay bớt` dòng ở ô lưới kề bên side
def sign_ic(ic): #ic ~ index_cell; left_cell: ic = 0; right_cell: ic = 1
    if ic == 0: return -1.0 # res -= flux
    else: return 1.0        # res += flux

# Điều kiện dòng chảy vào: P_out = P_freestream
def supersonic_inflow(boundary, ic):
    for side in boundary:
        # F = P2F(P_freestream, side) # có thể tính trực tiếp F
        P_in = side.cells[ic].P
        # xác định P_left và P_right để tính flux
        if ic == 1: F = flux_roe_fortran(side, P_freestream, P_in)
        else: F = flux_roe_fortran(side, P_in, P_freestream)
        side.cells[ic].res += sign_ic(ic) * F

# Điều kiện dòng chảy ra: P_out = P_in
def supersonic_outflow(boundary, ic):
    for side in boundary:
        F = P2F(side.cells[ic].P, side)
        side.cells[ic].res += sign_ic(ic) * F

# Điều kiện biên no_slip: P_side: u=v=0, dp/dn=0, drho/dn=0
def no_slip(boundary, ic):
    for side in boundary:
        P_in = side.cells[ic].P
        P_out = [P_in[0], 0.0, 0.0, P_in[3]] # P trên biên
        # F = P2F(P_out, side)
        # P_out = [P_in[0], -P_in[1], -P_in[2], P_in[3]] # P ở ghost cell
        if ic == 1: F = flux_roe_fortran(side, P_out, P_in)
        else: F = flux_roe_fortran(side, P_in, P_out)
        side.cells[ic].res += sign_ic(ic)*F

# symmetry and slip boco are the same: V_out = V_in - 2(V_in.n)n
def symmetry(boundary, ic):
    for side in boundary:
        P_in = side.cells[ic].P
        V_in = P_in[1:3]
        Vn = side.normal.dot(V_in)
        # V_out = V_in - 2*Vn*side.normal # V ở ghost cell
        V_out = V_in -  Vn * side.normal  # V trên biên
        P_out = [P_in[0], V_out[0], V_out[1], P_in[3]]
        # F = P2F(P_out, side)
        if ic == 1: F = flux_roe_fortran(side, P_out, P_in)
        else: F = flux_roe_fortran(side, P_in, P_out)
        side.cells[ic].res += sign_ic(ic) * F

# Điều kiện biên farfield: P_out = P_freestream
def farfield(boundary, ic):
    for side in boundary:
        # F = P2F(P_freestream, side)
        P_in = side.cells[ic].P
        if ic == 1: F = flux_roe_fortran(side, P_freestream, P_in)
        else: F = flux_roe_fortran(side, P_in, P_freestream)
        side.cells[ic].res += sign_ic(ic) * F

# Điều kiện biên outflow xác định theo các đường đặc trưng
def outflow(boundary, ic):
    for side in boundary:
        P_in = side.cells[ic].P
        rho_b = P_in[0]*(p_exit/P_in[3])**(1.0/gamma)
        a_b = (gamma*p_exit/rho_b)**0.5
        Vn_in = side.normal.dot(P_in[1:3])
        a_in = VSound(P_in)
        R_p = Vn_in + 2*a_in/(gamma_m1)
        Vn_b = R_p - 2*a_b/(gamma_m1)
        V_b = P_in[1:3] + (Vn_b - Vn_in)*side.normal
        P_b = [rho_b, V_b[0], V_b[1], p_exit]
        # F = P2F(P_b, side)
        if ic == 1: F = flux_roe_fortran(side, P_b, P_in)
        else: F = flux_roe_fortran(side, P_in, P_b)
        side.cells[ic].res += sign_ic(ic) * F

# Điều kiện biên inlow xác định theo các đường đặc trưng
def inflow(boundary, ic):
    P_e = P_freestream
    for side in boundary:
        P_in = side.cells[ic].P
        if (Mach(P_in) >= 1.0):
            P_b = P_e
        else:
            Vn_in = side.normal.dot(P_in[1:3])
            Vn_e  = side.normal.dot(P_e[1:3])
            a_in = VSound(P_in)
            a_e  = VSound(P_e)
            R_p = Vn_e + 2 * a_e / (gamma_m1)
            R_m = Vn_in - 2 * a_in / (gamma_m1)
            Vn_b = 0.5*(R_p+R_m)
            a_b = 0.25*(gamma_m1)*(R_p - R_m)
            V_b = P_e[1:3] + (Vn_b - Vn_e)*side.normal
            R = P_e[3]/(P_e[0]**gamma)
            rho_b = (a_b*a_b/(gamma*R))**(1.0/gamma_m1)
            p_b = R*rho_b**gamma
            P_b = [rho_b, V_b[0], V_b[1], p_b]

        # F = P2F(P_b, side)
        if ic == 1: F = flux_roe_fortran(side, P_b, P_in)
        else: F = flux_roe_fortran(side, P_in, P_b)
        side.cells[ic].res += sign_ic(ic)*F

# Điều kiện biên joint
def joint(boundary, ic):
    pass

# Điều kiện biên null
def null(boundary, ic):
    pass
