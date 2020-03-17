# coding: utf-8
# Nguyên mẫu FORTRAN - Katate Masatsuka
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

from numpy import zeros
from .constants import gamma, gamma_m1
from .fluxes_fortran import roe

def flux_roe_fortran(side, PL, PR):
    njk = side.normal
    flux, wsn = roe(PL, PR, njk)
    return flux * side.area

def flux_roe_python(side, PL, PR):
    nx = side.normal[0]
    ny = side.normal[1]
    tx = -ny
    ty = nx

    # Left state
    rhoL, vxL, vyL, pL = PL[0], PL[1], PL[2], PL[3]
    vnL = vxL * nx + vyL * ny
    vtL = vxL * tx + vyL * ty
    aL = (gamma * pL / rhoL) ** 0.5
    HL = pL / rhoL * gamma / gamma_m1 + 0.5 * (vxL ** 2 + vyL ** 2)

    # Right state
    rhoR, vxR, vyR, pR = PR[0], PR[1], PR[2], PR[3]
    vnR = vxR * nx + vyR * ny
    vtR = vxR * tx + vyR * ty
    aR = (gamma * pR / rhoR) ** 0.5
    HR = pR / rhoR * gamma / gamma_m1 + 0.5 * (vxR ** 2 + vyR ** 2)

    # First compute the Roe Averages
    RT = (rhoR / rhoL) ** 0.5
    rho = RT * rhoL
    vx = (vxL + RT * vxR) / (1.0 + RT)
    vy = (vyL + RT * vyR) / (1.0 + RT)
    H = (HL + RT * HR) / (1.0 + RT)
    a = (gamma_m1 * (H - 0.5 * (vx * vx + vy * vy))) ** 0.5
    vn = vx * nx + vy * ny
    vt = vx * tx + vy * ty

    # Wave Strengths
    drho = rhoR - rhoL
    dp = pR - pL
    dvn = vnR - vnL
    dvt = vtR - vtL

    dV = [0.0, 0.0, 0.0, 0.0]
    dV[0] = (dp - rho * a * dvn) / (2.0 * a * a)
    dV[1] = rho * dvt / a
    dV[2] = drho - dp / (a * a)
    dV[3] = (dp + rho * a * dvn) / (2.0 * a * a)

    # Wave Speed
    ws = [0.0, 0.0, 0.0, 0.0]
    ws[0] = abs(vn - a)
    ws[1] = abs(vn)
    ws[2] = abs(vn)
    ws[3] = abs(vn + a)

    # Harten's Entropy Fix JCP(1983), 49, pp357-393:
    # only for the nonlinear fields.
    dws = [0.2, 0.0, 0.0, 0.2]
    if (ws[0] < dws[0]): ws[0] = 0.5 * (ws[0] * ws[0] / dws[0] + dws[0])
    if (ws[3] < dws[3]): ws[3] = 0.5 * (ws[3] * ws[3] / dws[3] + dws[3])

    # Right Eigenvectors
    Rv = zeros((4, 4))
    Rv[0, 0] = 1.0
    Rv[1, 0] = vx - a * nx
    Rv[2, 0] = vy - a * ny
    Rv[3, 0] = H - vn * a

    Rv[0, 1] = 0.0
    Rv[1, 1] = a * tx
    Rv[2, 1] = a * ty
    Rv[3, 1] = vt * a

    Rv[0, 2] = 1.0
    Rv[1, 2] = vx
    Rv[2, 2] = vy
    Rv[3, 1] = 0.5 * (vx * vx + vy * vy)

    Rv[0, 3] = 1.0
    Rv[1, 3] = vx + a * nx
    Rv[2, 3] = vy + a * ny
    Rv[3, 3] = H + vn * a

    # Dissipation Term
    diss = zeros(4)
    for i in range(4):
        for j in range(4):
            diss[i] += ws[j] * dV[j] * Rv[i, j]

    # Compute the flux.
    fL = zeros(4)
    fL[0] = rhoL * vnL
    fL[1] = fL[0] * vxL + pL * nx
    fL[2] = fL[0] * vyL + pL * ny
    fL[3] = fL[0] * HL

    fR = zeros(4)
    fR[0] = rhoR * vnR
    fR[1] = fR[0] * vxR + pR * nx
    fR[2] = fR[0] * vyR + pR * ny
    fR[3] = fR[0] * HR

    Roe = 0.5 * (fL + fR - diss) * side.area
    return Roe
