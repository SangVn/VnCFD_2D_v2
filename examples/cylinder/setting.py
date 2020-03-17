# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

from lib.fluxes import flux_roe_python
from lib.boco import *
from lib.functions import Pm2P

# Cylinder
mesh_file = 'examples/cylinder/cylinder.mesh'
P_freestream = Pm2P(M=0.2, T=293.15, p=101325, alf=0.0) # pt = 104190.585

# điều kiện biên
def boco():
    set_boco_const(P_freestream)
    boco_0 = [(joint, None, None)]
    boco_1 = [(joint, None, None)]
    boco_2 = [(no_slip, None, None)]
    boco_3 = [(farfield, None, None)]
    boco_list = [[boco_0, boco_1, boco_2, boco_3]]
    joint_list = [(0, 0, None, None, 0, 1, None, None)]
    return boco_list, joint_list

boco_list, joint_list = boco()

# các thông số khác
CFL = 1.0
time_target = None
iter_target = 10000

# thời điểm ghi kết quả giữa chừng
write_field_frequency_time = None
write_field_frequency_iter = None

# thời điểm hiển thị các thông số cơ bản của một bước
print_frequency_iter = 1000

# lựa chọn hàm tính flux
flux_func = flux_roe_fortran