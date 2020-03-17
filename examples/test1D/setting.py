# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

from lib.boco import *
from numpy import array

mesh_file = 'examples/test1D/test1D.mesh'

P_left  = array([1.0, 0.0, 0.0, 101325.0])
P_right = array([0.1, 0.0, 0.0, 10132.50])

# điều kiện biên: boco = [(name, start_index, end_index), (...)]
def set_boco():
    blk1_bc_0 = [(supersonic_outflow, None, None)]
    blk1_bc_1 = [(supersonic_outflow, None, None)]
    blk1_bc_2 = [(null, None, None)]
    blk1_bc_3 = [(null, None, None)]
    blk1_bc_list  = [blk1_bc_0, blk1_bc_1, blk1_bc_2, blk1_bc_3]

    boco_list = [blk1_bc_list]
    joint_list = None
    return boco_list, joint_list

boco_list, joint_list = set_boco()

# các thông số khác
CFL = 1.0
time_target = 0.0005
iter_target = None

# thời điểm ghi kết quả giữa chừng
write_field_frequency_time = None
write_field_frequency_iter = None

# thời điểm hiển thị các thông số cơ bản của một bước
print_frequency_iter = 100

# lựa chọn hàm tính flux
flux_func = flux_roe_fortran