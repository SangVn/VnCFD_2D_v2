# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

from lib.fluxes import flux_roe_python
from lib.boco import *
from lib.functions import Pm2P

# wedge
mesh_file = 'examples/wedge/wedge.mesh'
P_freestream = Pm2P(M=2.0, T=293.15, p=101325, alf=0.0)

# điều kiện biên: boco = [(name, start_index, end_index), (...)]
def boco_2blk():
    set_boco_const(Pfree=P_freestream)
    blk1_bc_0 = [(supersonic_inflow, None, None)]
    blk1_bc_1 = [(joint, None, None)]
    blk1_bc_2 = [(symmetry, None, 30), (no_slip, 30, None)]
    blk1_bc_3 = [(symmetry, None, None)]
    blk1_bc_list  = [blk1_bc_0, blk1_bc_1, blk1_bc_2, blk1_bc_3]

    blk2_bc_0 = [(joint, None, None)]
    blk2_bc_1 = [(supersonic_outflow, None, None)]
    blk2_bc_2 = [(no_slip, None, None)]
    blk2_bc_3 = [(symmetry, None, None)]
    blk2_bc_list  = [blk2_bc_0, blk2_bc_1, blk2_bc_2, blk2_bc_3]

    #[left block, left bound, id_start, id_end,
    # right block, right bound, id_start, id_end]
    joint_list = [(0, 1, None, None, 1, 0, None, None)]
    boco_list = [blk1_bc_list, blk2_bc_list]
    return boco_list, joint_list

boco_list, joint_list = boco_2blk()


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