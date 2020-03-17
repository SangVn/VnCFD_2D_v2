# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

# Bài toán dòng chảy trong Nozzle

from lib.functions import Pm2P
from lib.fluxes import flux_roe_python
from lib.boco import *

# Trước hết cần chỉ rõ đường dẫn của file lưới
mesh_file = 'examples/nozzle/cdnozzle.mesh'

# Thiết lập thông số dòng chảy tới:
# Mach = 0.2, T = 293.15 K, p = 101325.0 Pa, góc tấn alfa = 0.0 độ
# Hàm Pm2P xác định các biến nguyên thủy (rho, u, v, p) từ các biến đầu vào trên
P_freestream = Pm2P(M=0.2, T=293.15, p=101325.0, alf=0.0)


# Đăt điều kiện biên: boco_i = [(name, start_side_index, end_side_index), (...)]
def set_boco():
    # hàm set_boco_const thiết lập thông số dòng chảy tới và áp suất p_exit tại mặt ra (nếu cần).
    # áp suất toàn phần: pt = 104190.585, xét 3 trường hợp p_exit/pt = 0.89, 0.6, 0.16
    set_boco_const(Pfree=P_freestream, pexit=104190.585 * 0.6)

    blk1_bc_0 = [(inflow, None, None)]  # Điều kiện biên bound_0 là inflow
    blk1_bc_1 = [(outflow, None, None)]  # Điều kiện biên bound_1 là outflow
    blk1_bc_2 = [(symmetry, None, None)]  # Điều kiện biên bound_2 là symmetry
    blk1_bc_3 = [(symmetry, None, None)]  # Điều kiện biên bound_3 là symmetry (hãy thử no_slip)
    # (..., Non, None) nghĩa là bắt đầu từ side đầu tiên tới side cuối cùng trên biên
    blk1_bc_list = [blk1_bc_0, blk1_bc_1, blk1_bc_2, blk1_bc_3]  # list các điều kiện biên của block 1

    boco_list = [blk1_bc_list]
    joint_list = None
    return boco_list, joint_list


boco_list, joint_list = set_boco()

# các thông số khác
CFL = 0.85
time_target = None
iter_target = 12000

# thời điểm ghi kết quả giữa chừng
write_field_frequency_time = None
write_field_frequency_iter = 5000

# thời điểm hiển thị các thông số cơ bản của một bước
print_frequency_iter = 500

# lựa chọn hàm tính flux, hàm flux_roe_fortran đã được include trong lib.boco
flux_func = flux_roe_fortran