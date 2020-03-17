# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

import setting as set
# import example_setting_outflow as set

# Hàm eu_solver thực hiện các bước lặp để tìm nghiệm
# Biến đầu vào gồm có: các ô lưới, các mặt, số vòng lặp, thời gian lúc ban đầu
def eu_solver(blocks):
    # đọc iter và time từ file status
    with open(blocks.state_file, 'r') as f: line = f.readlines()[1].split()
    iter = int(line[0])
    time = float(line[1])

    # tính theo thời gian
    if set.time_target is None: set.time_target = 1.e10
    if set.iter_target is None: set.iter_target = 1e10
    if set.write_field_frequency_time is None: set.write_field_frequency_time = 1.e10
    if set.write_field_frequency_iter is None: set.write_field_frequency_iter = 1e10
    if set.print_frequency_iter is None: set.display_iter = 1
    while(time < set.time_target and iter < set.iter_target):
        iter += 1
        # tính bước thời gian
        dt = blocks.set_time_step(set.CFL)  # có thể thiết lập dt trong set: dt = set.dt
        if (time + dt > set.time_target): dt = set.time_target - time
        time += dt

        blocks.iteration(set.flux_func) # gọi bước lặp
        if not(iter % set.print_frequency_iter): print('iteration: %d, dt: %f, time: %f' % (iter, dt, time))

        # ghi lại kết quả giữa chừng
        period = 1
        if (iter ==  period*set.write_field_frequency_iter) or (time >= period*set.write_field_frequency_time > time-dt):
            period += 1
            print('\nwrite_field at iteration: %d, time: %f\n' % (iter, time))
            blocks.write_field()

    # ghi lại kết quả cuối cùng
    blocks.write_field()
    with open(blocks.state_file, 'w') as f: f.write('iter time:\n%d %f' % (iter, time))


