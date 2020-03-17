# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

from lib.data import Blocks
# from setting import P_left, P_right # test1D

def run():
    blocks = Blocks()
    # blocks.init_field_test1D(P_left, P_right) # test1D
    blocks.init_field()
    blocks.run()
    blocks.export_block_data('field_bd.dat')
    blocks.plot_field(field='p', pfunc='pcolor')

if __name__ == '__main__':
    run()
