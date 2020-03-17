# coding: utf-8
# Copyright (C) 2019  Nguyen Ngoc Sang, <https://github.com/SangVn>

from time import time as timer
from numpy import cross, zeros, array, fromfile
import matplotlib.pyplot as plt
from .functions import P2U, U2P, Mach, Temperature, import_mesh
from .solver import eu_solver
from setting import P_freestream, mesh_file, joint_list, boco_list, path_dir

'''
    ------------------------------------
    Phần I: Lớp dữ liệu "Cells"
    ------------------------------------
'''
def center(vertices):
    '''
    :param vertices: vertices: array tọa độ bốn đỉnh ô lưới theo thứ tự ngược chiều KĐH
    :return: center - tọa độ tâm ô lưới bằng trung bình cộng tọa độ bốn đỉnh lưới
    '''
    return sum(vertices) / 4.

def volume(vertices):
    '''
    :param vertices: array tọa độ bốn đỉnh ô lưới theo thứ tự ngược chiều KĐH
    :return: volume - thể tích ô lưới, trong trường hợp 2 chiều - diện tích ô lưới:
                      bằng một nửa độ lớn tích có hướng hai vector đường chéo
    '''
    return abs(cross(vertices[0] - vertices[2], vertices[1] - vertices[3]) / 2.)

def cell_size(vertices):
    '''
    :param vertices: array tọa độ bốn đỉnh ô lưới theo thứ tự ngược chiều KĐH
    :return: size  - kích thước ô lưới bằng đường trung bình nhỏ nhất
    '''
    dx_vec = (vertices[1]-vertices[0] + vertices[2] - vertices[3])/2.
    dy_vec = (vertices[2]-vertices[1] + vertices[3] - vertices[0])/2.
    dx = dx_vec.dot(dx_vec)**0.5
    dy = dy_vec.dot(dy_vec)**0.5
    return min(dx, dy)


class Cell:
    '''
    Lớp dữ liệu Cell chứa  các thông số cơ bản của một ô lưới.
    Parameters
    ----------
    vertices : 4 đỉnh lưới theo thứ tự ngược chiều KĐH dùng để xác định center, volume, size

    Attributes
    ----------
    center: tạo độ tâm
    volume: thể tích
    size:   kích thước
    P:   biến nguyên thủy (rho, u, v, p)
    U:   biến bảo toàn (rho, rho*u, rho*v, e), e = p/(gamma-1) + rho*(u^2 + v^2)/2
    res: biến tổng dòng qua các bề mặt sum(F*S)
    dt:  bước thời gian trong ô lưới.

    Notes
    -----
    Tại thời điểm khởi tạo P, U, res, dt bằng 0
    '''
    def __init__(self, vertices):
        'Khởi tạo Cell từ 4 đỉnh vertices.'
        self.center   = center(vertices)
        self.volume   = volume(vertices)
        self.size     = cell_size(vertices)
        self.P = zeros(4)
        self.U = zeros(4)
        self.res = zeros(4)
        self.dt  = 0.0


class Cells():
    '''
    Lớp dữ liệu các ô lưới trong một block (zone).

    Parameters
    ----------
    nodes : mảng 2D các tọa độ các điểm lưới có kích thước (Nj+1)*(Ni+1)

    Attributes
    ----------
    size:  kích thước lưới 2D ([Nj, Ni])
    len:   tổng số ô lưới (Nj*Ni)
    cells: dãy các ô lưới "class Cell".

    Notes
    -----
    Ô lưới thứ (j,i) (hay thứ j*i) gồm 4 đỉnh [(j,i), (j,i+1), (j+1,i+1), (j+1,i)]
    '''
    def __init__(self, nodes):
        '''Khởi tạo Cells từ mảng nodes.'''
        Nj, Ni = nodes.shape[0]-1, nodes.shape[1]-1
        self.size  = [Nj, Ni]
        self.len   = Nj*Ni
        self.cells = []
        for j in range(Nj):
            for i in range(Ni):
                vers = (nodes[j, i], nodes[j, i + 1], nodes[j + 1, i + 1], nodes[j + 1, i])
                self.cells.append(Cell(vers))

    def __getitem__(self, item):
        '''
        Các phương thức cơ bản để lấy các phần tử của dãy các ô lưới "cells":
        :param item: có thể là kiểu tuple (j, i), int j*i hay để getslice [start:stop]
        :return: một hoặc một đoạn các ô lưới
        '''
        if isinstance(item, tuple): # Lấy ô lưới thứ (j,i): cells[j,i]
            j, i = item
            if (j < 0): j += self.size[0]
            if (i < 0): i += self.size[1]
            return self.cells[j * self.size[1] + i]
        else: # Lấy một đoạn các ô lưới: cells[start:stop] (getslice). Lấy ô lưới thứ j*i: cells[j*i]
            return self.cells[item]

    def time_step_cell(selfs):
        '''Tính bước thời gian cục bộ trong từng ô lưới.'''
        for cell in selfs.cells:
            a = (1.4 * cell.P[3] / cell.P[0]) ** 0.5      # vận tốc âm thanh
            v = (cell.P[1] ** 2 + cell.P[2] ** 2) ** 0.5  # vận tốc dòng chảy
            cell.dt = cell.size/(v + a)

    def time_step_global(self, CFL):
        '''
        Xác định bước thời gian cho toàn vùng tính toán block (zone).
        :param CFL: số CFL, phụ thuộc sơ đồ tính toán.
        :return: dt - bước thời gian cho mỗi iteration
        '''
        self.time_step_cell() # trước hết cần xác định bước thời gian trong từng ô lưới
        dt = 1e6
        for cell in self.cells: # sau đó tìm bước thời gian nhỏ nhất trong toàn block
            dt = min(dt, cell.dt)          
        return CFL*dt

    def new_U(self, dt):
        '''Thực hiện bước lặp: xác định U ở bước thời gian tiếp theo.'''
        for cell in self.cells:
            cell.U += dt/cell.volume*cell.res # công thức: U^{n+1} = U^{n}  + dt/dx*RES
            cell.res[:] = 0.0                 # sau khi xác định U, đưa giá trị res về 0.0

    def new_P(self):
        '''Thực hiện bước lặp: xác định P ở bước thời gian tiếp theo, sử dụng hàm U2P.'''
        for cell in self.cells:
            U2P(cell.U, cell.P)

'''
    ------------------------------------
    Phần II: Lớp dữ liệu "Sides"
    ------------------------------------
'''

def area(side_vec):
    '''
    Tính diện tích mặt ô lưới, xét trong trường hợp 2D - chiều dài cạnh.
    :param side_vec: vector side = (point2 - point1)
    :return: chiều dài cạnh ô lưới
    '''
    return (side_vec.dot(side_vec))**0.5

def normal(side_vec):
    '''xác định vector pháp tuyến của bề mặt: S*n, với n - vector pháp tuyến đơn vị'''
    return array([side_vec[1], -side_vec[0]])

#định nghĩa lớp bề mặt
class Side:
    '''
    Lớp Side chứa các thông số cơ bản của một bề mặt của ô lưới.

    Parameters
    ----------
    side_vec : vector side = (point2 - point1)

    Attributes
    ----------
    area:   diện tích bề mặt, 2D - chiều dài cạnh
    normal: vector pháp tuyến đơn vị
    cells:  hai ô lưới hai bên trái phải.
            Quy ước: self.cells[0] - ô bên trái, self.cells[1] - ô bên phải
    '''
    def __init__(self, side_vec):
        '''Khởi tạo Side từ vector side_vec (node_2 - node_1).'''
        self.area   = area(side_vec)
        self.normal = normal(side_vec)/self.area
        self.cells = None


class Sides:
    '''
    Lớp dữ liệu các bề mặt trong một vùng tính toán block (zone).
    Sides gồm hai loại:
        + Side bên trong vùng tính toán
        + Side trên biên
    Quy ước:
        Mỗi vùng tính toán block có 4 biên:
                        bound_3
                       <--------
                   ^               ^
            bound_0|  inner_sides  |bound_1

                       <--------
                        bound_2
    Parameters
    ----------
    nodes : mảng chứa tọa độ các điểm lưới
    cells : dãy các ô lưới

    Attributes
    ----------
    bounds : list
            [bound_0, bound_1, bound_2, bound_3]
    inner_sides : list
            tất cả các sides bên trong
    boco_list : list
            các điều kiện biên trên 4 biên [boco_bound_0, boco_bound_1, boco_bound_2, boco_bound_3],
            mỗi boco_bound_i là một list các điều kiện biên trên biên i.
    '''

    def __init__(self, nodes, cells):
        '''Khởi tạo Sides từ nodes và cells.'''
        #xác định các vector bề mặt từ các điểm lưới 
        sides_i = nodes[:, :-1] - nodes[:, 1:]   # sides nằm ngang
        sides_j = nodes[1:] - nodes[:-1]         # sides thẳng đứng

        # Biên_0 gồm các mặt ở cột đầu sides_j
        # Ô lưới bên trái không xác định, bên phải là các ô ở cột đầu tiên
        bound_0 = []
        for j in range(sides_j.shape[0]):
            side = Side(sides_j[j, 0])
            side.cells = [None, cells[j, 0]]
            bound_0.append(side)

        # Biên_1 gồm các mặt ở cột cuối sides_j
        # Ô lưới bên phải không xác định, bên trái là các ô ở cột cuối
        bound_1 = []
        for j in range(sides_j.shape[0]):
            side = Side(sides_j[j, -1])
            side.cells = [cells[j, -1], None]
            bound_1.append(side)

        # Biên_2 gồm các mặt ở hàng đầu sides_i
        # Ô lưới bên phải không xác định, bên trái là các ô ở hàng đầu
        bound_2 = []
        for i in range(sides_i.shape[1]):
            side = Side(sides_i[0, i])
            side.cells = [None, cells[0, i]]
            bound_2.append(side)

        # Biên_3 gồm các mặt ở hàng cuối sides_i
        # Ô lưới bên trái không xác định, bên phải là các ô ở hàng cuối
        bound_3 = []
        for i in range(sides_i.shape[1]):
            side = Side(sides_i[-1, i])
            side.cells = [cells[-1, i], None]
            bound_3.append(side)

        self.bounds = [bound_0, bound_1, bound_2, bound_3]

        # Những hàng còn lại bên trong sides_i
        self.inner_sides = []
        for j in range(1, sides_i.shape[0] - 1):
            for i in range(sides_i.shape[1]):
                side = Side(sides_i[j, i])
                side.cells = [cells[j - 1, i], cells[j, i]]
                self.inner_sides.append(side)

        # Những cột còn lại bên trong sides_j
        for i in range(1, sides_j.shape[1] - 1):
            for j in range(sides_j.shape[0]):
                side = Side(sides_j[j, i])
                side.cells = [cells[j, i - 1], cells[j, i]]
                self.inner_sides.append(side)

    def set_boco_list(self, boco_list):
        '''Thiết lập giá trị thuộc tính boco_list.'''
        self.boco_list = boco_list

    def set_joint(self, boundary_0, boundary_1, ic=0):
        '''Kết nối các mặt với điều kiện biên joint.'''
        for side_0, side_1 in zip(boundary_0, boundary_1):
            side_0.cells[ic] = side_1.cells[ic]
            self.inner_sides.append(side_0)

    def flux_bound_sides(self):
        '''Tính dòng qua các side trên biên.'''
        for i, bound in enumerate(self.bounds):
            for boco in self.boco_list[i]: boco[0](bound[boco[1]:boco[2]], (i+1)%2)

    def flux_inner_sides(self, flux_func):
        '''
        Tính dòng qua các sides bên trong vùng tính.
        :param flux_func: hàm tính dòng qua side
        '''
        for side in self.inner_sides:
            F = flux_func(side, side.cells[0].P, side.cells[1].P)
            side.cells[0].res -= F # ô bên trái -
            side.cells[1].res += F # ô bên phải +

'''
    ------------------------------------
    Phần III: Lớp dữ liệu "Blocks"
    ------------------------------------
'''

class Block:
    '''
    Lớp dữ liệu Block chứa các dữ liệu của một vùng tính toán (block, zone).

    Parameters
    ----------
    name  : tên block
    nodes : mảng tọa độ các điểm lưới

    Attributes
    ----------
    BField : file chứa trường khí động 'name.field'
    BNodes : các điểm lưới
    BCells : các ô lưới "class Cells"
    BSides : các mặt "class Sides"
    BSize  : kích thước lưới
    '''
    def __init__(self, name, nodes):
        '''Khởi tạo Block có tên "name", có tọa độ điểm lưới "nodes".'''
        self.BField = path_dir+name+'.field'
        self.BNodes = nodes
        self.BCellS = Cells(nodes)
        self.BSides = Sides(nodes, self.BCellS)
        self.BSize = self.BCellS.size

    def iteration(self, flux_func, dt):
        '''
        Mỗi bước lặp bao gồm:
            set_boco() : xác định các Ghost Ceels
            flux_bound_sides() : tính dòng qua các sides trên biên
            flux_inner_sides() : tính dòng qua các sides trong vùng tính
            new_U : xác định U ở bước thời gian tiếp theo
            new_P : xác định P ở bước thời gian tiếp theo

        :param flux_func: hàm tính dòng qua side
        :param dt: bước thời gian
        '''
        self.BSides.flux_bound_sides()
        self.BSides.flux_inner_sides(flux_func)
        self.BCellS.new_U(dt)
        self.BCellS.new_P()

    def write_field(self):
        '''Ghi trường khí động dạng vào file binary BField.'''
        BData = array([cell.P for cell in self.BCellS])
        f = open(self.BField, 'wb')
        BData.tofile(f)
        f.close()

    def read_field(self):
        '''Đọc trường khí động từ file binary BField.'''
        f = open(self.BField, 'rb')
        BData = fromfile(f).reshape((self.BCellS.len, 4))
        f.close()
        for i, cell in enumerate(self.BCellS):
            cell.P = BData[i]
            cell.U = P2U(BData[i])

    def init_field(self, P_t0):
        '''
        Thiết lập trường khí động tại thời điểm ban đầu.
        :param P_t0 : numpy.array(rho, u, v, p)
        '''
        U_t0 = P2U(P_t0)
        for cell in self.BCellS:
            cell.P = P_t0.copy()
            cell.U = U_t0.copy()


class Blocks():
    '''
    Lớp dữ liệu Blocks chứa tất cả các blocks trong vùng tính toán và có chức năng thực hiên
    tất cả các bước tính toán khí động.

    Parameters
    ----------
    meshfile : file lưới ở định dạng Tecplot

    Attributes
    ----------
    len:
            số lượng các blocks
    blocks:
            dãy các blocks
    time_step_global:
            bước thời gian trong toàn bộ vùng tính toán
    state_file:
            file trạng thái tính toán, chứ hai thông số của bước lặp cuối cùng - iter, time
    '''
    def __init__(self, meshfile=mesh_file):
        '''Khởi tạo Blocks từ file lưới "meshfile".'''
        start_time = timer()
        self.state_file = path_dir+'solver.state'
        self.time_step_global = 1e6
        # số lượng block; tên block, tọa độ điểm lưới trong mỗi block
        zone_n, zone_names, zone_nodes = import_mesh(meshfile)
        self.len = zone_n
        self.blocks = []
        for n in range(zone_n):
            block = Block(zone_names[n], zone_nodes[n])
            self.blocks.append(block)
        print('The time taken by init_block is %f seconds!' % (timer() - start_time))

    def __getitem__(self, item):
        '''Lấy phần tử của dãy Blocks.'''
        return self.blocks[item]

    def read_field(self):
        '''Đọc trường khí động.'''
        start_time = timer()
        for block in self.blocks: block.read_field()
        print('The time taken by read_field is %f seconds!\n' % (timer() - start_time))


    def write_field(self):
        '''Ghi trường khí động.'''
        for block in self.blocks: block.write_field()

    def init_field_test1D(self, P_left, P_right):
        '''Thiết lập điều kiện ban đầu cho test1D, lưới 1 block.'''
        start_time = timer()
        with open(self.state_file, 'w') as f: f.write('iter time:\n%d %f' % (0, 0.0))

        U_left  = P2U(P_left)
        U_right = P2U(P_right)

        cells = self.blocks[0].BCellS
        Ni = cells.len
        Nih = int(Ni/2)
        for cell in cells[: Nih]:
            cell.P = P_left.copy()
            cell.U = U_left.copy()
        for cell in cells[Nih:]:
            cell.P = P_right.copy()
            cell.U = U_right.copy()

        self.write_field()
        print('The time taken by init_field is %f seconds!' % (timer() - start_time))

    def init_field(self, P_t0=P_freestream):
        '''Thiết lập điều kiện ban đầu, ghi trường khí động.'''
        start_time = timer()
        with open(self.state_file, 'w') as f: f.write('iter time:\n%d %f' % (0, 0.0))
        for block in self.blocks: block.init_field(P_t0)
        self.write_field()
        print('The time taken by init_field is %f seconds!' % (timer() - start_time))

    def iteration(self, flux_func):
        '''Thực hiện bước lặp thời gian.'''
        for block in self.blocks: block.iteration(flux_func, self.time_step_global)

    def joint(self, joints = joint_list):
        '''
        Kết nối các biên có điều kiện biên joint.
        :param  joint_list : [joint_0, joint_1, ...]
                joint_0 = [blk1_id, bound1_id, start_side1_id, end_side1_id,
                           blk2_id, bound2_id, start_side2_id, end_side2_id]
        '''
        if joints is not None:
            for joint in joints:
                left_BSides = self.blocks[joint[0]].BSides
                left_bound = left_BSides.bounds[joint[1]][joint[2]:joint[3]]

                right_BSides = self.blocks[joint[4]].BSides
                right_bound = right_BSides.bounds[joint[5]][joint[6]:joint[7]]

                out_cell_id = joint[1]%2
                left_BSides.set_joint(left_bound, right_bound, ic=out_cell_id)

    def set_time_step(self, CFL):
        '''Xác định bước thời gian trong toàn bộ vùng tính.'''
        for block in self.blocks:
            dt = block.BCellS.time_step_global(CFL)
            self.time_step_global = min(self.time_step_global, dt)
        return self.time_step_global

    def run(self, solver=eu_solver, joints=joint_list, bocos=boco_list):
        '''
        Toàn bộ các bước tính toán trường khí động:
                        joint -> set_boco_list -> read_field -> solve
        :param solver : hàm solver thực hiện các bước lặp  Block.iteration.
        '''
        print('**************************************************')
        start_time = timer()
        self.joint(joints)
        for i, block in enumerate(self.blocks): block.BSides.set_boco_list(bocos[i])
        self.read_field()
        solver(self)
        print('\nThe time taken by solver is %f seconds!' % (timer() - start_time))
        print('**************************************************')

    def export_block_data(self, filename):
        '''Xuất trường khí động vào file định dạng Tecplot block data.'''
        filename = path_dir + filename
        print('Write block data to: %s' % filename)
        f = open(filename, 'w')
        f.write('TITLE = "vncfd field"\n')
        f.write('VARIABLES = "X", "Y", "rho", "u", "v", "p", "Mach", "T"\n')

        n = 0
        for block in self.blocks:
            n += 1
            nodes = block.BNodes
            cells = block.BCellS
            f.write('ZONE T="%d", I=%d, J=%d, DATAPACKING=BLOCK, VARLOCATION=([3,4,5,6,7,8]=CELLCENTERED)\n' % (
            n, nodes.shape[1], nodes.shape[0]))

            X_p, Y_p = nodes[:, :, 0].ravel(), nodes[:, :, 1].ravel()
            for x in X_p: f.write('%f ' % x)
            f.write('\n')
            for y in Y_p: f.write('%f ' % y)
            f.write('\n')

            for i in range(4):
                for cell in cells: f.write('%f ' % cell.P[i])
                f.write('\n')
            for cell in cells:
                M = Mach(cell.P)
                f.write('%f ' % M)
            for cell in cells:
                T = Temperature(cell.P)
                f.write('%f ' % T)
            f.write('\n')
        f.close()

    def plot_field(self, field='rho', pfunc='pcolor'):
        '''Plot trường khí động. Tốt nhất hãy sử dụng Paraview!'''
        id = {'rho':0, 'u':1, 'v':2, 'p':3}
        for block in self.blocks:
            nodes = block.BNodes
            cells = block.BCellS
            Nj, Ni = cells.size[0], cells.size[1]
            value_c = array([cell.P[id[field]] for cell in cells]).reshape((Nj, Ni))

            if pfunc == 'pcolor':
                X_p, Y_p = nodes[:, :, 0], nodes[:, :, 1]
                pcm = plt.pcolor(X_p, Y_p, value_c)
            else:
                centers = array([cell.center for cell in cells]).reshape((Nj, Ni, 2))
                X_c, Y_c = centers[:, :, 0], centers[:, :, 1]
                pcm = plt.contourf(X_c, Y_c, value_c)

        plt.title(field)
        plt.show()