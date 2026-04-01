import math
from enum import Enum
from scipy.optimize import fsolve

import scipy.special as ellip
import scipy.integrate as inte
import scipy.optimize as opt
import numpy


# 决定需要的数据类型
class PlotOption(Enum):
    X_Field = 0,
    Y_Field = 1,
    Z_Field = 2,
    NV_Field = 3,
    Intensity_Field = 4


class Magnet:
    r0 = 1
    l0 = 1
    theta = 0
    phi = 0
    Magnetization = 1
    __Miu = 4 * numpy.pi * 10 ** (-7)

    def __init__(self, r0, l0, theta, phi, magnetization):
        self.r0 = r0
        self.l0 = l0
        self.theta = theta
        self.phi = phi
        self.Magnetization = magnetization

    # 获取对应PlotOption类型的
    def _GetPlotData(self, B: list, dataType: PlotOption):
        # NV分量
        if dataType == PlotOption.NV_Field:
            return B[0] * numpy.cos(self.phi) * numpy.sin(self.theta) + B[1] * numpy.sin(self.phi) * numpy.sin(
                self.theta) + B[2] * numpy.cos(self.theta)
        # X分量
        if dataType == PlotOption.X_Field:
            return B[0]
        # Y分量
        if dataType == PlotOption.Y_Field:
            return B[1]
        # Z分量
        if dataType == PlotOption.Z_Field:
            return B[2]
        # 强度
        if dataType == PlotOption.Intensity_Field:
            return numpy.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)

    # 第一类椭圆积分
    @staticmethod
    def __GetE(k):
        return ellip.ellipe(k)

    # 第二类椭圆积分
    @staticmethod
    def __GetK(k):
        return ellip.ellipk(k)

    @staticmethod
    def __Getrou2(x, y):
        return x * x + y * y

    @staticmethod
    def __Getr2(x, y, z, h):
        return x * x + y * y + (z - h) * (z - h)

    @staticmethod
    def __Getalpha2(r0, x, y, z, h):
        return r0 ** 2 + Magnet.__Getr2(x, y, z, h) - 2 * r0 * numpy.sqrt(Magnet.__Getrou2(x, y))

    @staticmethod
    def __Getbeta2(r0, x, y, z, h):
        return Magnet.__Getalpha2(r0, x, y, z, h) + 4 * r0 * numpy.sqrt(Magnet.__Getrou2(x, y))

    @staticmethod
    def __Getk2(r0, x, y, z, h):
        return 1 - Magnet.__Getalpha2(r0, x, y, z, h) / Magnet.__Getbeta2(r0, x, y, z, h)

    def __getiG(self, x, y, z, h):
        alpha2 = Magnet.__Getalpha2(self.r0, x, y, z, h)
        beta2 = Magnet.__Getbeta2(self.r0, x, y, z, h)
        rou2 = Magnet.__Getrou2(x, y)
        r2 = Magnet.__Getr2(x, y, z, h)
        e = Magnet.__GetE(Magnet.__Getk2(self.r0, x, y, z, h))
        k = Magnet.__GetK(Magnet.__Getk2(self.r0, x, y, z, h))
        return ((self.r0 ** 2 + r2) * e - alpha2 * k) * (z - h) / (2 * alpha2 * numpy.sqrt(beta2) * rou2)

    def __getiBz(self, x, y, z, h):
        alpha2 = Magnet.__Getalpha2(self.r0, x, y, z, h)
        beta2 = Magnet.__Getbeta2(self.r0, x, y, z, h)
        rou2 = Magnet.__Getrou2(x, y)
        r2 = Magnet.__Getr2(x, y, z, h)
        e = Magnet.__GetE(Magnet.__Getk2(self.r0, x, y, z, h))
        k = Magnet.__GetK(Magnet.__Getk2(self.r0, x, y, z, h))
        return ((self.r0 ** 2 - r2) * e + alpha2 * k) / (2 * alpha2 * numpy.sqrt(beta2))

    # 得到目标磁场(磁感应强度沿Z轴正向,单位为高斯)坐标定义：磁铁堆成轴为z
    def GetField(self, x, y, z):
        bG = inte.quad(lambda h: self.__getiG(x, y, z, h), -self.l0 / 2, self.l0 / 2)[0]
        bx = x * bG
        by = y * bG
        bz = inte.quad(lambda h: self.__getiBz(x, y, z, h), -self.l0 / 2, self.l0 / 2)[0]
        temp = self.__Miu / numpy.pi * self.Magnetization
        return bx * temp, by * temp, bz * temp

    def GetNVField(self, x, y, z):
        B = self.GetField(x, y, z)
        return self._GetPlotData(B, PlotOption.NV_Field)

    def GetIntensityBField(self, x, y, z):
        B = self.GetField(x, y, z)
        return self._GetPlotData(B, PlotOption.Intensity_Field)

    def GetG(self, x, y, z):
        temp = self.__Miu / numpy.pi * self.Magnetization
        return temp * inte.quad(lambda h: self.__getiG(x, y, z, h), -self.l0 / 2, self.l0 / 2)[0]

    def GetBZ(self, rou, z):
        B = self.GetField(rou, 0, z)
        return B[2]

    # 获取在rou处，z=0的G的一阶导数(数值解)
    def GetD1GAtz0(self, x, y):
        det = 1e-5
        vp1 = self.GetG(x, y, det)
        vn1 = self.GetG(x, y, -det)
        vp2 = self.GetG(x, y, 2 * det)
        vn2 = self.GetG(x, y, -2 * det)
        return (-vp2 + 8 * vp1 - 8 * vn1 + vn2) / (12 * det)

    def GetD1NVAtz0(self, x, y):
        det = 1e-5
        vp1 = self.GetNVField(x, y, det)
        vn1 = self.GetNVField(x, y, -det)
        vp2 = self.GetNVField(x, y, 2 * det)
        vn2 = self.GetNVField(x, y, -2 * det)
        return (-vp2 + 8 * vp1 - 8 * vn1 + vn2) / (12 * det)

    def GetD2NVAtz0(self, x, y):
        det = 1e-5
        v1 = self.GetNVField(x, y, det)
        v2 = self.GetNVField(x, y, -det)
        v0 = self.GetNVField(x, y, 0)
        return (v2 + v1 - 2 * v0) / det ** 2

    def GetD3GAtz0(self, x, y):
        det = 0.0001
        vp1 = self.GetG(x, y, det)
        vp2 = self.GetG(x, y, 2 * det)
        vp3 = self.GetG(x, y, 3 * det)
        vn1 = self.GetG(x, y, -det)
        vn2 = self.GetG(x, y, -2 * det)
        vn3 = self.GetG(x, y, -3 * det)
        return (-vp3 + 8 * vp2 - 13 * vp1 + 13 * vn1 - 8 * vn2 + vn3) / (8 * det ** 3)

    def GetD2BzAtz0(self, rou):
        det = 1e-5
        v1 = self.GetBZ(rou, det)
        v2 = self.GetBZ(rou, - det)
        v0 = self.GetBZ(rou, 0)
        return (v2 + v1 - 2 * v0) / det ** 2

    def GetD1RouBzAtz0(self, rou):
        det = 1e-5
        v1 = self.GetBZ(rou - det, 0)
        v2 = self.GetBZ(rou + det, 0)
        return (v2 - v1) / (2 * det)


def GetPillarField(radius, length, MIntensity, x, y, z):
    M = Magnet(radius, length, 0, 0, MIntensity)
    res = M.GetField(x, y, z)
    return [res[0], res[1], res[2]]


def GetPillarNVField(radius, length, theta, phi, MIntensity, x, y, z):
    M = Magnet(radius, length, theta, phi, MIntensity)
    res = M.GetNVField(x, y, z)
    return res


def GetPillarIntensityField(radius, length, MIntensity, x, y, z):
    M = Magnet(radius, length, 0, 0, MIntensity)
    res = M.GetNVField(x, y, z)
    return math.sqrt(res[0] ** 2 + res[1] ** 2 + res[2] ** 2)

#磁铁和NV在同一平面，坐标系：x为磁铁正向
def GetAngle(M, x, dis):
    r = list()
    for i in x:
        B = M.GetField(0, dis, i)
        an = math.atan2(math.sqrt(B[0] ** 2 + B[2] ** 2), B[1])
        an = an * 180 / math.pi
        r.append(an)
    return numpy.array(r)

# 指定：方向（targetthe，targetphi），距离d。默认：phi方向不变。计算：磁铁x，y
# 指定方向和给定的Z方向距离，找出对应方向上的磁场，相对于角度基准点（磁铁朝向为X轴）要转动的角度，以及XY平面上相对于原点的坐标
def FindDire(r0, l0, targetthe, targetphi, distance):
    M = Magnet(r0, l0, 0, 0, 1)
    # 计算以X轴为正向的给定theta角的距离
    res = opt.root(lambda x: GetAngle(M, x, distance) - targetthe, numpy.array([0]), tol=1e-5)
    r = res.x[0]
    B = M.GetField(0, distance, r)
    # 计算detphi。相对NV动磁铁，所以反方向-180度
    detphi = targetphi - 180
    while detphi > 360:
        detphi -= 360

    while detphi < -360:
        detphi += 360

    # 返回：转的phi，磁铁的x坐标，y坐标，磁场计算值
    return [detphi, r * math.cos(detphi / 180 * math.pi), r * math.sin(detphi / 180 * math.pi), math.sqrt(
        B[0] ** 2 + B[1] ** 2 + B[2] ** 2)]

#2025.12.21新增
# 对任意坐标xyz求磁场θφ（坐标系：x为磁铁正向）
######问题：为什么之前有list？？？？？？？？？？？？？？？？？？？？？？？？？/
def GetThePhi(r0,l0, x,y,z,theta1):#theta1是磁铁方向
    M = Magnet(r0, l0, 0, 0, 1)
    #坐标系转换：实验室坐标系(x,y,z)转到磁铁坐标系(xM,yM,z)
    xM = x*math.cos(theta1/ 180 * math.pi)+y*math.sin(theta1/ 180 * math.pi)
    yM = -x*math.sin(theta1/ 180 * math.pi)+y*math.cos(theta1/ 180 * math.pi)
    #在磁铁坐标系计算(thetaM,phiM)
    BM = M.GetField(yM,z, xM)
    thetaM = math.atan2(math.sqrt(BM[0] ** 2 + BM[2] ** 2), BM[1])* 180 / math.pi
    phiM = math.atan2(BM[0], BM[2])* 180 / math.pi
    # 坐标系转换：磁铁坐标系(thetaM,phiM)转到实验室坐标系(thetaM,phiM)
    theta = thetaM
    phi = phiM+theta1
    if phi < 0:phi+=360
    return (theta, phi)

#输入phi，输出x,y
#r0, l0, targetthe, targetphi, distance
def FindDireFixedPhi(r0, l0, targetTheta, targetPhi, distance, magnetPhi):
    M = Magnet(r0, l0, 0, 0, 1)
    targetPhi=targetPhi%360
    def equations(vars):
        x, y = vars
        theta_val, phi_val = GetThePhi(r0,l0, x, y, distance,magnetPhi)

        # 返回与目标值的差
        return (theta_val - targetTheta, phi_val - targetPhi)

    # 初始猜测值
    initial_guessArray = [[5, 5],[-5,-5],[5,-5],[-5,5]]
    # 求解
    for initial_guess in initial_guessArray:
        solution = fsolve(equations, initial_guess)
        x_solution, y_solution = solution
        if abs(x_solution)<20 and abs(y_solution)<20:
            B = M.GetField(y_solution, distance, x_solution)
            return (x_solution, y_solution, math.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2))
    # 如果都解不出来，返回nan
    return (float('nan'), float('nan'), float('nan'))

# 根据Z方向测到的两个值估算Z方向的距离
def MagFun(M, x, dist):
    res = list()
    for i in x:
        res.append(M.GetNVField(i, 0, 0) / M.GetNVField(i + dist, 0, 0))
    return numpy.array(res)


def FindRoot(r0, l0, z1, z2, value):
    miu = 4 * math.pi * 10 ** (-7)
    MIn = 1.19 / miu
    M = Magnet(r0, l0, 0, 0, MIn / miu)
    xs = [2 * r0]
    res = opt.root(lambda x: MagFun(M, x, abs(z1 - z2)) - value, numpy.array(xs), tol=1e-5)
    return res.x[0]



"""
    测试代码
    绘制 x 和 y 随 phi 变化的曲线

    参数:
        phi_range: phi 的范围 (起始值, 结束值)，单位：度
        num_points: 点的数量
        figsize: 图形大小
    """
