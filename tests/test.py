import unittest
import numpy as np
from spatial_math.SO3 import *
from spatial_math.base import Base

class BaseTest(unittest.TestCase):
    def test_quaternion_R(self):
        qtn = np.random.random(4)
        qtn = qtn/np.linalg.norm(qtn)
        R = Base._quaternion_to_R(qtn)
        qtn_ = Base._R_to_quaternion(R)
        R_ = Base._quaternion_to_R(qtn_)
        is_same = np.allclose(R, R_)
        self.assertTrue(is_same)

    def test_axisangle_R(self):
        qtn = np.random.random(4)
        qtn = qtn/np.linalg.norm(qtn)
        R = Base._quaternion_to_R(qtn)
        axis, angle = Base._R_to_axisangle(R)
        R_ = Base._axisangle_to_R(axis,angle)
        is_same = np.allclose(R, R_)
        self.assertTrue(is_same)
    
    def test_quaternion_axisangle(self):
        qtn = np.random.random(4)
        qtn = qtn/np.linalg.norm(qtn)
        axis, angle = Base._quaternion_to_axisangle(qtn)
        qtn_ = Base._axisangle_to_quaternion(axis, angle)
        axis_, angle_ = Base._quaternion_to_axisangle(qtn_)
        self.assertTrue(np.allclose(angle_, angle))
        self.assertTrue(np.allclose(axis_, axis))

class SO3Test(unittest.TestCase):
    def isSO3(self, SO3_):
        #it also tests @ operator
        I = np.eye(3)
        I_ = (SO3_ @ SO3_.inv()).R
        return np.allclose(I,I_)

    def test_construct_random(self):
        SO3_ = SO3.random()
        self.assertTrue(self.isSO3(SO3_))
    
    def test_construct_axisangle(self):
        axislist = np.eye(3)
        for omega in axislist:
            SO3_ = SO3.axisangle(omega,1)
        self.assertTrue(self.isSO3(SO3_))

    def test_construct_quaternion(self):
        sqrt2 = np.sqrt(2)
        SO3_ = SO3.qtn([sqrt2,sqrt2, 0, 0])
        self.assertTrue(self.isSO3(SO3_))

    def test_construct_Rxyz(self):
        SO3_ = SO3.Rx(np.pi/2)
        self.assertTrue(self.isSO3(SO3_))
        SO3_ = SO3.Ry(np.pi/2)
        self.assertTrue(self.isSO3(SO3_))
        SO3_ = SO3.Rx(np.pi/2)
        self.assertTrue(self.isSO3(SO3_))

    def test_inv(self):
        SO3_1 = SO3.random()
        SO3_2 = SO3_1.inv().inv()
        isEqual = np.allclose(SO3_1.R, SO3_2.R)
        self.assertTrue(isEqual)
    

# class SE3Test(unittest.TestCase):
#     def isSE3(self, T):
#         I = np.eye(4)
#         R = T.R
#         p = T.p
#         Tinv = np.block([[R.T, -R.T.dot(p[:,None])],
#                          [np.zeros([1,3]), 1]])
#         return np.allclose(I, T.T().dot(Tinv))

#     def test_construct_random(self):
#         T = SE3().random()
#         self.assertTrue(self.isSE3(T))

#     def test_construct_twistangle(self):
#         Slist = np.eye(6)
#         for S in Slist:
#             T = SE3().twistangle(S, 1)
#             self.assertTrue(self.isSE3(T))
    
#     def test_construct_qtn_trans(self):
#         qtn = [0.7071,0.7071, 0, 0]
#         trans = np.random.rand(3)
#         T = SE3().qtn_trans(qtn, trans)
#         self.assertTrue(self.isSE3(T))

#     def test_construct_Rxyz(self):
#         Rx = SE3().Rx(np.pi/2)
#         self.assertTrue(self.isSE3(Rx))
#         Ry = SE3().Ry(np.pi/2)
#         self.assertTrue(self.isSE3(Ry))
#         Rz = SE3().Rx(np.pi/2)
#         self.assertTrue(self.isSE3(Rz))
    
#     def test_construct_trans(self):
#         TransList = np.random.rand(5,3)
#         for trans in TransList:
#             T = SE3().trans(trans)
#             self.assertTrue(self.isSE3(T))

#     def test_inv(self):
#         T = SE3().random()
#         T_ = T.inv().inv()
#         isEqual = np.allclose(T.T(), T_.T())
#         self.assertTrue(isEqual)

#     def test_conversion_to_qtn_trans(self):
#         Tlist = []
#         Tlist.append(SE3().random()) #random R
#         Tlist.append(SE3().Rx(180,'deg')) #180
#         Tlist.append(SE3().Rz(0,'deg')) #0
#         Tlist.append(SE3().trans([3,2,4])) #pure translation

#         for T in Tlist:
#             qtn, trans = T.to_qtn_trans()
#             T_ = SE3().qtn_trans(qtn, trans)
#             isEqual = np.allclose(T.T(), T_.T())
#             self.assertTrue(isEqual)

#     def test_conversion_to_twistangle(self):
#         Tlist = []
#         Tlist.append(SE3().random()) #random R
#         Tlist.append(SE3().Rx(180,'deg')) #180
#         Tlist.append(SE3().Rz(0,'deg')) #0
#         Tlist.append(SE3().trans([3,2,4])) #pure translation

#         for T in Tlist:
#             twist, angle = T.to_twistangle()
#             T_ = SE3().twistangle(twist, angle)
#             isEqual = np.allclose(T.T(), T_.T())
#             self.assertTrue(isEqual)
    
#     def test_dot(self):
#         vlist = []
#         vlist.append(np.random.rand(4,1))
#         vlist.append(np.random.rand(3))
#         vlist.append(np.random.rand(3,1))
#         vlist.append(np.random.rand(4,4))
#         vlist.append(SE3())
#         for v in vlist:
#             SE3().dot(v)

if __name__ == '__main__':
    unittest.main()