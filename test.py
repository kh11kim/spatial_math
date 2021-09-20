import unittest
import numpy as np
from spatial_math import *

class SO3Test(unittest.TestCase):
    def isSO3(self, R):
        I = np.eye(3)
        I_ = R.R.dot(R.R.T)
        return np.allclose(I,I_)

    def test_construct_random(self):
        R = SO3().random()
        self.assertTrue(self.isSO3(R))

    def test_construct_axisangle(self):
        axislist = np.eye(3)
        for omega in axislist:
            R = SO3().axisangle(omega,1)
        self.assertTrue(self.isSO3(R))
    
    def test_construct_quaternion(self):
        R = SO3().qtn([0.7071,0.7071, 0, 0])

    def test_construct_Rxyz(self):
        Rx = SO3().Rx(np.pi/2)
        self.assertTrue(self.isSO3(Rx))
        Ry = SO3().Ry(np.pi/2)
        self.assertTrue(self.isSO3(Ry))
        Rz = SO3().Rx(np.pi/2)
        self.assertTrue(self.isSO3(Rz))

    def test_inv(self):
        R = SO3().random()
        R_ = R.inv().inv()
        isEqual = np.allclose(R.R, R_.R)
        self.assertTrue(isEqual)

    def test_conversion_to_quaternion(self):
        Rlist = []
        Rlist.append(SO3().random()) #random R
        Rlist.append(SO3().Rx(180,'deg')) #180
        Rlist.append(SO3().Rz(0,'deg')) #0

        for R in Rlist:
            R_ = SO3().qtn(R.to_qtn())
            isEqual = np.allclose(R.R, R_.R)
            self.assertTrue(isEqual)

    def test_conversion_to_axisangle(self):
        Rlist = []
        Rlist.append(SO3().random()) #random R
        Rlist.append(SO3().Rx(180,'deg')) #180
        Rlist.append(SO3().Rz(0,'deg')) #0

        for R in Rlist:
            axis, angle = R.to_axisangle()
            R_ = SO3().axisangle(axis, angle)
            isEqual = np.allclose(R.R, R_.R)
            self.assertTrue(isEqual)

    def test_dot(self):
        vlist = []
        vlist.append(np.random.rand(3,1))
        vlist.append(np.random.rand(3,3))
        vlist.append(SO3())
        for v in vlist:
            SO3().dot(v)

class SE3Test(unittest.TestCase):
    def isSE3(self, T):
        I = np.eye(4)
        R = T.R
        p = T.p
        Tinv = np.block([[R.T, -R.T.dot(p[:,None])],
                         [np.zeros([1,3]), 1]])
        return np.allclose(I, T.T().dot(Tinv))

    def test_construct_random(self):
        T = SE3().random()
        self.assertTrue(self.isSE3(T))

    def test_construct_twistangle(self):
        Slist = np.eye(6)
        for S in Slist:
            T = SE3().twistangle(S, 1)
            self.assertTrue(self.isSE3(T))
    
    def test_construct_qtn_trans(self):
        qtn = [0.7071,0.7071, 0, 0]
        trans = np.random.rand(3)
        T = SE3().qtn_trans(qtn, trans)
        self.assertTrue(self.isSE3(T))

    def test_construct_Rxyz(self):
        Rx = SE3().Rx(np.pi/2)
        self.assertTrue(self.isSE3(Rx))
        Ry = SE3().Ry(np.pi/2)
        self.assertTrue(self.isSE3(Ry))
        Rz = SE3().Rx(np.pi/2)
        self.assertTrue(self.isSE3(Rz))
    
    def test_construct_trans(self):
        TransList = np.random.rand(5,3)
        for trans in TransList:
            T = SE3().trans(trans)
            self.assertTrue(self.isSE3(T))

    def test_inv(self):
        T = SE3().random()
        T_ = T.inv().inv()
        isEqual = np.allclose(T.T(), T_.T())
        self.assertTrue(isEqual)

    def test_conversion_to_qtn_trans(self):
        Tlist = []
        Tlist.append(SE3().random()) #random R
        Tlist.append(SE3().Rx(180,'deg')) #180
        Tlist.append(SE3().Rz(0,'deg')) #0
        Tlist.append(SE3().trans([3,2,4])) #pure translation

        for T in Tlist:
            qtn, trans = T.to_qtn_trans()
            T_ = SE3().qtn_trans(qtn, trans)
            isEqual = np.allclose(T.T(), T_.T())
            self.assertTrue(isEqual)

    def test_conversion_to_twistangle(self):
        Tlist = []
        Tlist.append(SE3().random()) #random R
        Tlist.append(SE3().Rx(180,'deg')) #180
        Tlist.append(SE3().Rz(0,'deg')) #0
        Tlist.append(SE3().trans([3,2,4])) #pure translation

        for T in Tlist:
            twist, angle = T.to_twistangle()
            T_ = SE3().twistangle(twist, angle)
            isEqual = np.allclose(T.T(), T_.T())
            self.assertTrue(isEqual)
    
    def test_dot(self):
        vlist = []
        vlist.append(np.random.rand(4,1))
        vlist.append(np.random.rand(3))
        vlist.append(np.random.rand(3,1))
        vlist.append(np.random.rand(4,4))
        vlist.append(SE3())
        for v in vlist:
            SE3().dot(v)

if __name__ == '__main__':
    unittest.main()