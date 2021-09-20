#!/usr/bin/python 
# -*- coding: utf-8 -*-
"""spatial_math.py

This module demonstrates documentation as specified by the `Google Python
Style Guide`_. Docstrings may extend over multiple lines. Sections are created
with a section header and a colon followed by a block of indented text.

Todo:
    * For module TODOs
    * You have to also use ``sphinx.ext.todo`` extension

"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

def frame_plot(R, p, frame=None, dims=[0,10]*3, axeslength=1, axeswidth=5):
    """[summary]

    Args:
        R ([type]): [description]
        p ([type]): [description]
        frame ([type], optional): [description]. Defaults to None.
        dims ([type], optional): [description]. Defaults to [0,10]*3.
        axeslength (int, optional): [description]. Defaults to 1.
        axeswidth (int, optional): [description]. Defaults to 5.
    """
    # is_3d_plot = False
    # if len(plt.get_fignums()) != 0:
    #     fig = plt.gcf()
    #     if len(fig.axes) != 0:
    #         ax = plt.gca()
    #         if ax.__class__.__name__ == "Axes3DSubplot":
    #             is_3d_plot = True
    
    # if not is_3d_plot:
    #     fig = plt.figure()
    #     ax = plt.axes(projection='3d')    
    
    # xaxis = np.vstack([p, p+R[:3,0]]).T * axeslength
    # yaxis = np.vstack([p, p+R[:3,1]]).T * axeslength
    # zaxis = np.vstack([p, p+R[:3,2]]).T * axeslength

    # ax.plot(*xaxis, 'r', linewidth=axeswidth)
    # ax.plot(*yaxis, 'g', linewidth=axeswidth)
    # ax.plot(*zaxis, 'b', linewidth=axeswidth)

    # ax.set_xlim3d(dims[:2])
    # ax.set_ylim3d(dims[2:4])
    # ax.set_zlim3d(dims[4:])
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    # ax.set_zlabel('z')

    # frame_name_offset = [-0.5, -0.5, -0.5]
    # frame_name_pos = R@np.array(frame_name_offset).T+p
    # if frame is not None:
    #     ax.text(*frame_name_pos, "{{{}}}".format(frame))



class base:
    """[summary]
    Collection of Basic functions.
    """
    def __init__(self):
        self.R = np.eye(3)
        self.p = np.zeros(3)

    def _check_vector(self, vector_len, v):
        """[summary]
        Check if vector is specified shape.
        List will be converted to a np.ndarray
        Args:
            vector_len (int): The length of a vector.
            v (list or np.ndarray): target vector.

        Returns:
            v [np.ndarray]: vector.
        """
        if type(v) is list:
            v = np.array(v)
        elif (type(v) is np.ndarray)\
                & (v.shape == (vector_len,1)):
            v = v.flatten()
        assert v.shape == (vector_len,)
        return v

    def _check_matrix(self, matrix_len, m):
        """[summary]
        Check if matrix is specified shape.

        Args:
            matrix_len (int): The one-side length of a matrix.
            m (list or np.ndarray): target matrix.

        Returns:
            m [np.ndarray]: matrix.
        """
        assert m.shape == (matrix_len, matrix_len)
        return m
    
    def _skew(self, x):
        """[summary]
        Make a skew-symmetric matrix of a vector.
        Args:
            x (list or np.ndarray): vector

        Returns:
            [np.ndarray]: skew-symmetric matrix of a vector.
        """
        x = self._check_vector(3, x)
        return np.array([[0    , -x[2], x[1] ],
                        [x[2] , 0    , -x[0]],
                        [-x[1], x[0] , 0    ]])

    def plot(self, **kwargs):
        """[summary]
        """
        # frame_plot(self.R, self.p, **kwargs)

class SO3(base):
    """[summary]
    SO3 class
    """
    def __init__(self, *args):
        """[summary]
        Initialize SO3 Class
        Raises:
            Exception: [description]
        """
        if len(args) == 0:
            self.R = np.eye(3) # no input
        elif (len(args) == 1) \
                & (type(args[0]) is np.ndarray):
            if args[0].shape == (3,3):
                self.R = args[0] #input is 3x3 matrix
            else:
                Exception("SO3 Initialize Error!")
        else:
            raise Exception("SO3 Initialize Error!")
        
        # this is just for compatibility for plotting.
        self.p = np.zeros(3) 
    
    # --Construction of Rotation Matrix--
    def random(self):
        """[summary]
        Random Generation of SO3
        using QR decomposition

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        q, _ = np.linalg.qr(np.random.randn(3,3))
        return SO3(q)

    def axisangle(self, axis, angle):
        """[summary]
        Generate SO3 Class using an AxisAngle.

        Args:
            axis (size=3, np.ndarray or list): 
                Axis
            angle (float): Angle

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        R = self._axisangle_to_R(axis, angle)
        return SO3(R)
    
    def qtn(self, qtn):
        """[summary]
        Generate SO3 Class using an unit quaternion.

        Args:
            qtn (size=4, np.ndarray or list): 
                unit quaternion. 

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        R = self._quaternion_to_R(qtn)
        return SO3(R)
    
    def Rx(self, theta, unit='rad'):
        """[summary]
        Generate SO3 Class by rotating x axis.

        Args:
            theta (float): rotation angle along the x axis.
            unit (str, optional): [rad or deg]. Defaults to 'rad'.

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        return SO3(self._get_R_by_axis(theta, "x", unit))

    def Ry(self, theta, unit='rad'):
        """[summary]
        Generate SO3 Class by rotating y axis.

        Args:
            theta (float): rotation angle along the y axis.
            unit (str, optional): [rad or deg]. Defaults to 'rad'.

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        return SO3(self._get_R_by_axis(theta, "y", unit))

    def Rz(self, theta, unit='rad'):
        """[summary]
        Generate SO3 Class by rotating z axis.

        Args:
            theta (float): rotation angle along the z axis.
            unit (str, optional): [rad or deg]. Defaults to 'rad'.

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        return SO3(self._get_R_by_axis(theta, "z", unit))
    
    # --Conversion--   
    def inv(self):
        """[summary]
        Inverse of SO3. It's just transpose.

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        return SO3(self.R.T)
 
    def to_qtn(self):
        """[summary]
        Convert SO3 to unit quaternion. 

        Returns:
            qtn [size=4 np.ndarray]: unit quaternion array.
        """
        return self._R_to_quaternion(self.R)
    
    def to_axisangle(self):
        """[summary]
        Convert SO3 to AxisAngle. 

        Returns:
            axis [size=3 np.ndarray]: array.
            angle (float): angle.
        """
        return self._R_to_axisangle(self.R)

    #--private functions--
    def _axisangle_to_R(self, axis, angle):
        """[summary]
        Convert AxisAngle to Rotation Matrix.
        (Exponential of SO3)
        Args:
            axis (size3 np.ndarray or list): unit vector of rotation axis
            angle (float): rotation angle

        Returns:
            R [3x3 np.ndarray]: rotation matrix
        """
        axis = self._check_vector(3, axis)

        if angle == 0.:
            return np.eye(3)
        axis_norm = np.linalg.norm(axis)
        if axis_norm != 1.:
            axis = axis/axis_norm
        
        theta = angle
        omega_hat = axis
        return np.eye(3) \
           + np.sin(theta)*self._skew(omega_hat) \
           + (1-np.cos(theta))*self._skew(omega_hat).dot(self._skew(omega_hat))

    def _R_to_axisangle(self, R):
        """[summary]
        Convert rotation matrix to unit quaternion
        (logarithm of SO3)

        Args:
            R (3x3 np.ndarray): Rotation matrix.

        Returns:
            axis [size 3 np.ndarray]: Unit axis.
            angle (float): rotation angle
        """
        R = self._check_matrix(3, R)

        if np.allclose(R, np.eye(3)):
            # no rotation
            axis = np.array([1., 0., 0.])
            angle = 0.
        elif np.trace(R) == -1:
            # angle is 180 degrees
            angle = np.pi
            if R[0,0] != -1:
                axis = 1/np.sqrt(2*(1+R[0,0])) * np.array([1+R[0,0], R[1,0], R[2,0]])
            elif R[1,1] != -1:
                axis = 1/np.sqrt(2*(1+R[1,1])) * np.array([R[0,1], 1+R[1,1], R[2,1]])
            else:
                axis = 1/np.sqrt(2*(1+R[2,2])) * np.array([R[0,2], R[1,2], 1+R[2,2]])
        else:
            angle = np.arccos(1/2*(np.trace(R)-1))
            axis = 1/(2*np.sin(angle))*np.array([R[2,1]-R[1,2], R[0,2]-R[2,0], R[1,0]-R[0,1]])
        return axis, angle
        
    def _quaternion_to_R(self, qtn):
        """[summary]
        Convert Quaternion to Rotation matrix 
        Args:
            qtn (size3 np.ndarray): quaternion

        Returns:
            [3x3 np.ndarray]: Rotation Matrix
        """
        qtn = self._check_vector(4, qtn)
        qtn_norm = np.linalg.norm(qtn)
        if qtn_norm != 1:
            w, x, y, z = qtn/qtn_norm
        else:
            w, x, y, z = qtn
        return np.array([[1-2*y**2-2*z**2, 2*x*y-2*w*z, 2*x*z+2*w*y],
                         [2*x*y+2*w*z, 1-2*x**2-2*z**2, 2*y*z-2*w*x],
                         [2*x*z-2*w*y, 2*y*z+2*w*x, 1-2*x**2-2*y**2]])
    
    def _R_to_quaternion(self, R):
        """[summary]
        Convert rotation matrix to unit quaternion

        Args:
            R (3x3 np.ndarray): Rotation matrix.

        Returns:
            qtn [size 3 np.ndarray]: Unit quaternion.
        """
        tr = np.trace(R)
        if tr > 0.:
            s = np.sqrt(tr+1.)*2
            w = 0.25*s
            x = (R[2,1] - R[1,2])/s
            y = (R[0,2] - R[2,0])/s
            z = (R[1,0] - R[0,1])/s
        elif (R[0,0]>R[1,1]) & (R[0,0] > R[2,2]):
            s= np.sqrt(1. + R[0,0] - R[1,1] - R[2,2]) * 2
            w = (R[2,1] - R[1,2]) / s
            x = 0.25 * s
            y = (R[0,1] + R[1,0]) / s; 
            z = (R[0,2] + R[2,0]) / s; 
        elif R[1,1] > R[2,2]:
            s= np.sqrt(1. + R[1,1] - R[0,0] - R[2,2]) * 2 
            w = (R[0,2] - R[2,0]) / s
            x = (R[0,1] + R[1,0]) / s
            y = 0.25 * s
            z = (R[1,2] + R[2,1]) / s
        else:
            s = np.sqrt(1. + R[2,2] - R[1,1] - R[0,0]) * 2
            w = (R[1,0] - R[0,1]) / s
            x = (R[0,2] + R[2,0]) / s
            y = (R[1,2] + R[2,1]) / s
            z = 0.25 * s
        return np.array([w,x,y,z])

    def _get_R_by_axis(self, theta, axis, unit="rad"):
        """[summary]
        private function to make Rotation matrix
        using a specified axis and an angle.

        Args:
            theta (float): rotation angle
            axis (str): ["x", "y", "z"] axis name string
            unit (str, optional): ["rad", "deg"]. Defaults to "rad".

        Raises:
            Exception: if unit is neither "rad" nor "deg".

        Returns:
            R [3x3 np.ndarray]: Rotation Matrix
        """
        if unit == "rad":
            pass
        elif unit == "deg":
            theta = theta/180.*np.pi
        else: 
            raise Exception("Wrong unit!")
        
        if axis == "x":
            omega = np.array([1,0,0]) 
        elif axis == "y":
            omega = np.array([0,1,0]) 
        elif axis == "z":
            omega = np.array([0,0,1])
        return self._axisangle_to_R(omega, theta)

    #--SO3 operators--
    def __matmul__(self, X):
        return self.dot(X)

    def dot(self, X):
        if type(X) is np.ndarray:
            return self.R.dot(X)
        else:
            return SO3(self.R.dot(X.R))

    def __repr__(self):
        return "SO3 Class\n"+self.R.__repr__()


class SE3(SO3):
    def __init__(self, *args):
        """[summary]

        Raises:
            Exception: [description]
        """
        try:
            if len(args) == 0:
                self.R = np.eye(3)
                self.p = np.zeros(3)
            elif len(args) == 1:
                if (type(args[0]) is np.ndarray) \
                    & (args[0].shape == (4,4)):
                    self.R, self.p = \
                            self._T_to_Rp(args[0])   
                else:
                    raise Exception()
            elif len(args) == 2:
                if (args[0].shape == (3,3))   \
                        & (len(args[1] == 3)):
                    self.R = args[0]
                    self.p = args[1]
            else:
                raise Exception()
        except:
            raise Exception("SE3 Initialize Error!")

    # --Construction of Transformation Matrix--
    def random(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        R, _ = np.linalg.qr(np.random.randn(3,3))
        p = np.random.uniform(-1,1, size=3)
        return SE3(R, p)
    
    def twistangle(self, twist, angle):
        """[summary]
        Construction of Transformation Matrix
        using twist and angle

        Args:
            twist (size 3 np.ndarray): unit twist axis
            angle (float): angle

        Returns:
            [SE3]: SE3 Class
        """
        R, p = self._twistangle_to_Rp(twist, angle)
        return SE3(R,p)

    def trans(self, p):
        """[summary]
        Construction of Transformation Matrix
        using translation vector
        Args:
            p (size 3 np.ndarray): translation

        Returns:
            [SE3]: SE3 Class
        """
        R = np.eye(3)
        p = self._check_vector(3,p)
        return SE3(R,p)

    def qtn_trans(self, qtn, trans):
        """[summary]
        Construction of Transformation Matrix
        using quaternion, translation vector
        Args:
            qtn (size 4 np.ndarray) : quaternion
            trans (size 3 np.ndarray): translation

        Returns:
            [SE3]: SE3 Class
        """
        R = self._quaternion_to_R(qtn)
        p = self._check_vector(3,trans)
        return SE3(R,p)

    def Rx(self, theta, unit='rad'):
        """[summary]
        Generate SE3 Class by rotating x axis.

        Args:
            theta (float): rotation angle along the x axis.
            unit (str, optional): [rad or deg]. Defaults to 'rad'.

        Returns:
            T [SE3]: SE3 Translation Class
        """
        R = self._get_R_by_axis(theta, "x", unit)
        p = np.zeros(3)
        return SE3(R, p)

    def Ry(self, theta, unit='rad'):
        """[summary]
        Generate SE3 Class by rotating y axis.

        Args:
            theta (float): rotation angle along the y axis.
            unit (str, optional): [rad or deg]. Defaults to 'rad'.

        Returns:
            T [SE3]: SE3 Translation Class
        """
        R = self._get_R_by_axis(theta, "y", unit)
        p = np.zeros(3)
        return SE3(R, p)

    def Rz(self, theta, unit='rad'):
        """[summary]
        Generate SE3 Class by rotating z axis.

        Args:
            theta (float): rotation angle along the z axis.
            unit (str, optional): [rad or deg]. Defaults to 'rad'.

        Returns:
            T [SE3]: SE3 Translation Class
        """
        R = self._get_R_by_axis(theta, "z", unit)
        p = np.zeros(3)
        return SE3(R, p)

    # --Conversion--  
    def inv(self):
        """[summary]
        Inversion of Transformation Matrix

        Returns:
            [SE3]: Inversion of SE3 Class
        """
        R, p = self.R, self.p
        return SE3(R.T, -R.T.dot(p))
    
    def to_qtn_trans(self):
        """[summary]
        Convert SE3 to unit quaternion and translation vector. 

        Returns:
            qtn [size=4 np.ndarray]: unit quaternion array.
            trans [size=3 np.ndarray] translation vector array.
        """
        qtn = self._R_to_quaternion(self.R)
        trans = self.p
        return qtn, trans
    
    def to_twistangle(self):
        """[summary]
        Convert SE3 to TwistAngle. 

        Returns:
            Twist [size=6 np.ndarray]: array.
            angle (float): angle.
        """
        return self._Rp_to_twistangle(self.R, self.p)


    def _Rp_to_T(self, R, p):
        """[summary]

        Args:
            R ([type]): [description]
            p ([type]): [description]

        Returns:
            [type]: [description]
        """
        R = self._check_matrix(3, R)
        p = self._check_vector(3, p)
        return np.block([[R,       p[:,None]],
                         [np.zeros((1,3)), 1]])
                         
    def _T_to_Rp(self, T):
        """[summary]

        Args:
            T ([type]): [description]

        Returns:
            [type]: [description]
        """
        T = self._check_matrix(4, T)
        R, p = T[:3,:3], T[:3,3]
        return R, p

    def _Rp_to_twistangle(self, R, p):
        """[summary]

        Args:
            T ([type]): [description]

        Returns:
            [type]: [description]
        """
        T = self._Rp_to_T(R, p)
        if np.allclose(T, np.eye(4)):
            tw = np.array([1.,0.,0.,0.,0.,0])
            theta = 0
        elif np.allclose(R, np.eye(3)):
            # pure translation
            omega = np.zeros(3)
            theta = np.linalg.norm(p)
            v = p/theta
            tw = np.hstack([omega, v])
        else:
            omega, theta = self._R_to_axisangle(R)
            Ginv = 1/theta*np.eye(3) \
                - 1/2*self._skew(omega) \
                + (1/theta-1/2/np.tan(theta/2))*self._skew(omega).dot(self._skew(omega))
            v = Ginv.dot(p)
            tw = np.hstack([omega, v])
        return tw, theta
    
    def _twistangle_to_Rp(self, tw, angle):
        """[summary]
        Convert Twist-angle to R, p

        Args:
            tw (size:6 np.ndarray): twist
            angle (float) : angle

        Returns:
            [type]: [description]
        """
        tw = self._check_vector(6, tw)
        omega, v = tw[:3], tw[3:]
        theta = angle
        G = np.eye(3)*theta \
            + (1-np.cos(theta))*self._skew(omega) \
            + (theta-np.sin(theta))*self._skew(omega).dot(self._skew(omega))

        if np.linalg.norm(omega) == 0:
            R = np.eye(3)
            p = v*theta
        else:
            R = self._axisangle_to_R(omega, theta)
            p = G.dot(v)
        
        return R, p

    #--SO3 operators--
    def T(self):
        """[summary]
        Get Transformation Matrix

        Returns:
            T [4x4 np.ndarray]: transformation Matrix
        """
        return self._Rp_to_T(self.R, self.p)
    
    def __matmul__(self, X):
        return self.dot(X)

    def dot(self, X):
        if type(X) is np.ndarray:
            if X.shape == (3,):
                X = np.hstack([X,1])[:,None]
            elif X.shape == (3,1):
                X = np.vstack([X,1])
            return self.T().dot(X)
        else:
            return SE3(self.dot(X.T()))

    def __repr__(self):
        return "SE3 Class\n"+self.T().__repr__()

if __name__ == "__main__":
    a = SE3().random()
    print(a)

