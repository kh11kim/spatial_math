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
    is_3d_plot = False
    if len(plt.get_fignums()) != 0:
        fig = plt.gcf()
        if len(fig.axes) != 0:
            ax = plt.gca()
            if ax.__class__.__name__ == "Axes3DSubplot":
                is_3d_plot = True
    
    if not is_3d_plot:
        fig = plt.figure()
        ax = plt.axes(projection='3d')    
    
    xaxis = np.vstack([p, p+R[:3,0]]).T * axeslength
    yaxis = np.vstack([p, p+R[:3,1]]).T * axeslength
    zaxis = np.vstack([p, p+R[:3,2]]).T * axeslength

    ax.plot(*xaxis, 'r', linewidth=axeswidth)
    ax.plot(*yaxis, 'g', linewidth=axeswidth)
    ax.plot(*zaxis, 'b', linewidth=axeswidth)

    ax.set_xlim3d(dims[:2])
    ax.set_ylim3d(dims[2:4])
    ax.set_zlim3d(dims[4:])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    frame_name_offset = [-0.5, -0.5, -0.5]
    frame_name_pos = R@np.array(frame_name_offset).T+p
    if frame is not None:
        ax.text(*frame_name_pos, "{{{}}}".format(frame))



class base:
    """[summary]
    """
    def __init__(self):
        self.R = np.eye(3)
        self.p = np.zeros(3)

    def _check_vector(self, vector_len, v):
        """[summary]

        Args:
            vector_len ([type]): [description]
            v ([type]): [description]

        Returns:
            [type]: [description]
        """
        if type(v) is list:
            v = np.array(v)
        assert v.shape == (vector_len,)
        return v

    def _check_matrix(self, matrix_len, m):
        """[summary]

        Args:
            matrix_len ([type]): [description]
            m ([type]): [description]

        Returns:
            [type]: [description]
        """
        assert m.shape == (matrix_len, matrix_len)
        return m
    
    def _skew(self, x):
        """[summary]

        Args:
            x ([type]): [description]

        Returns:
            [type]: [description]
        """
        x = self._check_vector(3, x)
        return np.array([[0    , -x[2], x[1] ],
                        [x[2] , 0    , -x[0]],
                        [-x[1], x[0] , 0    ]])

    def plot(self, **kwargs):
        """[summary]
        """
        frame_plot(self.R, self.p, **kwargs)

class SO3(base):
    """[summary]
    """
    def __init__(self, *args):
        """[summary]

        Raises:
            Exception: [description]
        """
        if len(args)==0:
            self.R = np.eye(3)
        elif len(args)==1:
            if (type(args[0]) is np.ndarray)    \
                    & (args[0].shape == (3,3)):
                self.R = args[0]
            elif (len(args[0]) == 3):
                self.R = self._omega_to_R(args[0])
            elif (len(args[0]) == 3):
                self.R = self._quaternion_to_R(args[0])
        else:
            raise Exception("SO3 Initialize Error!")
        self.p = np.zeros(3)
    
    def random(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        omega = np.random.uniform(-1,1, size=3)
        return SO3(self._omega_to_R(omega))
    
    def Rx(self, theta, unit='rad'):
        """[summary]

        Args:
            theta ([type]): [description]
            unit (str, optional): [description]. Defaults to 'rad'.

        Returns:
            [type]: [description]
        """
        return SO3(self._get_R_by_axis(theta, "x", unit))

    def Ry(self, theta, unit='rad'):
        """[summary]

        Args:
            theta ([type]): [description]
            unit (str, optional): [description]. Defaults to 'rad'.

        Returns:
            [type]: [description]
        """
        return SO3(self._get_R_by_axis(theta, "y", unit))

    def Rz(self, theta, unit='rad'):
        """[summary]

        Args:
            theta ([type]): [description]
            unit (str, optional): [description]. Defaults to 'rad'.

        Returns:
            [type]: [description]
        """
        return SO3(self._get_R_by_axis(theta, "z", unit))

    def to_quaternion(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        return self._R_to_quaternion(self.R)

    def _quaternion_to_R(self, q):
        q = self._check_vector(4, q)
        w, x, y, z = q
        return np.array([[1-2*y**2-2*z**2, 2*x*y-2*w*z, 2*x*z+2*w*y],
                         [2*x*y+2*w*z, 1-2*x**2-2*z**2, 2*y*z-2*w*x],
                         [2*x*z-2*w*y, 2*y*z+2*w*x, 1-2*x**2-2*y**2]])
    
    def _R_to_quaternion(self, R):
        omega = self._R_to_omega(R)
        theta = np.linalg.norm(omega)
        if theta == 0:
            pass
        omega_x, omega_y, omega_z = omega/theta
        c = np.cos(theta/2)
        s = np.sin(theta/2)
        return np.array([c, s*omega_x, s*omega_y, s*omega_z])



    
    def _omega_to_R(self, omega):
        """[summary]

        Args:
            omega ([type]): [description]

        Returns:
            [type]: [description]
        """
        omega = self._check_vector(3, omega)

        if np.linalg.norm(omega) == 0:
            return np.eye(3)
        theta = np.linalg.norm(omega)
        omega_hat = omega/theta
        return np.eye(3) \
           + np.sin(theta)*self._skew(omega_hat) \
           + (1-np.cos(theta))*self._skew(omega_hat)@self._skew(omega_hat)

    def _R_to_omega(self, R):
        """[summary]

        Args:
            R ([type]): [description]

        Raises:
            Exception: [description]

        Returns:
            [type]: [description]
        """
        R = self._check_matrix(3, R)

        if np.array_equal(R, np.eye(3)):
            raise Exception("̂ω is undefined since R=I")
        elif np.trace(R) == -1:
            theta = np.pi
            if R[0,0] != -1:
                omega_hat = 1/np.sqrt(2*(1+R[0,0])) * np.array([1+R[0,0], R[1,0], R[2,0]])
            elif R[1,1] != -1:
                omega_hat = 1/np.sqrt(2*(1+R[1,1])) * np.array([R[0,1], 1+R[1,1], R[2,1]])
            else:
                omega_hat = 1/np.sqrt(2*(1+R[2,2])) * np.array([R[0,2], R[1,2], 1+R[2,2]])
        else:
            theta = np.arccos(1/2*(np.trace(R)-1))
            omega_hat = 1/(2*np.sin(theta))*np.array([R[2,1]-R[1,2], R[0,2]-R[2,0], R[1,0]-R[0,1]])
        return theta*omega_hat

    def _get_R_by_axis(self, theta, axis, unit="rad"):
        """[summary]

        Args:
            theta ([type]): [description]
            axis ([type]): [description]
            unit (str, optional): [description]. Defaults to "rad".

        Raises:
            Exception: [description]

        Returns:
            [type]: [description]
        """
        if unit == "rad":
            pass
        elif unit == "deg":
            theta = theta/180*np.pi
        else: 
            raise Exception("Wrong unit!")
        
        if axis == "x":
            omega = theta * np.array([1,0,0]) 
        elif axis == "y":
            omega = theta * np.array([0,1,0]) 
        elif axis == "z":
            omega = theta * np.array([0,0,1])
        return self._omega_to_R(omega)

    def __matmul__(self, X):
        return self.R @ X

    def __repr__(self):
        return self.R.__repr__()


class SE3(SO3):
    def __init__(self, *args):
        """[summary]

        Raises:
            Exception: [description]
        """
        if len(args) == 0:
            self.R = np.eye(3)
            self.p = np.zeros(3)
        elif len(args) == 1:
            if len(args[0]) == 3:
                self.R = np.eye(3)
                self.p = self._check_vector(3,args[0])
            elif (type(args[0]) is np.ndarray) \
                    & (args[0].shape == (4,4)):
                self.R, self.p = self._T_to_Rp(args[0])
        elif len(args) == 2:
            if (args[0].shape == (3,3))   \
                    & (len(args[1] == 3)):
                self.R = args[0]
                self.p = args[1]
        else:
            raise Exception("SE3 Initialize Error!")
        self.T = self._Rp_to_T(self.R, self.p)

    def Rx(self, theta, unit='rad'):
        """[summary]

        Args:
            theta ([type]): [description]
            unit (str, optional): [description]. Defaults to 'rad'.

        Returns:
            [type]: [description]
        """
        R = self._get_R_by_axis(theta, "x", unit)
        p = np.zeros(3)
        return SE3(R, p)

    def Ry(self, theta, unit='rad'):
        """[summary]

        Args:
            theta ([type]): [description]
            unit (str, optional): [description]. Defaults to 'rad'.

        Returns:
            [type]: [description]
        """
        R = self._get_R_by_axis(theta, "y", unit)
        p = np.zeros(3)
        return SE3(R, p)

    def Rz(self, theta, unit='rad'):
        """[summary]

        Args:
            theta ([type]): [description]
            unit (str, optional): [description]. Defaults to 'rad'.

        Returns:
            [type]: [description]
        """
        R = self._get_R_by_axis(theta, "z", unit)
        p = np.zeros(3)
        return SE3(R, p)

    def random(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        omega = np.random.uniform(-1,1, size=3)
        R = self._omega_to_R(omega)
        p = np.random.uniform(-1,1, size=3)
        return SE3(R, p)

    def inv(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        R, p = self.R, self.p
        return SE3(R.T, -R.T@p)
    
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

    def _T_to_twist(self, T):
        """[summary]

        Args:
            T ([type]): [description]

        Returns:
            [type]: [description]
        """
        R, p = self._T_to_Rp(T)
        if np.array_equal(R, np.eye(3)):
            omega = np.zeros(3)
            theta = np.linalg.norm(p)
            v = p/theta
        else:
            omega = self._R_to_omega(R)
            theta = np.linalg.norm(omega)
            Ginv = 1/theta*np.eye(3) \
                - 1/2*self._skew(omega) \
                + (1/theta-1/2/np.tan(theta/2))*self._skew(omega)@self._skew(omega)
            v = Ginv@p
        return np.hstack([omega, v])
    
    def _twist_to_T(self, tw):
        """[summary]

        Args:
            tw ([type]): [description]

        Returns:
            [type]: [description]
        """
        tw = self._check_vector(6, tw)
        omega, v = tw[:3], tw[3:]
        theta = np.linalg.norm(omega)
        G = np.eye(3)*theta \
            + (1-np.cos(theta))*self._skew(omega) \
            + (theta-np.sin(theta))*self._skew(omega)@self._skew(omega)
        return self._Rp_to_T(self._omega_to_R(omega), G@v)

    def __matmul__(self, X):
        return self.T@X

    def __repr__(self):
        return self.T.__repr__()
    
if __name__ == "__main__":
    a = SE3().random()
    print(a)

