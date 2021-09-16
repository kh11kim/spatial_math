import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

def frame_plot(R, p, frame=None, dims=[0,10]*3, axeslength=1, axeswidth=5):
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
    def __init__(self):
        self.R = np.eye(3)
        self.p = np.zeros(3)

    def _check_vector(self, vector_len, v):
        if type(v) is list:
            v = np.array(v)
        assert v.shape == (vector_len,)
        return v

    def _check_matrix(self, matrix_len, m):
        assert m.shape == (matrix_len, matrix_len)
        return m
    
    def _skew(self, x):
        x = self._check_vector(3, x)
        return np.array([[0    , -x[2], x[1] ],
                        [x[2] , 0    , -x[0]],
                        [-x[1], x[0] , 0    ]])

    def plot(self, **kwargs):
        frame_plot(self.R, self.p, **kwargs)

class SO3(base):
    def __init__(self, *args):
        if len(args)==0:
            self.R = np.eye(3)
        elif len(args)==1:
            if (type(args[0]) is np.ndarray)    \
                    & (args[0].shape == (3,3)):
                self.R = args[0]
            elif (len(args[0]) == 3):
                print(args[0])
                self.R = self._omega_to_R(args[0])
        else:
            raise Exception("SO3 Initialize Error!")
        self.p = np.zeros(3)
    
    def random(self):
        omega = np.random.uniform(-1,1, size=3)
        return SO3(self._omega_to_R(omega))
    
    def Rx(self, theta, unit='rad'):
        return SO3(self._get_R_by_axis(theta, "x", unit))

    def Ry(self, theta, unit='rad'):
        return SO3(self._get_R_by_axis(theta, "y", unit))

    def Rz(self, theta, unit='rad'):
        return SO3(self._get_R_by_axis(theta, "z", unit))
    
    def _omega_to_R(self, omega):
        omega = self._check_vector(3, omega)

        if np.linalg.norm(omega) == 0:
            return np.eye(3)
        theta = np.linalg.norm(omega)
        omega_hat = omega/theta
        return np.eye(3) \
           + np.sin(theta)*self._skew(omega_hat) \
           + (1-np.cos(theta))*self._skew(omega_hat)@self._skew(omega_hat)

    def _R_to_omega(self, R):
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
        R = self._get_R_by_axis(theta, "x", unit)
        p = np.zeros(3)
        return SE3(R, p)

    def Ry(self, theta, unit='rad'):
        R = self._get_R_by_axis(theta, "y", unit)
        p = np.zeros(3)
        return SE3(R, p)

    def Rz(self, theta, unit='rad'):
        R = self._get_R_by_axis(theta, "z", unit)
        p = np.zeros(3)
        return SE3(R, p)

    def random(self):
        omega = np.random.uniform(-1,1, size=3)
        R = self._omega_to_R(omega)
        p = np.random.uniform(-1,1, size=3)
        return SE3(R, p)

    def inv(self):
        R, p = self.R, self.p
        return SE3(R.T, -R.T@p)
    
    def _Rp_to_T(self, R, p):
        R = self._check_matrix(3, R)
        p = self._check_vector(3, p)
        return np.block([[R,       p[:,None]],
                         [np.zeros((1,3)), 1]])
                         
    def _T_to_Rp(self, T):
        T = self._check_matrix(4, T)
        R, p = T[:3,:3], T[:3,3]
        return R, p

    def _T_to_twist(self, T):
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

