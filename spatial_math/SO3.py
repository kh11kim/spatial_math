import numpy as np
import matplotlib.pyplot as plt
from spatial_math.base import Base
        
class SO3(Base):
    """[summary]
    SO3 class
    """
    def __init__(self, *args):
        """[summary]
        Initialize SO3 Class
        Raises:
            Exception: [description]
        """
        self._qtn = np.empty(4)

        if len(args) == 0:
            self._qtn = np.array([1,0,0,0])
            #self.R = np.eye(3) # no input
        elif (len(args) == 1) \
                & (type(args[0]) is np.ndarray):
            if args[0].shape == (3,3):
                self.R = args[0] #input is 3x3 matrix
            elif args[0].shape == (4,):
                self._qtn = args[0] / np.linalg.norm(args[0])
            else:
                Exception("SO3 Initialize Error! : input is not R or qtn")
        else:
            raise Exception("SO3 Initialize Error! : the number of input should be one")
        
        # this is just for compatibility for plotting.
        #self.p = np.zeros(3) 
    
    # -- properties --
    @property
    def R(self):
        return self._quaternion_to_R(self._qtn)
    
    @R.setter
    def R(self, R):
        self._qtn = self._R_to_quaternion(R)
    
    def qtn(self):
        return self._qtn / np.linalg.norm(self._qtn)

    # --Construction of Rotation Matrix--
    @staticmethod
    def random():
        """[summary]
        Random Generation of SO3
        using QR decomposition

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        R, _ = np.linalg.qr(np.random.randn(3,3))
        return SO3(R)

    @staticmethod
    def axisangle(axis, angle):
        """[summary]
        Generate SO3 Class using an AxisAngle.

        Args:
            axis (size=3, np.ndarray or list): 
                Axis
            angle (float): Angle

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        qtn = SO3._axisangle_to_quaternion(axis, angle)
        return SO3(qtn)
    
    @staticmethod
    def qtn(self, qtn):
        """[summary]
        Generate SO3 Class using an unit quaternion.

        Args:
            qtn (size=4, np.ndarray or list): 
                unit quaternion. 

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        return SO3(qtn)
    
    @staticmethod
    def Rx(self, theta, unit='rad'):
        """[summary]
        Generate SO3 Class by rotating x axis.

        Args:
            theta (float): rotation angle along the x axis.
            unit (str, optional): [rad or deg]. Defaults to 'rad'.

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        return SO3(Base._get_quaternion_by_axis(theta, "x", unit))

    @staticmethod
    def Ry(self, theta, unit='rad'):
        """[summary]
        Generate SO3 Class by rotating y axis.

        Args:
            theta (float): rotation angle along the y axis.
            unit (str, optional): [rad or deg]. Defaults to 'rad'.

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        return SO3(Base._get_quaternion_by_axis(theta, "x", unit))

    @staticmethod
    def Rz(self, theta, unit='rad'):
        """[summary]
        Generate SO3 Class by rotating z axis.

        Args:
            theta (float): rotation angle along the z axis.
            unit (str, optional): [rad or deg]. Defaults to 'rad'.

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        return SO3(Base._get_quaternion_by_axis(theta, "x", unit))
    
    # --Conversion--   
    def inv(self):
        """[summary]
        Inverse of SO3. It's just transpose.

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        w, x, y, z = self._qtn
        return SO3(np.array([w, -x, -y, -z]))
 
    def to_qtn(self):
        """[summary]
        Convert SO3 to unit quaternion. 

        Returns:
            qtn [size=4 np.ndarray]: unit quaternion array.
        """
        return self._qtn
    
    def to_axisangle(self):
        """[summary]
        Convert SO3 to AxisAngle. 

        Returns:
            axis [size=3 np.ndarray]: array.
            angle (float): angle.
        """
        return self._quaternion_to_axisangle(self._qtn)

    #--private functions--
    @staticmethod
    def _get_quaternion_by_axis(self, theta, axis, unit="rad"):
        """[summary]
        private function to make quaternion
        using a specified axis and an angle.

        Args:
            theta (float): rotation angle
            axis (str): ["x", "y", "z"] axis name string
            unit (str, optional): ["rad", "deg"]. Defaults to "rad".

        Raises:
            Exception: if unit is neither "rad" nor "deg".

        Returns:
            qtn [np.ndarray(4)]: quaternion
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
        return Base._axisangle_to_quaternion(omega, theta)

    #--SO3 operators--
    def __matmul__(self, X):
        return self.dot(X)

    def dot(self, X):
        if type(X) is np.ndarray:
            if len(X) == 3:
                return self.R @ X
        if type(X) is SO3:
            qtn = np.multiply(self._qtn, X._qtn)
            return SO3(qtn)
            
        # if type(X) is np.ndarray:
        #     return self.R.dot(X)
        # else:
        #     return SO3(self.R.dot(X.R))

    def __repr__(self):
        return "SO3 Class\n"+self.R.__repr__()