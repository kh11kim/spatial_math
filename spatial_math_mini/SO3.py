import numpy as np
from spatial_math_mini.base import Base

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
        if len(args) == 0:
            self._qtn = np.array([1,0,0,0])
            
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
    def qtn(qtn):
        """[summary]
        Generate SO3 Class using an unit quaternion.

        Args:
            qtn (size=4, np.ndarray or list): 
                unit quaternion. 

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        qtn = Base._check_vector(4,qtn)
        return SO3(qtn)
    
    @staticmethod
    def Rx(theta, unit='rad'):
        """[summary]
        Generate SO3 Class by rotating x axis.

        Args:
            theta (float): rotation angle along the x axis.
            unit (str, optional): [rad or deg]. Defaults to 'rad'.

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        return SO3(SO3._get_quaternion_by_axis(theta, "x", unit))

    @staticmethod
    def Ry(theta, unit='rad'):
        """[summary]
        Generate SO3 Class by rotating y axis.

        Args:
            theta (float): rotation angle along the y axis.
            unit (str, optional): [rad or deg]. Defaults to 'rad'.

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        return SO3(SO3._get_quaternion_by_axis(theta, "x", unit))

    @staticmethod
    def Rz(theta, unit='rad'):
        """[summary]
        Generate SO3 Class by rotating z axis.

        Args:
            theta (float): rotation angle along the z axis.
            unit (str, optional): [rad or deg]. Defaults to 'rad'.

        Returns:
            R [SO3]: SO3 Rotation Class
        """
        return SO3(SO3._get_quaternion_by_axis(theta, "x", unit))
    
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

    

    #--SO3 operators--
    def __matmul__(self, X):
        return self.dot(X)

    def dot(self, X):
        if type(X) is np.ndarray:
            if len(X) == 3:
                return self.R @ X
                
        elif type(X) is SO3:
            w1, v1 = self._qtn[0], self._qtn[1:]
            w2, v2 = X._qtn[0], X._qtn[1:]
            qtn = np.array([w1+w2-v1@v2, *(w1*v2 + w2*v1 + np.cross(v1, v2))])
            return SO3(qtn)

    def __repr__(self):
        return "SO3 Class\n"+self.R.__repr__()