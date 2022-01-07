# spatial_math_mini

### Description
Python Module for spatial vector algebra with minimum dependency (numpy, matplotlib).  
SO3, SE3 class and its arithmatic, visualization is implemented.   
  

### Install
```shell
pip install spatial_math_mini
```

### Usage
```python
from spatial_math_mini import SE3, SO3, Viz

# make SE3/SO3 easy
T0 = SE3(np.eye(4))
T1 = SE3.Rx(180, 'deg')
T2 = SE3.axisangle(axis, angle)

# calculate arithmatic using @(matmul), inv()
T = T1.inv() @ T2 @ T0
R, t = T.R, T.t
p = T @ np.array([1, 2, 3])

# conversion to matrix/axisangle/qtn
R = SO3.random()
qtn = R.to_qtn()
axis, angle = R.to_axisangle()
qtn, t = T.to_qtn_trans()
twist, angle = T.to_twist_angle()

# visualize results
v = Viz(axeslength=0.1, axeswidth=1, dims=[-0.5,0.5]*3)
v.plot(T0, name, color="r")
v.plot(T1, name)
v.clear()
```

### TODO
make utility function module (interpolation...)
