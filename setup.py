from setuptools import setup, find_packages
 
setup(
    name                = 'spatial_math_mini',
    version             = '0.1.3.5',
    description         = 'Python Module for spatial math with minimum dependency',
    author              = 'Kanghyun Kim',
    author_email        = 'kh11kim@kaist.ac.kr',
    url                 = 'https://github.com/kh11kim/spatial_math_mini',
    install_requires    =  ["numpy", "matplotlib>=3.5"],
    packages            = find_packages(exclude = []),
    keywords            = ['spatial_math_mini'],
    python_requires     = '>=3',
    package_data        = {},
    zip_safe            = False,
    classifiers         = [
        'Programming Language :: Python'
    ],
)
"""build process
1. python setup.py bdist_wheel
2. twine upload dist/spatial_math_mini-x.x.x-py3-none-any.whl
3. pip install spatial_math_mini -U
"""