from setuptools import setup

setup(
    name='KiLCA_interface',
    version='1.0.0',    
    description='Python interface for KiLCA',
    url='',
    author='Markus Markl',
    author_email='markl@tugraz.at',
    license='',
    packages=['KiLCA_interface'],
    install_requires=['numpy',
                      'matplotlib',
                      'scipy'                     
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)


