from setuptools import setup

setup(
    name='KiLCA_QB',
    version='1.0.0',    
    description='Python interface for KiLCA and QL-Balance',
    url='',
    author='Markus Markl',
    author_email='markl@tugraz.at',
    license='',
    py_modules=['pyABC', 'fieldpy', 'AnLoBiCr', 'neo2_for_Er', 'KQ_processor', 'tMHD_current', 'device_config', 'misalign_field', 'susc_functions', 'postproc_class', 'KiLCA_interface', 'profile_processor', 'balance_interface', 'analytical_local_bif_criterion'],
    install_requires=['numpy',
                      'matplotlib',
                      'scipy',
                      'h5py',
                      'f90nml'           
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        'Operating System :: OS Independent',        
        'Programming Language :: Python :: 3',
    ],
    python_requires='>=3.6',
)


