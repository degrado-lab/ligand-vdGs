from setuptools import setup, find_packages

# # Read the contents of your README file
# def read_readme():
#     with open("README.md", "r") as f:
#         return f.read()

setup(
    name='ligand_vdgs',
    version='0.1.0',  # Update this to the version of your package
    description='A package for generating and using vdGs for docking.',  # Short description of your package
    # long_description=read_readme(),  # Long description read from README.md

    #author='Your Name',  # Replace with your name
    #author_email='your.email@example.com',  # Replace with your email address
    #url='https://github.com/yourusername/smart_vdms',  # Replace with the URL of your package repository
    packages=find_packages(),  # Automatically find packages in the current directory
    # classifiers=[
    #     'Development Status :: 3 - Alpha',
    #     'Intended Audience :: Developers',
    #     'License :: OSI Approved :: MIT License',
    #     'Programming Language :: Python :: 3',
    #     'Programming Language :: Python :: 3.8',  # Specify Python versions you support
    #     'Programming Language :: Python :: 3.9',
    #     'Programming Language :: Python :: 3.10',
    # ],
    python_requires='>=3.7',  # Minimum Python version required
    install_requires=[
        "prody==2.0.1",
    ],
    #entry_points={
    #    'console_scripts': [
            # If you have command-line scripts, list them here, e.g.,
            #'svdm-trim-database=ligand_vdgs.scripts.trim_database:main',
    #    ],
    #},
    #include_package_data=True,  # Include non-code files specified in MANIFEST.in
    #zip_safe=False,
)
