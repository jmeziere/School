from setuptools import setup, find_packages
import os

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            if not filename.endswith('.pyc'):
                paths.append(os.path.join('..', path, filename))
    return paths

extra_files = package_files('MooseTutorial')
setup(
        name='MOOSETutorial',
        version='0.1dev',
        packages=find_packages(),
        package_data ={'': extra_files},
        install_requires=['kivy>=1.11','os','sys','subprocess','collection','image','setuptools'],
        )
