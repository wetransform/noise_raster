# -*- coding: utf-8 -*-

# Resource object code
#
# Created by: The Resource Compiler for PyQt5 (Qt v5.15.2)
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore

qt_resource_data = b"\
\x00\x00\x01\x3a\
\x89\
\x50\x4e\x47\x0d\x0a\x1a\x0a\x00\x00\x00\x0d\x49\x48\x44\x52\x00\
\x00\x00\x1c\x00\x00\x00\x1c\x08\x03\x00\x00\x00\x45\xd3\x2f\xa6\
\x00\x00\x00\x01\x73\x52\x47\x42\x00\xae\xce\x1c\xe9\x00\x00\x00\
\x04\x67\x41\x4d\x41\x00\x00\xb1\x8f\x0b\xfc\x61\x05\x00\x00\x00\
\x27\x50\x4c\x54\x45\x00\x00\x00\x07\x08\x08\x18\x1c\x1b\xb8\xd6\
\xd1\x40\x40\x40\x08\x02\x02\x1b\x09\x08\xcd\x46\x3e\x09\x07\x05\
\x20\x1a\x11\xf3\xc6\x83\x50\x50\x50\x00\x00\x00\x9e\xb3\x83\x2f\
\x00\x00\x00\x0d\x74\x52\x4e\x53\xff\xff\xff\xff\xff\xff\xff\xff\
\xff\xff\xff\xff\x00\x3d\xe8\x22\x86\x00\x00\x00\x09\x70\x48\x59\
\x73\x00\x00\x0e\xc3\x00\x00\x0e\xc3\x01\xc7\x6f\xa8\x64\x00\x00\
\x00\x83\x49\x44\x41\x54\x38\x4f\x95\x92\xd1\x0e\xc0\x10\x0c\x45\
\x31\x8c\xcd\xff\x7f\xef\xa6\xbb\x28\x95\x89\xf3\xa2\x7a\xb8\x91\
\x86\x4a\x3f\xa8\xa4\x18\x68\x16\x5e\xa9\x0d\xd0\xf5\x64\x93\xe6\
\x00\xa6\x35\xcb\x2a\xa4\x7d\xaf\x5a\xaa\x66\x37\x6b\xee\x6e\x6c\
\x63\x94\xf4\x58\x20\xa5\xf3\xc0\x4d\xa4\x3f\x81\x27\xc9\x87\x20\
\x24\x43\xc8\x9e\x1d\x89\x44\x42\xca\x10\x41\x98\xc8\x78\x81\xb8\
\x94\xc3\x10\x7a\xc9\x10\xb2\x63\x2a\xdb\x2a\x65\x4e\xfc\x0a\x29\
\xf3\x37\xb9\xa9\x5a\xc5\xf2\x21\xd0\x2b\x4b\x6a\x96\x0c\x34\x0b\
\xe3\x9e\x91\xd2\x03\x28\x49\x16\x58\x70\x0e\x1a\x5a\x00\x00\x00\
\x00\x49\x45\x4e\x44\xae\x42\x60\x82\
"

qt_resource_name = b"\
\x00\x07\
\x07\x3b\xe0\xb3\
\x00\x70\
\x00\x6c\x00\x75\x00\x67\x00\x69\x00\x6e\x00\x73\
\x00\x0c\
\x08\x24\xb9\xa2\
\x00\x6e\
\x00\x6f\x00\x69\x00\x73\x00\x65\x00\x5f\x00\x72\x00\x61\x00\x73\x00\x74\x00\x65\x00\x72\
\x00\x08\
\x0a\x61\x5a\xa7\
\x00\x69\
\x00\x63\x00\x6f\x00\x6e\x00\x2e\x00\x70\x00\x6e\x00\x67\
"

qt_resource_struct_v1 = b"\
\x00\x00\x00\x00\x00\x02\x00\x00\x00\x01\x00\x00\x00\x01\
\x00\x00\x00\x00\x00\x02\x00\x00\x00\x01\x00\x00\x00\x02\
\x00\x00\x00\x14\x00\x02\x00\x00\x00\x01\x00\x00\x00\x03\
\x00\x00\x00\x32\x00\x00\x00\x00\x00\x01\x00\x00\x00\x00\
"

qt_resource_struct_v2 = b"\
\x00\x00\x00\x00\x00\x02\x00\x00\x00\x01\x00\x00\x00\x01\
\x00\x00\x00\x00\x00\x00\x00\x00\
\x00\x00\x00\x00\x00\x02\x00\x00\x00\x01\x00\x00\x00\x02\
\x00\x00\x00\x00\x00\x00\x00\x00\
\x00\x00\x00\x14\x00\x02\x00\x00\x00\x01\x00\x00\x00\x03\
\x00\x00\x00\x00\x00\x00\x00\x00\
\x00\x00\x00\x32\x00\x00\x00\x00\x00\x01\x00\x00\x00\x00\
\x00\x00\x01\x7f\x8d\x62\xcf\x0a\
"

qt_version = [int(v) for v in QtCore.qVersion().split('.')]
if qt_version < [5, 8, 0]:
    rcc_version = 1
    qt_resource_struct = qt_resource_struct_v1
else:
    rcc_version = 2
    qt_resource_struct = qt_resource_struct_v2

def qInitResources():
    QtCore.qRegisterResourceData(rcc_version, qt_resource_struct, qt_resource_name, qt_resource_data)

def qCleanupResources():
    QtCore.qUnregisterResourceData(rcc_version, qt_resource_struct, qt_resource_name, qt_resource_data)

qInitResources()
