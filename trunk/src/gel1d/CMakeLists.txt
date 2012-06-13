PROJECT(gel1d)

ADD_EXECUTABLE(gel1d main.cpp config.lua)
TARGET_LINK_YOCTO(gel1d aqueous lua lua++ swamp mk)

