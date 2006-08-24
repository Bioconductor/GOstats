PKG=GOstats
TOP=../../..
RUNIT_DIR=$(shell pwd)

all: test

testonly:
	R --slave < runalltests.R

test:
	cd ${TOP}; \
	R CMD INSTALL ${PKG} && \
	cd ${RUNIT_DIR}
	R --slave < runalltests.R

