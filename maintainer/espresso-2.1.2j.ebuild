# Copyright 2010,2011,2012,2013,2014,2016 The ESPResSo project
# Copyright 1999-2009 Gentoo Foundation
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.
#
# $Header$

EAPI="2"

inherit autotools savedconfig

DESCRIPTION="Extensible Simulation Package for Research on Soft matter"
HOMEPAGE="http://espressomd.org"
SRC_URI="http://espressowiki.mpip-mainz.mpg.de/wiki/uploads/4/43/Espresso-2.1.2j.tar.gz"

LICENSE="GPL-2"
SLOT="0"
KEYWORDS="~x86"
IUSE="X doc examples fftw mpi packages test tk"

DEPEND="dev-lang/tcl
	X? ( x11-libs/libX11 )
	doc? ( app-doc/doxygen
		virtual/tex-base
		virtual/latex-base )
	fftw? ( sci-libs/fftw:3.0 )
	mpi? ( virtual/mpi )
	tk? ( >=dev-lang/tk-8.4.18-r1 )"

RDEPEND="${DEPEND}"

src_prepare() {
	AT_M4DIR="config" eautoreconf
	restore_config myconfig.h
}

src_configure() {
	#disable processor-optimization, we have make.conf
	econf \
		--disable-processor-optimization \
		$(use_with fftw) \
		$(use_with mpi) \
		$(use_with tk) \
		$(use_with X x)
}

src_compile() {
	emake || die "emake failed"
	use doc && emake doc || die "emake doc failed"
}

src_install() {
	cd "${S}"
	emake DESTDIR="${D}" install || die "Installing failed"

	dodoc INSTALL README RELEASE_NOTES

	insinto /usr/share/${PN}
	doins myconfig-sample.h

	if [ -f myconfig.h ]; then
		save_config myconfig.h
	else
		save_config config/myconfig.h
	fi

	if use doc; then
		newdoc doc/ug/ug.pdf user_guide.pdf
		dohtml -r doc/dg/html/*
		newdoc doc/tutorials/tut2/tut2.pdf tutorial.pdf
	fi

	if use examples; then
		insinto /usr/share/${PN}/examples
		doins samples/*
		#the testsuite are also good examples
		rm testsuite/Makefile* testsuite/test.sh.in
		insinto /usr/share/${PN}/testsuite
		doins testsuite/*
	fi

	if use packages; then
		insinto /usr/share/${PN}/packages
		doins -r packages/*
	fi

	echo "ESPRESSO_SOURCE=/usr/bin" > "${T}/80${PN}"
	echo "ESPRESSO_SCRIPTS=/usr/share/espresso/scripts" >> "${T}/80${PN}"
	doenvd "${T}/80${PN}"

	cd "${D}"
	#remove Espresso_wrapper
	rm -f usr/bin/Espresso
	#install Espresso directly
	newbin usr/libexec/Espresso_bin Espresso
	rm -f usr/libexec/Espresso_bin
}

pkg_postinst() {
	env-update && source /etc/profile
	elog
	elog Please read and cite:
	elog ESPResSo, Comput. Phys. Commun. 174\(9\) ,704, 2006.
	elog http://dx.doi.org/10.1016/j.cpc.2005.10.005
	elog
	elog If you need more features change
	elog /etc/portage/savedconfig/${CATEGORY}/${PF}
	elog and reemerge with USE=savedconfig
	elog
	elog For a full feature list see: 
	elog /usr/share/${PN}/myconfig-sample.h 
	elog
}
