# Copyright 1999-2009 Gentoo Foundation
# Distributed under the terms of the GNU General Public License v2
# $Header$

EAPI="1"

LIBTOOLIZE="true"

inherit autotools savedconfig

DESCRIPTION="Extensible Simulation Package for Research on Soft matter"
HOMEPAGE="http://www.espresso.mpg.de"
SRC_URI="http://espressowiki.mpip-mainz.mpg.de/wiki/uploads/4/43/Espresso-2.1.2j.tar.gz"

LICENSE="GPL-2"
SLOT="0"
KEYWORDS="x86"
IUSE="X doc examples fftw mpi test packages"

DEPEND="dev-lang/tcl
	X? ( x11-libs/libX11 
		tk? ( dev-lang/tk ) )
	doc? ( app-doc/doxygen
		virtual/tex-base
		virtual/latex-base )
	fftw? ( >=sci-libs/fftw-3.0.1 )
	mpi? ( virtual/mpi )"

RDEPEND="${DEPEND}"

src_unpack() {
	unpack ${A}
	cd "${S}"

	AT_M4DIR="config" eautoreconf
}

src_compile() {
	restore_config myconfig.h

	#disable tk bug #225999, add tk to IUSE when fixed

	econf \
		--enable-processor-optimization \
		$(use_with fftw) \
		$(use_with mpi) \
		$(use_with tk) \
		$(use_with X x)

	emake || die "emake failed"
	use doc && emake doc || die "emake doc failed"
}

src_install() {
	cd "${S}"
	emake DESTDIR="${D}" install || die "Installing failed"

	dodoc INSTALL README RELEASE_NOTES

	insinto /usr/share/${PN}
	doins myconfig-sample.h

	save_config config/myconfig.h

	if use doc; then
		dodoc doc/ug/ug.pdf
		dohtml -r doc/dg/html/*
	fi

	if use examples; then
		insinto /usr/share/${PN}/examples
		doins samples/*
		#the testsuite are also good examples
		rm testsuite/Makefile*
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
}

pkg_postinst() {
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
